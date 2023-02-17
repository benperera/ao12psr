
import logging
import argparse, os

import numpy as np
import astropy.io.fits as pyfits
import multiprocessing as mp
from iqrm import iqrm_mask
from astropy.io.fits import getdata



class YourArgparseFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
    """
    Allows both Raw Text Formatting and Default Args
    """
    pass


def read_data(fn):
    """

    Read in the data of a psrfits file.

    Args:
        fn (str): Psrfits file name

    Returns:
        freqs (np.ndarray): center frequencies of channels
        data (np.ndarray): Data array
        df (float): channels bandwidth
        tbin (float): sample size of the data in sec

    """
    logging.info(f'Read in fits file {fn}....')
    fits = pyfits.open(fn, mode="readonly", memmap=True)
    npol = fits["SUBINT"].header["NPOL"]
    nchan = fits["SUBINT"].header["NCHAN"]
    tbin = fits["SUBINT"].header["TBIN"]
    obsfreq = fits["PRIMARY"].header["OBSFREQ"]
    df = fits["SUBINT"].header["CHAN_BW"]
    bw = fits["PRIMARY"].header["OBSBW"]
    dd, hdr = getdata(fn, 3, header=True)
    data = dd['DATA']#.astype(np.float32)
    fits.close()
    #seems not required to consider scale,offset, and weights
    #Needs to be checked it is required for future data!
    if df < 0.:
        data = np.flip(data, axis=3)
        df = np.abs(df)
    freqs = np.zeros(nchan)
    freq_str = obsfreq - (bw/2.) + df/2.
    for i in range(nchan):
        freqs[i] = freq_str + (df*i)
    return [freqs, data, df, tbin]


def iqrm_rfi(data, tchn):
    """
    Excise Radio Frequency Inteference using IQRM algorithm.

    Args:
        data (np.ndarray): Data array to clean. Structure is similar to fits file table structure.
        tchn (int): Data chunk length in sec to run iqrm

    Returns:
        np.ndarray, np.ndarray: Clean data array, masked elements of the cleaned dara array
        
    """
    logging.info(f'Running IQRM cleaning....')
    npols = data.shape[2]
    data0 = np.mean(data, axis=1)
    data0_zap = np.full(data0.shape, np.nan)
    chunks = int(data0.shape[0]/tchn*1.)
    if data0.shape[0] > chunks*tchn: chunks += 1
    for z in range(npols):
        for j in range(chunks):
            chn_str = tchn*j
            chn_stp = tchn*(j+1)
            if chn_stp > data0.shape[0]: chn_stp = data0.shape[0] - 1
            spec_std = data0[chn_str:chn_stp,z,:,0].std(axis=0)
            maski, votes = iqrm_mask(spec_std, radius=3, threshold = 8.0)
            mask = np.full(spec_std.shape, False, dtype=bool)
            for caster in sorted(votes.keys()):
                rec = votes[caster]
                rec0 = str(rec).replace('{','')
                rec0 = rec0.replace('}','')
                if rec0.find(',') != -1:
                    rec0 = rec0.split(',')
                    zap = []
                    zap.append(int(caster))
                    for ele in rec0:
                        zap.append(int(ele))
                    for i in range(len(zap)):
                        if mask[zap[i]] == False: mask[zap[i]] = True
                else:
                    zap = [int(caster),int(rec0)]
                    for i in range(len(zap)):
                        if mask[zap[i]] == False: mask[zap[i]] = True
            cc = 0
            block = 0
            for k in range(len(mask)):
                if mask[k] == True and cc == 0:
                    s1 = data0[chn_str:chn_stp,z,k-1,0]
                    ch_str = k-1
                    cc = 1
                if mask[k] == False and cc == 1:
                    s2 = data0[chn_str:chn_stp,z,k,0]
                    ch_stp = k
                    cc = 0
                    block = 1
                if block == 1:
                    for r in range(len(s1)):
                        chans = np.arange(ch_str+1, ch_stp)
                        pwr_intrp = np.interp(chans, [ch_str,ch_stp], [s1[r],s2[r]])
                        if len(chans) == 1:
                            data[chn_str+r,:,z,chans[0],:] = pwr_intrp[0]
                            data0_zap[chn_str+r,z,chans[0],:] = 1.
                        else:
                            for q in range(len(chans)):
                                data[chn_str+r,:,z,chans[q],:] = pwr_intrp[q]
                                data0_zap[chn_str+r,z,chans[q],:] = 1.
                        block = 0
    return data, data0_zap


def get_data_spec(data):
    """

    Specs of the data array.

    Args:
        data (np.ndarray): Data array to get specs

    Returns:
        int: number of rows, samples, polarizations, and frequency channels

    """
    nrows, nsamp, npols, nfreq, qq = data.shape
    return nrows, nsamp, npols, nfreq

    
def bandshape_corr(data):
    """

    Correct for the band shape.

    Args:
        data (np.ndarray): Data array for the band shape correction

    Returns:
        np.ndarray: Band shape corrected array   

    """
    logging.info(f'Correct for the band shape....')
    npols = get_data_spec(data)[2]
    nfreq = get_data_spec(data)[3]
    if npols == 1:
        bp = np.mean(np.mean(data[:,:,0,:,0], axis=1), axis=0)
        for i in range(nfreq):
            data[:,:,0,i,0] = data[:,:,0,i,0]/bp[i]
    else:
        bp_aa = np.mean(np.mean(data[:,:,0,:,0], axis=1), axis=0)
        bp_bb = np.mean(np.mean(data[:,:,1,:,0], axis=1), axis=0)
        for i in range(nfreq):
            data[:,:,0,i,0] = data[:,:,0,i,0]/bp_aa[i]
            data[:,:,1,i,0] = data[:,:,1,i,0]/bp_bb[i]
    return data


def tot_int(data):
    """
    
    Sum Polarization channels AA and BB to obtain Total Intensity

    Args:
        data (np.ndarray): Data array to get total intensity

    Returns:
        np.ndarray: Total intensity data array 

    """
    logging.info(f'Summ Pol channels and get Total intensity only....')
    data = np.sum(data, axis=2, keepdims=True)
    return data


def filenames(fn):
    """

    Collect the names of file sequence.

    Args:
        str: Name of the first data file

    Retunrs:
        list (str): Names of files of the sequence 

    """
    proj, date, psr, beam, seq, dd = fn.split('.')
    seqf = float(seq)
    #extract the info of raw fits files
    max_slice = 100  #maximum time slices
    files = []
    for i in range(max_slice):
        seq1 = str(int(seqf+i)).zfill(5)
        ff = proj+'.'+date+'.'+psr+'.'+beam+'.'+seq1+'.fits'
        check_file = os.path.isfile('./'+ff)
        if check_file:
            files = np.append(files, ff)
    if len(files) == 0:
        raise ValueError(f"Input file {fn} does not exist!")
    return files


def run_dspsr(fn, psr, subint_length=10, nbins=256, prefix='psrchive'):
    """

    Process the psrfits files using dspsr package

    Args:
        fn (str): psrfits file name
	psr (str): Pulsar name for the par file (par file should name as psrname.par)
	subint_length (int): Length of the sub-integration in sec 
	nbins (int): Number of bins required across the profile
	prefix (str): Prefix of the output archive file

    Returns:
	None

    """
    logging.info(f'Running DSPSR and creating the archive file {prefix}.ar')
    cmd = 'dspsr -E '+psr+'.par'+' -A -L '+str(subint_length)+' -b '+str(nbins)+' -O '+prefix+' '+fn
    os.system(cmd)


def get_snr(y, prop=False):
    """

    Calculate the Signal-to-noise of the pulse profile

    Args:
	y (float): Pulse profile data
	prop (bool): Estimate the S/N using Eq 7.1 in HBPA, otherwise Amp/off_rms

    Returns:
	snr (float): S/N of the profile	

    """
    nbin = len(y)
    shift = int(nbin/2. - np.argmax(y))
    y = np.roll(y, shift)
    off = np.concatenate((y[0:int(nbin/4.)],y[int(3.*nbin/4.):nbin-1]))
    y -= np.mean(off)
    weq = np.sum(y[int(nbin/4.):int(3.*nbin/4.)])/np.amax(y)
    off = np.concatenate((y[0:int(nbin/4.)],y[int(3.*nbin/4.):nbin-1]))
    if prop:
        snr = np.round(np.sum(y-np.mean(off))/np.std(off)/np.sqrt(weq), 2)
    else:
	#If the profile is noisy (e.g. J0437), weq becomes -ve. So decided to use 
	#snr = peak/off_std (May 27, 2022)
    	snr = np.amax(y)/np.std(off)
    return snr


def skip_chan(freqs1, freqs2, df1):
    """

    Frequency channels to skip in the overlap region of bands when combining them

    Args:
	freqs1 (np.ndarray): Frequencies of the lower band
	freqs2 (np.ndarray): Frequencies of the upper band
	df1 (float): Bandwidth of a channel

    Returns:
	upchanskip (float): Number of channels to skip in the lower band
	lowchanskip (float): Number of channels to skip in the upper band

    """
    upperfreqoflower = freqs1.max()
    lowerfreqofupper = freqs2.min()
    nextfromlower = upperfreqoflower + np.abs(df1)
    numchandiff = int(np.round((nextfromlower - lowerfreqofupper) / np.abs(df1)))
    chanskip = numchandiff if numchandiff > 0 else 0
    (upchanskip, lowchanskip) = (chanskip // 2, chanskip // 2 + 1)
    return upchanskip, lowchanskip


def read_sband(files, nproc=8):
    """

    Read S-band data of the 12-m telescope

    Args:
        files (str): List of file names of the 7band data
        nproc (int): Number of cpu cores to use in parallel

    Returns:
        freq (np.ndarray): Frequencies of the band
        data (np.ndarray): Data array
        tbin (float): Sample time of the data

    """
    if mp.cpu_count() < nproc:
        logging.info('Num of cpu cores < nproc: set nrpoc = cpu_cores')
        nproc = mp.cpu_count()
    pool = mp.Pool(ncpus)
    dd = pool.map(utils.read_data, values.file)
    freq = dd[0][0]
    tbin0 = dd[0][3]
    for i in range(len(values.file)):
        if i > 0:
            data = np.concatenate((data, dd[i][1]),axis=0)
        else:
            data = dd[i][1]
    return [freq, data, tbin]


def read_xband(files, nproc=8, pwr_scaling=False):
    """

    Combine all 7 bands of the Mocks data to get the full bandwidth

    Args:
        files (str): List of file names of the 7band data
        nproc (int): Number of cpu cores to use in parallel
        pwr_scaling (bool): Enable scaling the power of overlap regions
                            of adjecent bands

    Returns:
        freq (np.ndarray): Frequencies of the combined full band
        data (np.ndarray): Data array of the combined data
        tbin (float): Sample time of the data

    """
    if mp.cpu_count() < nproc:
        logging.info('Num of cpu cores < nproc: set nrpoc = cpu_cores')
        nproc = mp.cpu_count()
    pool = mp.Pool(nproc)
    dd = pool.map(read_data, files)
    df0 = dd[0][2]
    tbin = dd[0][3]
    for w in range(len(files)):
        globals()[f'freq{w}'] = dd[w][0]
        globals()[f'data{w}'] = dd[w][1]
    nrows, nsamp, npols, nfreq, tmp = data0.shape 

    #***combining all seven bands to get the full bandwidth***
    logging.info('Combining all 7 bands')
    freq_sband0 = freq0
    data_sband0 = data0
    for j in range(len(files)-1):
        freq_sband1 = globals()[f"freq{j+1}"]
        data_sband1 = globals()[f"data{j+1}"]
        uband_chanskip, lband_chanskip = skip_chan(freq_sband0, freq_sband1, df0)
        freq = np.concatenate((freq_sband0[:-lband_chanskip], freq_sband1[uband_chanskip:]))
        #print('Nchan = ', len(freq))
        if pwr_scaling is True:
            for i in range(npols):
                pwr_d0 = np.mean(np.mean(data_sband0[:,:,i,:,0], axis=1), axis=0)
                pwr_d1 = np.mean(np.mean(data_sband1[:,:,i,:,0], axis=1), axis=0)
                pwr_scl = np.mean(pwr_d0[-lband_chanskip:])/np.mean(pwr_d1[:uband_chanskip])
                data_sband1[:,:,i,:,:] = data_sband1[:,:,i,:,:]*pwr_scl
        datatmp = np.concatenate((data_sband0[:,:,:,:-lband_chanskip,:], 
                            data_sband1[:,:,:,uband_chanskip:,:]),axis=3)
        #print('\nSHAPE : ', datatmp.shape)
        freq_sband0 = freq
        data_sband0 = datatmp
    #padding extral freq channels if the No. of channels is not a 
    #multiple of 32. This condition is required for heimdall
    check_multi = len(freq) % 32.
    if check_multi != 0:
        chan_pad = int(32 - check_multi)
        for i in range(chan_pad):
            value = freq[len(freq)-1] + df0
            freq = np.append(freq,value)
    #print('Nchan pad = ', chan_pad)
    data_pad = np.full((nrows,nsamp,npols,chan_pad,tmp), 0.0001)    #padding zero messed up with 
                                                                    #band shape correction
    data = np.concatenate((datatmp,data_pad),axis=3)
    return [freq, data, tbin]

