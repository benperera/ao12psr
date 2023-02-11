
import logging
import argparse

import numpy as np
import astropy.io.fits as pyfits
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
        fn (str): input fits file
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
    return [fn, freqs, data, df, tbin]


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


