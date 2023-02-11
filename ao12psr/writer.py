
import logging

import numpy as np
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.io import fits



def write_fits(data, freq, fn, outfile=None):
    """

    Write a new psrfits file in 8-bit format

    Args:
        data (np.ndarray): Data array to write
        freq (float): Center frequencies of all channels
        fn (str): name of the original 16-bit fits file
        outfile (str): name of the output fits file

    Returns:
        None

    """
    if outfile is None:
        outfile = 'newrawdata.fits'
    logging.info(f'Writing 8-bit converted fits file {outfile}  ....')
    nrows = data.shape[0]
    npols = data.shape[2]
    freq_tbl = np.zeros((nrows,len(freq)))

    #convert the 16 bit data into 8 bit data
    data = data.astype(np.float32)
    for fsub in range(nrows):
        data[fsub,:,:,:,:] -= np.mean(data[fsub,:,:,:,:])
        data[fsub,:,:,:,:] /= np.std(data[fsub,:,:,:,:])
        data[fsub,:,:,:,:] *= np.sqrt(2 ** 8)
        data[fsub,:,:,:,:] += 2 ** 7
        np.round(data[fsub,:,:,:,:], out=data[fsub,:,:,:,:])
        np.clip(data[fsub,:,:,:,:], 0, 2 ** 8 - 1, out=data[fsub,:,:,:,:])
    data = data.astype('uint8')

    #flip the band before writing
    freq = np.flip(freq)
    data = np.flip(data, axis=3)

    for i in range(nrows):
        freq_tbl[i,:] = freq
    hdr0 = getheader(fn, 0)    
    hdr0['HDRVER'] = '3.4'
    hdr0['OBSFREQ'] = (freq[0] + freq[-1])/2.
    hdr0['OBSBW'] = freq[-1] - freq[0]
    hdr0['OBSNCHAN'] = len(freq)
    hdr0['TRK_MODE'] = 'TRACK'
    hdr0['TELESCOP'] = 'ARECIBO_12'
    hdr0['BACKEND'] = 'PUPPI'
    data3, hdr3 = getdata(fn, 3, header=True)
    hdr3['NCHAN'] = len(freq)
    hdr3['NPOL'] = npols
    if npols == 1:
        hdr3['POL_TYPE'] = 'AA+BB'
    else:
        hdr3['POL_TYPE'] = 'AABB'

    primary_hdu = fits.PrimaryHDU(header=hdr0)
    col1 = data3['TSUBINT']
    col2 = data3['OFFS_SUB']
    col3 = data3['LST_SUB']
    col4 = data3['RA_SUB']
    col5 = data3['DEC_SUB']
    col6 = data3['GLON_SUB']
    col7 = data3['GLAT_SUB']
    col8 = data3['FD_ANG']
    col9 = data3['POS_ANG']
    col10 = data3['PAR_ANG']
    col11 = data3['TEL_AZ']
    col12 = data3['TEL_ZEN']
    col13 = data3['DAT_FREQ']
    col14 = data3['DAT_WTS']
    col15 = data3['DAT_OFFS']
    col16 = data3['DAT_SCL']
    col17 = data3['DATA']
    col18 = data3['STAT']
    dat_wts = np.full((nrows,len(freq)), 1)
    dat_offs = np.full((nrows,len(freq)*npols), 0)
    dat_scl = np.full((nrows,len(freq)*npols), 1)
    c1 = fits.Column(name='TSUBINT', array=col1, unit='s', format='1D')
    c2 = fits.Column(name='OFFS_SUB', array=col2, unit='s', format='1D')
    c3 = fits.Column(name='LST_SUB', array=col3, unit='s', format='1D')
    c4 = fits.Column(name='RA_SUB', array=col4, unit='deg', format='1D')
    c5 = fits.Column(name='DEC_SUB', array=col5, unit='deg', format='1D')
    c6 = fits.Column(name='GLON_SUB', array=col6, unit='deg', format='1D')
    c7 = fits.Column(name='GLAT_SUB', array=col7, unit='deg', format='1D')
    c8 = fits.Column(name='FD_ANG', array=col8, unit='deg', format='1E')
    c9 = fits.Column(name='POS_ANG', array=col9, unit='deg', format='1E')
    c10 = fits.Column(name='PAR_ANG', array=col10, unit='deg', format='1E')
    c11 = fits.Column(name='TEL_AZ', array=col11, unit='deg', format='1E')
    c12 = fits.Column(name='TEL_ZEN', array=col12, unit='deg', format='1E')
    c13 = fits.Column(name='DAT_FREQ', array=freq_tbl, unit='MHz', format=str(hdr3['NCHAN'])+'E')
    c14 = fits.Column(name='DAT_WTS', array=dat_wts, format=str(hdr3['NCHAN'])+'E')
    c15 = fits.Column(name='DAT_OFFS', array=dat_offs, format=str(hdr3['NCHAN']*hdr3['NPOL'])+'E')
    c16 = fits.Column(name='DAT_SCL', array=dat_scl, format=str(hdr3['NCHAN']*hdr3['NPOL'])+'E')
    c17 = fits.Column(name='DATA', array=data, unit='Jy', dim='(1,'+str(hdr3['NCHAN'])+','+
        str(hdr3['NPOL'])+','+str(hdr3['NSBLK'])+')', format=str(hdr3['NCHAN']*hdr3['NPOL']*hdr3['NSBLK'])+'B')
    cc = fits.ColDefs([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17])
    table_hdu = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15,
        c16, c17], name='SUBINT')
    hdul = fits.HDUList([primary_hdu, table_hdu])
    hh = hdul[1].header
    hh.set('INT_TYPE', 'TIME', before='EXTNAME')
    hh.set('INT_UNIT', 'SEC', before='EXTNAME')
    hh.set('SCALE', 'FluxDen', before='EXTNAME')
    hh.set('NPOL', hdr3['NPOL'], before='EXTNAME')
    hh.set('POL_TYPE', hdr3['POL_TYPE'], before='EXTNAME')
    hh.set('TBIN', hdr3['TBIN'], before='EXTNAME')
    hh.set('NBIN', 1, before='EXTNAME')
    hh.set('NBIN_PRD', 0, before='EXTNAME')
    hh.set('PHS_OFFS', 0., before='EXTNAME')
    hh.set('NBITS', 8, before='EXTNAME')
    hh.set('NSUBOFFS', 0, before='EXTNAME')
    hh.set('NCHAN', hdr3['NCHAN'], before='EXTNAME')
    hh.set('CHAN_BW', -1.*np.abs(hdr3['CHAN_BW']), before='EXTNAME')
    hh.set('NCHNOFFS', 0, before='EXTNAME')
    hh.set('NSBLK', hdr3['NSBLK'], before='EXTNAME')
    hdul.writeto(outfile)
    
