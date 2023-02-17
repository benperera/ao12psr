
import logging

import psrchive
import numpy as np
import matplotlib.pyplot as plt



def plot_raw(data, freq, tbin, rfi_zap=None, outfile=None, pshow=False):
    """

    Plot the raw data

    Args:
        data (np.ndarray): Data array to plot. Array should include full resolution, not scrunched data.
        freq (np.ndarray): Center frequencies of all channels
        tbin (float): Sampling time in sec
        rfi_zap (np.ndarray): maksed data from IQRM algorithm.
        outfile (str): Output png plot file name. If it is not set, outfile is "rawdata.png"
        pshow (bool): Show the plot

    Returns:
        None 

    """
    if outfile is None:
        outfile = 'rawdata.png'
    logging.info(f'Creating raw data plot {outfile} ....')
    tsub = tbin*data.shape[1] #time of a subint
    nrows = data.shape[0]
    npols = data.shape[2]
    plt.figure(0,figsize=(17,8))
    grid = plt.GridSpec(2, 4, wspace=0.3, hspace=0.3)   
    spectra_aa = np.mean(data[:,:,0,:,0], axis=1)
    pwrt_aa = np.mean(spectra_aa, axis=1)
    bp_aa = np.mean(spectra_aa, axis=0)
    if npols != 1 :
        spectra_bb = np.mean(data[:,:,1,:,0], axis=1)
        pwrt_bb = np.mean(spectra_bb, axis=1)
        bp_bb = np.mean(spectra_bb, axis=0)
    if npols == 1:
        plt.subplot(grid[:,0]).set_title('Spectra: Total I')
    else:
        plt.subplot(grid[:,0]).set_title('Spectra: Pol A')
    plt.ylabel('Time (s)')
    plt.xlabel('Frequency (MHz)')
    plt.imshow(spectra_aa,origin='lower', aspect='auto',extent=[freq[0],freq[len(freq)-1],
                    0,spectra_aa.shape[0]*tsub],vmin=0, vmax = 2.*np.median(bp_aa))
    plt.colorbar()
    #if rfi_zap.shape != 0:
    if rfi_zap is not None:
        for h in range(rfi_zap.shape[0]):
            yval = rfi_zap[h,0,:,0] + h*tsub
            plt.plot(np.linspace(freq[0],freq[len(freq)-1],len(freq)),yval,linestyle='',marker='o',markersize=.3,color='r')
        plt.gcf().text(0.77,0.35,'Cleaned by IQRM', fontsize=10)
    if npols != 1:
        plt.subplot(grid[:,1]).set_title('Spectra: Pol B')
        plt.xlabel('Frequency (MHz)')
        plt.imshow(spectra_bb,origin='lower', aspect='auto',extent=[freq[0],freq[len(freq)-1],
                        0,spectra_bb.shape[0]*tsub],vmin=0, vmax = 2.*np.median(bp_bb))
        #if rfi_zap.shape != 0:
        if rfi_zap is not None:
            for h in range(rfi_zap.shape[0]):
                yval = rfi_zap[h,1,:,0] + h*tsub
                plt.plot(np.linspace(freq[0],freq[len(freq)-1],len(freq)),yval,linestyle='',marker='o',markersize=.3,color='r')
        plt.colorbar()

    plt.subplot(grid[0,2]).set_title('Power vs Time')
    plt.xlabel('Time (s)')
    plt.xlim(0,nrows*tsub)
    if npols == 1:
        plt.plot(np.arange(nrows)*tsub, pwrt_aa, label='Total I',linestyle=':',color='g')
    else:
        plt.plot(np.arange(nrows)*tsub, pwrt_aa, label='Pol A',linestyle=':',color='g')
        plt.plot(np.arange(nrows)*tsub, pwrt_bb, label='Pol B',linestyle=':',color='orange')
    plt.legend(loc='best')

    if npols != 1:
        plt.subplot(grid[1,2])
        plt.xlabel('Time (s)')
        plt.xlim(0,nrows*tsub) 
        plt.plot(np.arange(nrows)*tsub, (pwrt_aa+pwrt_bb), label='Tot I',linestyle=':',color='b')
        plt.legend(loc='best')

    plt.subplot(grid[0,3]).set_title('Bandpass')
    plt.xlabel('Frequency (MHz)')
    plt.xlim(freq[0], freq[len(freq)-1])
    if npols == 1:
        plt.plot(freq, bp_aa, label='Total I',color='g')
    else:
        plt.plot(freq, bp_aa, label='Pol A',color='g')
        plt.plot(freq, bp_bb, label='Pol B',color='orange')       
    plt.legend(loc='best')

    plt.gcf().text(0.77,0.25,'Arecibo Observatory', fontsize=10)
    plt.gcf().text(0.77,0.23,'12 m Telescope', fontsize=10)


    plt.savefig(outfile,bbox_inches="tight",dpi=200)
    if pshow is True:
        plt.show()


