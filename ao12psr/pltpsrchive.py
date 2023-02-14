
import logging
import psrchive

import numpy as np
import matplotlib.pyplot as plt
from ao12psr import utils


def plot_psr(fn, outfile=None, pshow=False):
    """

    Plotting pulsar data using the psrchive processed files

    Args:
        fn (str): psrchive file name (unscrunched data file)
        outfile (str): Output png plot file name. If it is not set, outfile is "psrdata.png"
        pshow (bool): Show the plot
	
    Returns:
	None

    """
    if outfile is None:
        outfile = 'psrdata.png'
    plt.figure(0,figsize=(14,8))
    grid = plt.GridSpec(2, 3, wspace=0.3, hspace=0.3)
    obs = psrchive.Archive_load(fn)
    obs.pscrunch()
    obs.dedisperse()
    obs.remove_baseline()
    nbins = obs.get_nbin()
    freq = obs.get_frequencies()
    data0 = obs.get_data()
    weights = obs.get_weights()
    for ww in range(nbins):
        data0[:,0,:,ww] *= weights
    data0 = np.flip(data0, axis=2)
    freq = np.flip(freq)
    freq_phase = np.mean(data0,axis=0)
    time_phase = np.mean(data0,axis=2)
    prof = np.mean(time_phase, axis=0)

    snr = utils.get_snr(prof[0,:])

    plt.subplot(grid[0,0])
    plt.xlim(0,1)
    plt.xlabel('Pulse Phase')
    plt.ylabel('Amp (arb. units)')
    plt.plot(np.arange(nbins)/(nbins*1.), prof[0,:])

    plt.subplot(grid[0,1])
    plt.xlabel('Pulse Phase')
    plt.ylabel('Frequency (MHz)')
    plt.imshow(freq_phase[0,:,:], origin='lower', aspect='auto', 
                extent=[0,1,freq[0],freq[-1]], cmap='hot')
    plt.colorbar()

    plt.subplot(grid[1,1])
    plt.xlabel('Pulse Phase')
    plt.ylabel('Time (s)')
    plt.imshow(time_phase[:,0,:], origin='lower', aspect='auto', 
                extent=[0,1,0,time_phase.shape[0]*10], cmap='hot')
    plt.colorbar()

    shift =  int(nbins/2. - np.argmax(prof))
    data0 = np.roll(data0, shift, axis=3)
    off_pulse = np.concatenate((data0[:,:,:,0:int(nbins/3.)],
                data0[:,:,:,int(2.*nbins/3.):nbins-1]),axis=3)
    freq_time = np.std(off_pulse, axis=3)

    plt.subplot(grid[:,2])
    plt.xlabel('Frequency (MHz)')
    plt.ylabel('Time (s)')
    plt.imshow(freq_time[:,0,:], origin='lower', aspect='auto', 
                extent=[freq[0],freq[-1],0,time_phase.shape[0]*10],cmap='hot')
    plt.colorbar()

    plt.gcf().text(0.17,0.35,'Profile S/N = '+str(snr), fontsize=10)
    plt.gcf().text(0.17,0.25,'Arecibo Observatory', fontsize=10)
    plt.gcf().text(0.17,0.23,'12 m Telescope', fontsize=10)

    plt.savefig(outfile,bbox_inches="tight",dpi=200)
    if pshow is True:
        plt.show()
