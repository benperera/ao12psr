# AO12PSR
Arecibo Observatory 12-m telescope Pulsar processing pipeline.

The 12-m telescope is currently equipped with S-band (2.3 GHz) and X-band (8.6 GHz) room-temperature receivers, and the astronomical data are recorded using the Mock Spectrometer. The S-band and X-band receivers have bandwidths of ~135 MHz and ~1 GHz, respectively. The Mock Spectrometer has seven boxes, and each box can be configured to operate with 172 MHz bandwidth. Therefore, the S-band receiver is configured to use a single Mock box, and the X-band receiver is configured to use all seven boxes to cover its ~1 GHz bandwidth. The Mock data are recorded with the PSRFITS format at 16-bit.



This pipeline is developed to process the pulsar data obtained with the 12-m. The pulsar observations are mainly conducted at S-band. The pipeline pre-processes the data by cleaning for radio frequency interference and downscaling to 8-bit as most pulsar data processing packages are not compatible with 16-bit format. It can produce raw data plots (including spectra over time, power fluctuation over time, and band shape plots for all polarization channels) and rewrite the downscaled data as PSRFITS files. These new data files can be processed directly using pulsar data processing packages. The X-band pulsar observations use seven spectrometer boxes, and the data are written as seven separate PSRFITS files, one for each sub-band box. The pipeline can combine the seven bands together to produce a total bandwidth of 1 GHz during the pre-processing. 



## Install
```
  git clone git@github.com:benperera/ao12psr.git
  cd ao12psr
  pip install -e .
```
The installation will add `psr12run.py` to your `PYTHONPATH`.


## Usage
The command line script `psr12run.py` has various options to handle the raw data files, such as reading in the original 12-bit psrfits format data files, cleaning the data for radio frequency interference, plotting the data, outputting 8-bit psrfits files, correcting for the band shape, summing polarization channels, etc. 

To read the raw data, make plots, and output the 8-bit-converted psrfits file, you can do:
```
  psr12run.py -f filename.fits -p -w
```

The tools in the pipeline can also be used manually as shown in [this example notebook](https://github.com/benperera/ao12psr/blob/main/example/plot_write_data.ipynb).
