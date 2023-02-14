# AO12PSR
Arecibo Observatory 12-m telescope Pulsar processing pipeline.

This pipeline is used to process the pulsar data obtained with the 12-m telescope of the Arecibo Observatory. The telescope is currently equipped with S-band (2.3 GHz) and X-band (8.6 GHz) room-temperature receivers, and the data are recorded using the Mock Spectrometer. The current pulsar observations are mainly conducted at 2.3 GHz frequency. 


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
