# AO12PSR
Arecibo Observatory 12-m telescope Pulsar processing pipeline.

This the pipeline uses to process the pulsar data obtained daily with the 12-m telescope of the Arecibo Observatory. The telescope is currently equipped with S-band (2.3 GHz) and X-band (8.6 GHz) room-temperature receivers. The pulsar observations are mainly conducted using the 2.3 GHz receiver, and the data are recorded using the Mocks spectrometer.


## Install
```
  git clone git@github.com:benperera/ao12psr.git
  cd ao12psr
  pip install -e .
```
The installation will add `psr12run.py` to your `PYTHONPATH`. The 


## Usage
The command line script `psr12run.py` has various options to handle the raw data files, such as reading in the original 12-bit psrfits format data files, cleaning the data for radio frequency interference, plotting the data, outputting 8-bit psrfits files, correcting for the band shape, summing polarization channels, etc. 

Plot the raw data and write the 8-bit-converted psrfits file as follow:
```
  psr12run.py -f filename.fits -p -w
```

The tools in the pipeline can also be used manually as shown in [this example notebook](https://github.com/benperera/ao12psr/blob/main/example/plotting_data.ipynb).
