# AO12PSR
Arecibo Observatory 12-m telescope Pulsar processing pipeline.

This the pipeline uses to process the pulsar data obtained daily with the 12-m telescope of the Arecibo Observatory. The telescope is currently equipped with S-band (2.3 GHz) and X-band (8.6 GHz) room-temperature receivers. The pulsar observations are mainly conducted using the 2.3 GHz receiver, and the data are recorded using the Mocks spectrometer.


## Install
```
  git clone git@github.com:benperera/ao12psr.git
  cd ao12psr
  pip install -e .
```
The installation will add `psr12run.py` to your `PYTHONPATH`.


## Usage
The script `psr12run.py` has various options to handle the raw data files such as read in the original 12-bit pstfits data files, clean the data for radio frequency interference, plot the data, output 8-bit psrfits files, correct for the band shape, sum polarization channels, etc.
