#!/usr/bin/env python3

import logging
import argparse
import numpy as np
import multiprocessing as mp
from ao12psr import utils, writer, plotter
from ao12psr.utils import YourArgparseFormatter




if __name__ == "__main__":
    __spec__ = None
    #__spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    parser = argparse.ArgumentParser(
        prog="psr12run.py",
        description="AO 12-m telescope Pulsar data processing pipeline",
        formatter_class=YourArgparseFormatter,
    )
    parser.add_argument("-f", "--file", help="psrfits file or the fits file sequence as a wild card", nargs="+")
    parser.add_argument(
        "-band",
        "--band",
        help="Band to be processed ('sband' or 'xband')",
        type=str,
        required=False,
        default='sband',
    )
    parser.add_argument(
        "-nproc",
        "--nproc",
        help="Number of CPUs for fits file reading",
        type=int,
        required=False,
        default=8,
    )
    parser.add_argument(
        "-no_iqrm",
        "--no_iqrm",
        help="Disable RFI excision via IQRM algorithm",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-p",
        "--p",
        help="Make the plots of raw data",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-w",
        "--w",
        help="Write 8-bit converted new psrfits file (newraw.fits)",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-wo",
        "--wout",
        help="Output psrfits file name",
        type=str,
        required=False,
        default=None,
    )
    parser.add_argument(
        "-b",
        "--band_corr",
        help="Correct for the bandshape",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-tot",
        "--tot_int",
        help="Total internsity - sum pol channels",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-pshow",
        "--pshow",
        help="Show the raw data plot",
        required=False,
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-po",
        "--pout",
        help="Output plot file name. This should be .png file",
        type=str,
        required=False,
        default=None,
    )
    #parser.add_argument("-v", "--verbose", help="Be verbose", action="store_true")
    values = parser.parse_args()

    logger = logging.getLogger()
    logging_format = '%(asctime)s - %(levelname)s - %(message)s'
    logging.basicConfig(level=logging.INFO,format=logging_format, datefmt='%Y/%m/%d %I:%M:%S %p')

    logging.info("Input Arguments:-")
    for arg, value in sorted(vars(values).items()):
        logging.info("%s: %r", arg, value)
    
    if not values.file:
        raise ValueError("[-f FILE] Input FITS file is missing!")
    else:
        logging.info(f'Input files {values.file}')

        if mp.cpu_count() < values.nproc:
                values.nproc = mp.cpu_count()
        if values.band == 'sband':
            freq, data, tbin0 = utils.read_sband(values.file)
        if values.band == 'xband':
            freq, data, tbin0 = utils.read_xband(values.file)
        if values.tot_int:
            data = utils.tot_int(data)
        if values.no_iqrm is False:
            data, rfi_zap = utils.iqrm_rfi(data,10) #10 s chuncks
        if values.band_corr:
            data = utils.bandshape_corr(data)
        if values.p:
            if values.no_iqrm is False:
                plotter.plot_raw(data, freq, tbin0, rfi_zap, outfile=values.pout, pshow=values.pshow)
            else:
                plotter.plot_raw(data, freq, tbin0, outfile=values.pout, pshow=values.pshow)
        if values.w:
            writer.write_fits(data, freq, values.file[0], outfile=values.wout)


