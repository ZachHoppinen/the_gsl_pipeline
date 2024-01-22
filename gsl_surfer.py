#!/home/rdcrlzh1s/miniforge3/envs/gamma/bin/python

import os
import argparse
from pathlib import Path
import subprocess
import shlex
import ast
import json

# progress bar stuck to bottom
# https://stackoverflow.com/questions/76063620/stick-top-progress-bar-in-python-to-stdout
import enlighten

# pbar = manager.counter(total = len(raws), desc = 'GSL Raw to SLC')
# for raw in raws:
#   pbar.update()

import numpy as np
import pandas as pd

from typing import Union
import logging

from gs_proc import raw_gsl_proc, parse_param
from gs_funcs import execute, assert_one_fp, parse_positions
from gamma_class import GSLRawData

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

# set up progress bar managet
manager = enlighten.get_manager()
status = manager.status_bar(status_format=u'GSL_Surfer{fill}Step: {step}{fill}{elapsed}',
                                    color='bold_underline_bright_white_on_lightslategray',
                                    justify=enlighten.Justify.CENTER, step='GSL Raw to SLC',
                                    autorefresh=True, min_delta=0.5)

if __name__ == "__main__":

    # Parse incoming command line arguments
    parser = argparse.ArgumentParser(description='Convert GSL data to SLCs and generate all available interferograms.')
    parser.add_argument('dir', metavar='in-dir', type=Path,
                        help='The location directory containing raw data, nav data, and DEMs')
    
    parser.add_argument('-p', metavar='positions', type=Path, nargs = '?',
                        help='[Optional] The file containing start and end locations')
    
    parser.add_argument('-d', metavar='dem_fp', type=Path, nargs = '?',
                        help='[Optional] DEM filepath - can be tif or GAMMA')
    
    parser.add_argument('--default_dem', metavar='default_dem_fp', type=Path, nargs = '?',
                        help='[Optional] DEM default filepath - can be tif or GAMMA. Used if no dem is found')

    parser.add_argument('-n', metavar='nav_fp', type=Path, nargs = '?',
                    help='[Optional] NAV filepath - atlans. Will add 18 seconds if it does not have extension ".txt_plus18sec"')

    parser.add_argument('-w', default = False,
                        action='store_true', help='[Optional] Should the dem be converted to WGS84?')
    
    parser.add_argument('--error_log', metavar='error_log', type=Path, nargs = '?',
                        help='[Optional] Path to save collected errors to. Will be created.')
    
    parser.add_argument('-o', metavar='out_dir', type=Path, nargs = '?',
                        help='[Optional] Output directory filename. Will be created.')
    
    parser.add_argument('-r',  default = False,
                        action='store_true', help='[Optional] Should existing products be reprocessed?')
    
    args = parser.parse_args()
    
    # get in directory
    in_dir = args.dir.expanduser().resolve()
    log.info(f"Starting directory {in_dir}")
    
    # get positions filepath
    pos = args.p
    log.info(f'Position file path: {pos}')

    # get dem and nav filepaths
    dem_fp = args.d
    if dem_fp: dem_fp = dem_fp.expanduser().resolve()
    nav_fp = args.n
    if nav_fp: nav_fp = nav_fp.expanduser().resolve()
    log.info(f'dem filepath: {dem_fp}\nNav filepath: {nav_fp}')

    default_dem_fp = args.default_dem
    if default_dem_fp: default_dem_fp = default_dem_fp.expanduser().resolve()
    log.info(f'Default dem: {default_dem_fp}')

    # get flag for wgs84
    wgs_flag = args.w
    log.info(f'Convert to WGS84: {wgs_flag}')
    
    # get output directory and error log file
    out_dir = args.o
    if out_dir: out_dir = out_dir.expanduser().resolve()
    log.info(f'Output directory: {out_dir}')

    error_log_fp = args.error_log
    log.info(f'Error log file: {error_log_fp}')
    if error_log_fp: error_log_fp = error_log_fp.expanduser().resolve()

    # get reprocessing flag
    reprocess = args.r
    if reprocess == None:
        reprocess = False
    log.info(f'Reprocessing flag: {reprocess}')

    # get position file other wise just use start and end of raw data
    if pos is not None:
        pos = pos.expanduser().resolve()
        log.info(f"Using position file {pos}.")
        start, end = parse_positions(pos)
        log.info(f"Using starting position {start} and end {end}")
        pos = {'start': start, 'end': end}
        
    # this is ugly - find better solution - currently checking in the directory this file is in
    m_path = assert_one_fp(Path(os.path.realpath(__file__)).parent, 'mfunctions')

    # find all raw files inside of in directory
    raws = [GSLRawData(raw_fp, raw_fp.parent.joinpath(f"{'_'.join(raw_fp.stem.split('_')[:-2])}.raw_par")) for raw_fp in in_dir.rglob('*.raw')]
    log.info(f'Found {len(list(raws))} raw files.')
    
    # parents = []
    errors = {}
    successes = []

    pbar = manager.counter(total = len(raws), desc = 'GSL Raw to SLC')

    for raw in raws:

        pbar.update()

        if raw.fp.parent.parent.stem == '20210214':
            continue
    
        # if raw.fp.parent.stem != 'albion':
            # continue
            
        # if raw.fp.parent.stem in parents:
            # continue
        # else:
            # parents.append(raw.fp.parent.stem)

        log.info(f'Starting on raw GSL file: {raw.fp}')
            
        pol_dict = {'R1': 'H', 'R2': 'V'}
        raw.find_polarization(pol_dict)
        log.info(f'Polarization: {raw.pol} using dictionary: {pol_dict}')

        try:
            raw.find_dem(dem_fp, in_dir)
        except ValueError as e:
            raw.find_dem(default_dem_fp, in_dir)

        log.info(f'Using DEM: {raw.dem.fp}')

        if wgs_flag and raw.dem.tif_fp:
            raw.dem.convert_to_wgs()
        
        raw.find_nav(nav_fp, in_dir)
        log.info(f'Using nav file {raw.nav.fp}')

        try:
            raw.get_start_end(raw.nav.fp, pos)
        except KeyError as ke:
            errors[raw.fp] = ke
            continue

        log.info(f'Using start: {raw.start_pos}. End: {raw.end_pos}')
        log.info(f"Starting GPS second {raw.start_gps_sec}. Ending GPS second {raw.end_gps_sec}")
        
        raw.setup_output(out_dir)
        log.info(f'Complex raw directory: {raw.cplx_dir}. Processed dir: {raw.processed_dir}.')

        logfn = raw.processed_dir.joinpath(f'gs_proc_{raw.fp.stem}.log')
        logf = open(logfn, 'w')
        log.info(f'This raw processing logged at {logfn}')
    #         """
    #         raw_gsl_proc(raw = raw.raw_fp, # GSL raw data 
    #                         raw_par = raw.raw_par, # GSL raw data parameter file
    #                         nav_data = nav_fp, # Navigation INS trajectory data file
    #                         dem_par = dem_par_fp, # DEM parameter file
    #                         cplx_dir = cplx_dir,  # complex intermediate directory
    #                         slc_dir = out_dir, # out directory to save slc out directory
    #                         mfunc_dir = m_path, # directory path to 
    #                         start_pos = raw.start_pos, # track start coordinates [lat, lon, altitude] in decimal degrees and meters, ellipsoidal height (m)
    #                         end_pos = raw.end_pos, # track end coordinates [lat, lon, altitude] in decimal degrees and meters, ellipsoidal height (m)
    #                         pol = raw.pol, # polarization of your data
    #                         geoslc = None, # if you want a specific file name for the slc. Otherwise derives from raw data name
    #                         logf = 'test.log', # name of log file
    #                         geoid_hgt = '-', # height of the geoid relative to the WGS84 reference ellipsoid, (if '-' use the egm2008 geoid model)
    #                         hflg = False, # recalculating Hilbert transform of the raw data
    #                         tflg = False, # recalculate the position and velocity at each pulse for the radar phase center
    #                         cflg = False, # recalculating (overwrite) the geo SLC using the TDBP processor
    #                         e_win = 1, # easting multi-look window dimension
    #                         n_win = 1, # northing multi-look window dimension
    #                         scale = 0.6, # power-law display exponent (enter - for default: 0.3)
    #                         exp = 0.3) # power-law display scale factor (enter - for default: 0.6)

    #         """

        try:
            raw_gsl_proc(raw = str(raw.fp), \
                raw_par = str(raw.par_fp), \
                nav_data = str(raw.nav.fp), \
                dem_par = str(raw.dem.par_fp), \
                cplx_dir = str(raw.cplx_dir), \
                slc_dir = str(raw.processed_dir), \
                mfunc_dir = str(m_path), \
                start_pos = raw.start_pos, \
                end_pos = raw.end_pos, \
                pol = raw.pol, \
                geoslc = None, \
                logf = logf, \
                geoid_hgt = 10, \
                hflg = reprocess, \
                tflg = reprocess, \
                cflg = reprocess, \
                e_win = 1, \
                n_win = 1, \
                scale = 0.6, \
                exp = 0.3)
            
            successes.append(raw.fp)
        except Exception as e:
            errors[raw.fp] = e
            log.warn(f'Error in raw file {raw.fp}: {e}')
    
    log.info(f'Executed with {len(successes)} processed and {len(errors)} errors')
    for raw, error in errors.items():
        log.info(f'Filepath: {raw} failed with error: {error}')
    
    
    if error_log_fp:
        with open(error_log_fp, 'w') as f:
            f.writelines(json.dumps({str(key): str(val) for key, val in errors.items()}, indent = 2))
        
        with open(error_log_fp.parent.joinpath('successes.json'), 'w') as f:
            f.writelines(json.dumps([str(s) for s in successes], indent = 2))
