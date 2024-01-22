#!/home/rdcrlzh1s/miniforge3/envs/gamma/bin/python

import os
import argparse
from pathlib import Path
import subprocess
import shlex
import ast

import numpy as np
import pandas as pd

from typing import Union
import logging

logging.basicConfig()
log = logging.getLogger(__name__)
log.setLevel(logging.INFO)

def check_path_depth(fp, parent):
    depth = 0
    while fp.parent != parent:
        depth += 1
        fp = fp.parent
        if fp == fp.parent:
            raise ValueError(f"Parent: {parent} not found in file path of {fp}")
    return depth

def execute(cmd):
  # execute a command in the shell as a subprocess
  # v1.0 29-Jan-2014 clw
    log.info('\n'+cmd+'\n')
    args2 = shlex.split(cmd)
    subprocess.call(args2)
    return 0

def add_18_secs(fp):
    """
    Add 18 seconds to navigation txt file to convert from UTC time to GPS time
    """

    if not fp.with_suffix('.txt_plus18sec').exists():

        nav_cols = ['time','latitude','longitude','altitude','roll','pitch','heading']
        df = pd.read_csv(fp, comment = '#', names = nav_cols, index_col= 0)

        df.index = df.index + 18

        df.to_csv(fp.with_suffix('.txt_plus18sec'), header = False)

    return fp.with_suffix('.txt_plus18sec')

def assert_one_fp(search_dir :Path, search: str, exclude: Union[str, None] = None):
    """
    Takes a directory and searches that path with a search string optionally excluding a different string.
    Returns path if there is only filepath otherwise raises ValueError.
    """

    if not exclude:
        fp_list = [fp for fp in search_dir.glob('*') if search in fp.name.lower()]
    else:
        fp_list = [fp for fp in search_dir.glob('*') if search in fp.name.lower() and exclude not in fp.name.lower()]

    if len(fp_list) == 1:
        return fp_list[0]

    elif (len(fp_list)) == 0:
        raise ValueError("No filepaths found...")

    else:
        log.debug(f"Found {len(fp_list)} number of filepaths. Listing:")
        for fp in fp_list:
            log.debug(fp)
        raise ValueError(f"{len(fp_list)} filepath found. Please specify correct filepath.")

def get_date_fps(date_dir: Path):
    """
    Processes a date directory to return Raw, Nav datafiles.
    """
    nav_dir = assert_one_fp(date_dir, 'nav')
    log.info(f"Nav directory: {nav_dir}")
    raw_dir = assert_one_fp(date_dir, 'raw', 'cplx')
    log.info(f"Raw directory: {raw_dir}")
    cplx_dir = date_dir.joinpath('raw_cplx')
    log.info(f"Cplx directory: {cplx_dir}")
    processed_dir = date_dir.joinpath('processed')
    log.info(f"Output directory: {processed_dir}")

    raw_fps = raw_dir.glob('*.raw')
    raws = [RawData(raw_fp, raw_dir.joinpath(f"{'_'.join(raw_fp.stem.split('_')[:-2])}.raw_par")) for raw_fp in raw_fps]

    nav_fp = assert_one_fp(nav_dir, '.xpf.txt', '_plus18sec')

    return raws, nav_fp, cplx_dir, processed_dir

def get_dem_fp(dem_directory):
    """
    Function to get dem fp and if only tifs then to convert and return dem and dem par
    """
    dem_pars = list(dem_dir.glob('*.dem_par'))
    dem_pars = [d for d in dem_pars if 'mli' not in d.name]

    if len(dem_pars) == 0:
        log.info(f"No .dem_par file found in {dem_dir}. Trying to find tifs to convert.")
        tif_fp = assert_one_fp(dem_dir, '.tif')
        execute(f'dem_import {tif_fp} {tif_fp.with_suffix('.dem')} {tif_fp.with_suffix('.dem_par')}')
        return tif_fp.with_suffix('.dem_par')
    
    if len(dem_pars) == 1:
        return dem_pars[0]
    
    else:
        log.debug(f"Found {len(dem_pars)} number of filepaths. Listing:")
        for fp in dem_pars:
            log.debug(fp)
        raise ValueError(f"{len(dem_pars)} filepath found. Please specify correct filepath.")

def utctoweekseconds(utc,leapseconds):
    """ 
    Returns the GPS week, the GPS day, and the seconds 
    and microseconds since the beginning of the GPS week 
    https://stackoverflow.com/questions/45422739/gps-time-in-weeks-since-epoch-in-python    
    """
    import datetime, calendar
    datetimeformat = "%Y-%m-%d %H:%M:%S"
    epoch = datetime.datetime.strptime("1980-01-06 00:00:00",datetimeformat)
    tdiff = utc - epoch  + datetime.timedelta(seconds=leapseconds)
    gpsweek = tdiff.days // 7 
    gpsdays = tdiff.days - 7*gpsweek         
    gpsseconds = tdiff.seconds + 86400* (tdiff.days -7*gpsweek) 
    return gpsweek, gpsdays, gpsseconds, tdiff.microseconds

def parse_positions(fp):
    """
    parse position file to start and end list
    """
    with open(fp, 'r') as f:
            lns = f.readlines()

    assert len(lns) == 2, f"Did not find two coordinates in position file {fp}"

    start, end = lns
    start, end = ast.literal_eval(start.split('=')[1]), ast.literal_eval(end.split('=')[1])

    return start, end