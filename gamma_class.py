import os
from abc import ABC

import numpy as np
import pandas as pd

from gs_proc import parse_param
from gs_funcs import utctoweekseconds, assert_one_fp, execute, check_path_depth

class GammaData(ABC):
    """
    Generic Gamma datafile. Contains file path and par (annotation) file path
    """

    def __init__(self, fp, par_fp):
        self.fp = fp
        self.par_fp = par_fp

    def find_dem(self, dem_fp, in_dir):

        # if we get the dem filepath just add it in
        if dem_fp:
            if dem_fp.suffix == '.tif':
                execute(f'dem_import {dem_fp} {dem_fp.with_suffix('.dem')} {dem_fp.with_suffix('.dem_par')}')
                self.dem = GammaDEM(dem_fp.with_suffix('.dem'), dem_fp.with_suffix('.dem_par'), dem_fp)
            elif dem_fp.suffix == '.dem_par':
                self.dem = GammaDEM(dem_fp.with_suffix('.dem'), dem_fp)
            elif dem_fp.suffix == '.dem':
                self.dem = GammaDEM(dem_fp, dem_fp.with_suffix('.dem_par'))
            else:
                raise ValueError(f"Unable to parse dem filepath of: {dem_fp}. Unknown extension: {dem_fp.suffix}")
            return
        
        loc_name = self.fp.parent.stem
        parent = self.fp.parent
        self.dem = None

        while self.dem is None:
            # check for dems in this directory with our name
            dems = list(parent.rglob(f'{loc_name}.dem_par'))
            # only let us search one levels up
            dems = [d for d in dems if check_path_depth(d, parent) < 2]
            # check for tifs and convert them to the right format
            tifs = list(parent.rglob(f'{loc_name}.tif'))
            # only search one level up
            tifs = [d for d in tifs if check_path_depth(d, parent) < 2]

            # if we find a dem_par file with the right name add it
            if len(dems) == 1:
                if len(tifs) == 1:
                    self.dem = GammaDEM(dems[0].with_suffix('.dem'), dems[0], tif_fp=tifs[0])
                else:
                    self.dem = GammaDEM(dems[0].with_suffix('.dem'), dems[0])
                return

            # if we find a tif only
            if len(tifs) == 1:
                tif_fp = tifs[0]

                # convert to Gamma format
                execute(f'dem_import {tif_fp} {tif_fp.with_suffix('.dem')} {tif_fp.with_suffix('.dem_par')}')

                # add to object
                self.dem = GammaDEM(tif_fp.with_suffix('.dem'), tif_fp.with_suffix('.dem_par'), tif_fp)
                return

            # if we get to our start directory end this loop
            if parent == in_dir:
                raise ValueError(f"No DEM found for {self.fp}. Searching for DEM: {loc_name}.")

            parent = parent.parent

class GSLRawData(GammaData):
    """
    GSL raw data file
    """
    def __init__(self, raw_fp, par_fp):
        super().__init__(raw_fp, par_fp)
        self.add_metadata()
    
    def add_metadata(self):

        params = parse_param(self.par_fp)
        for key in params:
            setattr(self, key, params[key])
        
        self.add_start_time()

        t_naive = pd.to_datetime(self.RAW_start_time.strftime("%Y-%m-%d %H:%M:%S"))
        self.gps_week, self.gps_day, self.start_gps_sec, _ = utctoweekseconds(t_naive, 18)

        self.add_end_time()
    
    def add_start_time(self):
        with open(self.par_fp) as f:
            lns = f.readlines()
        
        setattr(self, lns[0].split(':')[0] ,pd.to_datetime(':'.join(lns[0].replace('\n','').split(':')[1:])))
    
    def find_polarization(self, polarization_dict: dict):
        self.pol = polarization_dict[self.fp.stem.split('_')[-2]]

    def add_end_time(self):
        block_length = int(self.CGEN_num_samp)+1
        sizeof_data = 2
        num_channels = 1
        bytes_per_record = num_channels * sizeof_data * block_length
        filesize = os.stat(self.fp).st_size 
        total_raw_echoes = np.floor(filesize/bytes_per_record)

        # this assumes we don't use alternative TX 
        tcycle = block_length/float(self.ADC_sample_rate)
        prf = abs(1.0/tcycle)
        self.end_gps_sec = self.start_gps_sec + (total_raw_echoes - 1) / (prf)
    
    def get_start_end(self, nav_fp, positions):
        if positions != None:
            self.start_pos = positions['start']
            self.end_pos = positions['end']
        
        else:
            nav_cols = ['time','latitude','longitude','altitude','roll','pitch','heading']
            df = pd.read_csv(nav_fp, comment = '#', names = nav_cols, index_col= 0)
            sub = df.loc[self.start_gps_sec:self.end_gps_sec]
            if len(sub) < 4:
                raise KeyError(f"No data found in nav file: {nav_fp} with start: {self.start_gps_sec} and end {self.end_gps_sec}")
            start, end = sub.iloc[0], sub.iloc[-1]
            self.start_pos = [np.rad2deg(start.latitude), np.rad2deg(start.longitude), start.altitude]
            self.end_pos = [np.rad2deg(end.latitude), np.rad2deg(end.longitude), end.altitude]

    def find_nav(self, nav_fp, in_dir):

        # if we get the dem filepath just add it in
        if nav_fp:
            if '_plus18sec' in nav_fp.suffix:
                self.nav = AtlansNav(nav_fp, True)
            else:
                self.nav = AtlansNav(nav_fp, False)
            return
        
        assert self.RAW_start_time
        parent = self.fp.parent
        self.nav = None
        date = self.RAW_start_time.strftime('%Y%m%d')
        
        yesterday = (self.RAW_start_time - pd.Timedelta('1 day')).strftime('%Y%m%d')

        while self.nav == None:
            navs = list(parent.joinpath('nav').rglob(f'ATLANS*{date}*POSTPROCESSING*.txt_plus18sec'))
            # some are next day due to time zones and running navigation files for a long time
            navs = navs if len(navs) > 0 else list(parent.joinpath('nav').rglob(f'ATLANS*{yesterday}*POSTPROCESSING*.txt_plus18sec'))

            if len(navs) == 1:
                self.nav = AtlansNav(navs[0], True)
                return
            if len(navs) > 1:
                nav_dates = {abs(pd.to_datetime(date) - pd.to_datetime('T'.join(nav.stem.split('_')[0].split('-')[1:]))): nav for nav in navs}
                self.nav = AtlansNav(nav_dates[min(nav_dates)], True)
                return
            
            navs = list(parent.joinpath('nav').rglob(f'ATLANS*{date}*POSTPROCESSING*.txt'))
            # some are next day due to time zones and running navigation files for a long time
            navs = navs if len(navs) > 0 else list(parent.joinpath('nav').rglob(f'ATLANS*{yesterday}*POSTPROCESSING*.txt'))
            if len(navs) == 1:
                self.nav = AtlansNav(navs[0], False)
                return
            if len(navs) > 1:
                nav_dates = {abs(pd.to_datetime(date) - pd.to_datetime('T'.join(nav.stem.split('_')[0].split('-')[1:]))): nav for nav in navs}
                self.nav = AtlansNav(nav_dates[min(nav_dates)], False)
                return
                
            # if we get to our start directory end this loop
            if parent == in_dir:
                raise ValueError(f"No nav found for {self.fp}")
            
            parent = parent.parent
    
    def setup_output(self, out_dir):
        
        # if we get output directory just make it and set up cplx and processed
        if out_dir:
            self.out_dir = out_dir
        # if we have raw then do one out
        elif self.fp.parent.stem == 'raw':
                self.out_dir = self.fp.parent.parent
        # otherwise just put it with the raw files.
        else:
            self.out_dir = self.fp.parent

        self.cplx_dir = self.out_dir.joinpath('raw_cplx')
        self.processed_dir = self.out_dir.joinpath('processed')
        
        for d in [self.out_dir, self.cplx_dir, self.processed_dir]:
            d.mkdir(exist_ok = True)

class GammaDEM(GammaData):
    def __init__(self, dem_fp, par_fp, tif_fp = None):
        super().__init__(dem_fp, par_fp)
        self.tif_fp = tif_fp
    
    def convert_to_wgs(self):
        assert self.tif_fp
        
        self.wgs84_tif = self.tif_fp.with_suffix('.wgs.tif')
        
        if not self.wgs84_tif.exists():
            execute(f"gdalwarp {self.tif_fp} {self.wgs84_tif} -t_srs EPSG:4326 -r bilinear")
        
        if not self.wgs84_tif.with_suffix('.dem').exists():
            execute(f"dem_import {self.wgs84_tif} {self.wgs84_tif.with_suffix('.dem')} {self.wgs84_tif.with_suffix('.dem_par')}")

        self.fp = self.wgs84_tif.with_suffix('.dem')
        self.par_fp = self.wgs84_tif.with_suffix('.dem_par')

class AtlansNav(ABC):
    def __init__(self, fp, utm_2_gps):
        self.fp = fp
        self.gsp_time = utm_2_gps

        if self.gsp_time == False:
            self.add_18_secs()
    
    def add_18_secs(self, utm_2_gps_correction = 18):
            """
            Add 18 seconds to navigation txt file to convert from UTC time to GPS time
            """

            if self.fp.suffix == '.txt_plus18_sec': return

            nav_cols = ['time','latitude','longitude','altitude','roll','pitch','heading']
            df = pd.read_csv(self.fp, comment = '#', names = nav_cols, index_col= 0)

            df.index = df.index + utm_2_gps_correction

            df.to_csv(self.fp.with_suffix('.txt_plus18sec'), header = False)

            self.fp = self.fp.with_suffix('.txt_plus18sec')
            return

class GammaSLC(GSLRawData):
    """
    Gamma Single Look Complex Image.
    """
    pass