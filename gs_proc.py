#!/usr/bin/env python3
from cmath import log
import oct2py
import os,sys,getopt
from os.path import join, basename, dirname, abspath, expanduser, exists, splitext
import subprocess
import shlex
import time
import datetime
#import dateutil.parser
import pytz
import ast
import numpy as np
import py_gamma as pg
UTC = pytz.timezone('UTC')  #select time zone (CET, UTC...)

mfunc_dir = None
hflg = False  # default values of the recalculation flags for the hilbert transform, INS trajectory resampling , and TDBP processing 
tflg = False
cflg = False
geoslc = None
scale = 0.6   # raspwr() default scaling factor
exp = 0.3     # raspwr() defaykt exponent
lks = [1,1]   # default number of looks in e_win and n_win

def usage():
  print(f"""
  Calculates a geocoded SLC image using TDBP given a GS-L raw data file and raw_data parameter file, 
  INS data, DEM, and start and end coordinates along the radar track

Usage:
  gs_proc.py <RAW> <RAW_par> <NAV_data> <DEM_par> <SLC_dir> <start_pos> <end_pos> <geoid_hgt> <pol>  Options

  RAW         (input) GS-L raw data
  RAW_par     (input) RAW data parameter file
  NAV_data    (input) INS navigation trajectory data file
              Note: Atlans data requires an +18 second offset
  DEM_par     (input) DEM parameter file
  SLC_dir     directory for the ouput geocoded SLC image
  start_pos   coordinates of the start of the track segment, '[lat, long, elev]', use decimal deg.
  end_pos     coordinates of  the end of the track segment, '[lat, long, elev]', use decimal deg.
  geoid_hgt   height of the geoid relative to the WGS84 reference ellipsoid, (enter -: use the egm2008 geoid model)
  pol         polarization ID included in the output SLC and MLI filenames, e.g., HH, HV, VH, VV

Options:
  -m mfunc_dir       directory containing Matlab/octave functions, the search path can be setup in .octaverc
  -o SLC_file        SLC filename (without path), default: constructed from the RAW file name, and polarization
  -l [e_win, n_win]  easting and northing multi-look factors used in multi_look_geo2(), default: {lks}
  -s scale           power-law display scale factor (enter - for default: 0.6)
  -e exp             power-law display exponent (enter - for default: 0.3)
  -h                 recalculate the Hilbert transform of the raw data
  -t                 recalculate position and velocity and moco files
  -c                 recalculate the TDBP focused radar image
""")
  sys.exit(-1)

def execute(cmd, logf):
  # execute a command in the shell as a subprocess
  # v1.0 29-Jan-2014 clw
    print(cmd)
    logf.write('\n'+cmd+'\n')
    logf.flush()
    args2 = shlex.split(cmd)
    subprocess.call(args2, stdout=logf)
    return 0

def parse_param(fp):
  with open(fp) as f:
    lns = f.readlines()
  d = {}
  ## remove lines with # in front
  for ln in lns:
    if ln[0] != "#":
        vs = ln.split(':')
        if len(vs) == 2:
            d[vs[0].strip()] = vs[1].strip()
  return d

def get_geoid_hgt(clat, clon,):
# 
# get the interpolated value of the the height of the geoid relative to the WGS84 reference ellipsoid from the egm2008-5 geoid model
# clw 24-Oct-2022
  pg.create_dem_par('test.dem_par', '-', '-', '-', '-', 4326, 0, stderr_flag=False, stdout_flag = False)   # use WGS84 ellipsoid
  coord = [clat, clon, 0.000001]
  pg.write_tab(coord, 'coord_in.txt')
  pg.coord_trans_list('coord_in.txt', 'geoid_hgt.txt', 'test.dem_par', 'test.dem_par', pg.which('egm2008-5.dem'), pg.which('egm2008-5.dem_par'), stdout_flag=False)
  ghgt = pg.read_tab('geoid_hgt.txt', dtype=np.float32)
  os.remove('test.dem_par')   # delete temporary files
  os.remove('coord_in.txt')
  return(ghgt[2])   # height of the geoid relative to the ellipsoid (m), e.g. geoidal undulation

def raw_gsl_proc(raw, raw_par, nav_data, dem_par, cplx_dir, slc_dir, geoslc, mfunc_dir, start_pos, end_pos, logf, pol, geoid_hgt, hflg, tflg, cflg, e_win, n_win, scale, exp):
#
# raw           GSL raw data
# raw_par       GSL raw data parameter file
# nav_data      Navigation INS trajectory data file
# dem_par       DEM parameter file
# cplx_dir      directory to contain the Hilbert transforms of the raw data
# slc_dir       directory to contain the processed SLC data
# mfunc_dir     directory containing the Matlab/Octave processing functions
# start_pos     track start coordinates [lat, lon, altitude] in decimal degrees and meters, ellipsoidal height (m)
# end_pos       track end   coordinates [lat, lon, altitude] in decimal degrees and meters, ellipsoidal height (m)
# logf          processing log file object
# pol           polarization ID of the image: HH HV VH VV  (RxTx), the second letter is the polarization of the transmit pulse
# geoid_hgt     height of the geoid relative to the reference ellipsoid
# hflg          recalculate the Hilbert transform of the raw data, if already calculated
# tflg          recalculate the position and velocity at each pulse for the radar phase center
# e_win         easting multi-look window dimension
# n_win         northing multi-look window dimension
# scale         intensity scale factor for raster image generation with raspwr
# exp           intensity exponent for raster image generation with raspwr()
# 
  oc = oct2py.Oct2Py()
  if mfunc_dir != None:    # Octave permits the mfunc_dir to be defined in .octaverc
    oc.addpath(mfunc_dir)
  
  # Constants for the complex Hilbert transform
  # oc.extractGammaSARL2SARData()
  rg_3dB_bw = 40
  az_3dB_bw = 40
  zerosAtNearRg = 300
  
  # Constants for calculation of the position and velocity of the radar phase center for each pulse  
  # oc.prepareGAMMASARHoneyWellHGuideMoco()
  leverArmToAPC = [0,0,0]
  verbose = 0

  # Constants for selecting echo indices
  # oc.selectFirstAndLastEchoBasedOnPositions()
  subaperture_length = 0.0
  minVelThreshold = 4

  # TDBP processing command (requires NVIDIA GPU)
  AZ_PROC_CMD = "az_proc_tdbp_gpu"
  
  # TDBP processing parameters
  elec_delay = 8.375   # range offset including cables used in azimuth focusing
  nav_format = 0  #  0 : two ascii files (1. for pos., 2. for vel.)
  grid_type = 1   #  type of image reconstruction grid: 
                  #    0: SLC geometry (default)
                  #    1: DEM reconstruction grid
                  #    2: custom reconstruction grid (used to process directly
                  #        into the geometry of another SLC-type of reconstruction)
                  
  GEOSLC_format = 0 # SLC_format    SLC output format (enter - for default, from proc_par):
                  #    0: FCOMPLEX, pairs of 4-byte float values
                  #    1: SCOMPLEX, pairs of 2-byte short integers         
  # FFT upsampling ratio
  upsamp_ratio = 4
  # Default settings
  dopcen_azi = "-"
  phase_corr = "-"
  
  # flag indicating whether start-stop-approximation holds
  startstop_flg=0     # 0: start-stop approximation invalid (FMCW case)
                      # 1: start-stop approximation valid (default)                 
  # use defaults
  plist = "-"
  pRC_bandpass1 = "-"
  pRC_bandpass2 = "-"
  # Set fraction of azimuth bandwidth to be processed
  # AZBWFrac = 0.1

  # multi-look parameters, gaussian window, now decimation multi_look_geo2()
  wflg = 2   # circular/ellipsoidal gaussian window
  e_dec = 1  # no decimation in easting and northing
  n_dec = 1  

  # Create sar_par and proc_par, and raw_cplx filenames
  proc_par = join(cplx_dir, basename(raw).replace('.raw','.proc_par'))
  sar_par = join(cplx_dir, basename(raw).replace('.raw','.sar_par'))
  raw_cplx = join(cplx_dir, basename(raw).replace('.raw','.raw_cplx'))   # Hilbert transform
  
  dem = dem_par.replace('.dem_par','.dem')                 # construct DEM filename
  dem_par2 = dem_par.replace('.dem_par','_mli.dem_par')    # DEM par of the MLI image
  nav_dir = dirname(nav_data)    # directory containing the navigation data

  # construct velocity, position, and mocom file names, and calculate those values for every echo in the raw data
  pos_dat = join(nav_dir, basename(raw).replace('.raw', '_pos.dat'))
  vel_dat = join(nav_dir, basename(raw).replace('.raw', '_vel.dat'))
  outputMocoFile = join(cplx_dir, basename(raw).replace('.raw', '.moco'))
  
  # construct idxTxtFile filename that contains the indices of the first and last echo used to calculate the SLC
  idxTxtFile = join(cplx_dir, basename(pos_dat).replace('pos.dat', 'pos.az_ind'))

  # construct SLC file name
  if geoslc is None:      # SLC filename not specified on the command line using the -o flag
    gbn = splitext(basename(raw))[0]   # constract SLC file name from the raw file, remove .raw extension
    
    GEOSLC = join(slc_dir, f'{gbn}_{pol}.slc')
  else:
    gbn = basename(geoslc)  # use the SLC filename specified on the command line using the -o flag
    GEOSLC = join(slc_dir, geoslc)
  
  mlifn = GEOSLC.replace('.slc','.mli')   # construct MLI filename

  print(f'\nMSP SAR_par:     {sar_par}')
  print(f'MSP PROC_par:    {proc_par}')
  print(f'RAW cplx data:   {raw_cplx}')
  print(f'position data:   {pos_dat}')
  print(f'velocity data:   {vel_dat}')
  print(f'index text file: {idxTxtFile}')
  print(f'geocoded SLC:    {GEOSLC}')
  print(f'MLI image:       {mlifn}\n')

  logf.write(f'\nMSP SAR_par:     {sar_par}\n')
  logf.write(f'MSP PROC_par:    {proc_par}\n')
  logf.write(f'RAW cplx data:   {raw_cplx}\n')
  logf.write(f'position data:   {pos_dat}\n')
  logf.write(f'velocity data:   {vel_dat}\n')
  logf.write(f'index text file: {idxTxtFile}\n')
  logf.write(f'geocoded SLC:    {GEOSLC}\n')
  logf.write(f'MLI image:       {mlifn}\n\n')

  print("*** Calculating the Hilbert transform of each radar pulse ***")
# calculate Hilbert transform of the raw data
  if hflg == True:
    for fp in [raw_cplx, proc_par, sar_par]:
      if exists(fp):
        os.remove(fp)

  if exists(raw_cplx) and exists(proc_par) and exists(sar_par):  
    print('Hilbert transform of the raw data hase been preiously calculated, to recalculate use the -h option')
    logf.write("The Hilbert transform of the raw data already exists, to recalculate, use the -h option\n")
  else:
    oc.extractGammaSARL2SARData(raw, raw_par, raw_cplx, proc_par, sar_par, rg_3dB_bw, az_3dB_bw, pol, nout = 0)

# calculate position and velocity of each pulse in the raw data
  print("\n*** Calculating radar position and velocity of each pulse in the raw data ***")
  if tflg == True:
    for fp in [pos_dat, vel_dat, outputMocoFile]: 
      if exists(fp):
        os.remove(fp)
    
  if exists(pos_dat) and exists(vel_dat):    # outputMocoFile may not be calculated by prepareGAMMASARHoneyWellHGuideMoco
    print('The position and velocity of each pulse has been previously calculated, to recalculate use the -t option')
    logf.write('The position and velocity of each pulse has been calculated, to recalculate use the -t option\n') 
  else:
    print(f'oc.prepareGAMMASARHoneyWellHGuideMoco({proc_par}, {nav_data}, {leverArmToAPC}, {outputMocoFile}, {pos_dat}, {vel_dat}, {verbose}, {geoid_hgt}, nout = 0)') 
    oc.prepareGAMMASARHoneyWellHGuideMoco(proc_par, nav_data, leverArmToAPC, outputMocoFile, pos_dat, vel_dat, verbose, geoid_hgt, nout = 0)
    logf.write(f'oc.prepareGAMMASARHoneyWellHGuideMoco({proc_par}, {nav_data}, {leverArmToAPC}, {outputMocoFile}, {pos_dat}, {vel_dat}, {verbose}, {geoid_hgt}, nout = 0)\n') 

  print("\n*** Calculating echo index for the start and end of the selected track ***")
  # calculate the first and last echo index for the synthetic aperture based on the start and end coordinates
  # If this gives warnings of matlab shortcircuit. Then change & to && in all matlab scripts.
  print(f'arrayOfIndices = oc.selectFirstAndLastEchoBasedOnPositions({start_pos}, {end_pos}, {pos_dat}, {vel_dat}, {minVelThreshold}, {subaperture_length}, {idxTxtFile})')
  arrayOfIndices = oc.selectFirstAndLastEchoBasedOnPositions(start_pos, end_pos, pos_dat, vel_dat, minVelThreshold, subaperture_length, idxTxtFile)
  logf.write(f'arrayOfIndices = oc.selectFirstAndLastEchoBasedOnPositions({start_pos}, {end_pos}, {pos_dat}, {vel_dat}, {minVelThreshold}, {subaperture_length}, {idxTxtFile})\n')
  oc.exit()
  
  first_az_ind = int(arrayOfIndices[0][0])
  end_az_ind = int(arrayOfIndices[0][1])
  print(f'indices of the first and last echo in the synth. aperture: {arrayOfIndices}')
  assert first_az_ind < end_az_ind, f"First index {first_az_ind} is larger than end index {end_az_ind}"

  # TDBP focusing to calcuate the SAR image in map coordinates 
  print("\n*** TDBP focusing to calculate the SAR image in map coordinates ***")
  logf.write(f'indices of the first and last echo in the synth. aperture: {arrayOfIndices}\n')
  logf.write("*** TDBP focusing to calculate the SAR image in map coordinates ***\n")
  
  if cflg == True:
    for fp in [GEOSLC,]: 
      if exists(fp):
        os.remove(fp)

  TDBP_cmd = f'{AZ_PROC_CMD} {sar_par} {proc_par} {raw_cplx} {nav_format} {pos_dat} {vel_dat}\
  {grid_type} {dem} {dem_par} {GEOSLC} {GEOSLC_format} {upsamp_ratio} {dopcen_azi} {phase_corr} {startstop_flg}\
  {first_az_ind} {end_az_ind} {plist} {pRC_bandpass1} {pRC_bandpass2} {elec_delay}'

  if exists(GEOSLC):    # check if the GEOSLC already exists
    print('The TDBP focused SLC already exists, set the -c option to recalculate')
    logf.write('NOTE: TDBP focused SLC already exists, to recalculate use the -c option\n')
  else:
    start_time = time.perf_counter()
    execute(TDBP_cmd, logf)
    stop_time = time.perf_counter()
    print(f"{AZ_PROC_CMD} processing time: {stop_time - start_time} s.***")
    logf.write(f"{AZ_PROC_CMD} processing time: {stop_time - start_time} s.***\n")

  dp = parse_param(dem_par)    # get the width of the focused SLC
  dem_width = dp['width']
  
  print("\n*** Calculating geocoded MLI image ***")
  print(f'MLI averaging window e_win: {e_win}  n_win: {n_win}')
  logf.write(f'MLI circular gaussian averaging window e_win: {e_win} n_win: {n_win}\n')
  execute(f'multi_look_geo2 {GEOSLC} {dem_par} {mlifn} {dem_par2} {e_dec} {n_dec} {e_win} {n_win} {wflg} 1 ',logf)
  execute(f'raspwr {mlifn} {dem_width} - - 1 1 {scale} {exp} grey.cm {mlifn}.bmp', logf)   # raster image is in BMP format
  if dp["DEM_projection"] == "EQA":
    execute(f'mk_kml {dem_par} {mlifn}.bmp {mlifn}.kml', logf)   # create a KML file in addition to the geotiff for display with googleearth
  execute(f'data2geotiff {dem_par} {mlifn}.bmp 0 {mlifn}_geo.tif', logf)
  return 0

if __name__ == '__main__': 
  print("*** GS-L processing of a single scene v1.5 25-Oct-2022 cw/zk ***")
  cmd_str = " ".join(sys.argv)
  #print('sys.argv:',sys.argv, len(sys.argv))
  #print(sys.argv[10:])

  opts, args = getopt.gnu_getopt(sys.argv[10:], "htco:m:s:e:l:")

  if len(sys.argv) < 10:
    usage()
  # except getopt.GetoptError as err:
  #   print(str(err))     #prints something like "option -a not recognized")
  #   usage()
  gs_start_time = datetime.datetime.now(UTC)
  print('*** gs_proc.py processing started %s *** '%(gs_start_time))

  raw = sys.argv[1]
  raw_par = sys.argv[2]
  nav_data = sys.argv[3]
  dem_par = sys.argv[4]
  slc_dir = sys.argv[5]
  start_pos = str(sys.argv[6])
  end_pos =str(sys.argv[7])
  geoid_hgt = str(sys.argv[8])
  pol = str(sys.argv[9])
  rbn = splitext(basename(raw))[0]

  for o, a in opts:
    if o == '-m':
      mfunc_dir =  abspath(expanduser(str(a)))
      print(f'Matlab/Octave TDBP processing functions directory: {mfunc_dir}')
    elif o == '-h':
      hflg = True
      print('recalculating Hilbert transform of the raw data')
    elif o == '-t':
      tflg = True
      print('recalculating position, velocity, and moco files from the INS trajectory data')
    elif o == '-c':
      cflg = True
      print('recalculating the geo SLC using the TDBP processor')
    elif o == '-o':
      geoslc = str(a)
    elif o == '-s':
      scale = float(a)
    elif o == '-e':
      exp = float(exp)
    elif o == '-l':
      try:
        lks = ast.literal_eval(a)
      except:
        print("\nERROR: invalid format for specifing MLI easting and northing looks:",a)
        sys.exit(-1)
    else:
      print("ERROR: invalid option on the command line: ",o)
      usage()
  
  e_win = lks[0]
  n_win = lks[1]

  if tflg is True:   # if the track is recalculated, then recalculate the SLC image
    print('\n NOTE: recalculating the SLC after the track is recalculated!')
    cflg = True

  start_pos = ast.literal_eval(start_pos)
  end_pos = ast.literal_eval(end_pos)

  #start_pos = start_pos.replace('[','').replace(']','').split(',')
  #end_pos = end_pos.replace('[','').replace(']','').split(',')

  start_pos = [float(c) for c in start_pos]
  end_pos = [float(c) for c in end_pos]
  
 # check if files required for processing exist
  for fp in [raw, raw_par, nav_data, dem_par]:
    assert exists(fp), f'File or directory does not exist {fp}'

  base_dir = dirname(dirname(raw))   # base directory
  # create directory for Hilbert transform of the raw data, same level as the raw directory
  
  cplx_dir = join(base_dir, 'raw_cplx')
  print(f'directory for Hilbert transforms of the raw data: {cplx_dir}')
  if not os.path.isdir(cplx_dir):
    try:
      os.makedirs(cplx_dir, exist_ok= True)
    except:
      print('ERROR: cannot create directory for the raw_cplx data: %s'%cplx_dir)
      sys.exit(-1)  
  
  print('output directory for geocoded SLC images : %s'%(slc_dir,))  # SLC directory
  
  if not os.path.isdir(slc_dir):
    try:
      os.makedirs(slc_dir, exist_ok=True)
    except:
      print('ERROR: cannot create output directory for the SLC images: %s'%slc_dir)
      sys.exit(-1)
    print('created output SLC directory: %s'%slc_dir)
  
  logfn = join(slc_dir, f'gs_proc_{rbn}.log')   # use raw data basename to create the log file name 

  if geoid_hgt != '-':
    try:
      geoid_hgt = float(geoid_hgt)
    except:
      print("ERROR: cannot interpret geoidal height entered on the command line as a floating point number: %s",geoid_hgt)
      sys.exit(-1)
    egm_model=False
  else:
    rp = pg.ParFile(raw_par)
    gc = rp.get_value('geographic_coordinates')
    clat = float(gc[0])
    clon = float(gc[1])
    geoid_hgt = get_geoid_hgt(clat, clon)  # no screen output
    egm_model = True

  print(f'{cmd_str}\n')
  print('*** gs_proc.py processing started %s *** '%(gs_start_time))
  print(f'RAW DATA:             {raw}')
  print(f'RAW_cplx_dir:         {cplx_dir}')
  print(f'DEM_par:              {dem_par}')
  print(f'NAV INS data:         {nav_data}')
  print(f'gs_proc.py log file:  {logfn}')
  print(f'track start position: {start_pos}')
  print(f'track end position:   {end_pos}')
  if egm_model is False:
    print(f'geoidal height (m):   {geoid_hgt:.3f}')
  else:
    print(f'geoidal height (m):   {geoid_hgt:.3f} (EGM2008)')
  print(f'gs_proc.py log file:  {logfn}')

  logf = open(logfn, 'w') 
  logf.write(f'{cmd_str}\n\n')
  logf.write(f'RAW DATA:             {raw}\n')
  logf.write(f'RAW_cplx_dir:         {cplx_dir}\n')
  logf.write(f'DEM_par:              {dem_par}\n')
  logf.write(f'NAV INS data:         {nav_data}\n')
  logf.write(f'gs_proc.py log file:  {logfn}')
  logf.write(f'track start position: {start_pos}\n')
  logf.write(f'track end position:   {end_pos}\n')
  if egm_model is False:
    logf.write(f'geoidal height (m):   {geoid_hgt:.3f}\n')
  else:
    logf.write(f'geoidal height (m):   {geoid_hgt:.3f} (EGM2008)\n')
  logf.write(f'gs_proc.py log file:  {logfn}\n')
  logf.flush()

  raw_gsl_proc(raw=raw, raw_par=raw_par, nav_data=nav_data, dem_par=dem_par, cplx_dir=cplx_dir, slc_dir=slc_dir, geoslc=geoslc,\
  mfunc_dir=mfunc_dir, start_pos=start_pos, end_pos=end_pos, logf=logf, \
  pol=pol, geoid_hgt=geoid_hgt, hflg=hflg, tflg=tflg, cflg=cflg, e_win=e_win, n_win=n_win, scale=scale, exp=exp)
  
  gs_end_time=datetime.datetime.now(UTC)

  logf.write('\n*** gs_proc.py processing completed %s ***\n'%(gs_end_time,))
  logf.write(f"elapsed time: {gs_end_time - gs_start_time}\n")
  print('\n*** gs_proc.py processing completed %s *** '%(gs_end_time,))
  print(f"*** elapsed time: {gs_end_time - gs_start_time} ***")
