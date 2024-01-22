#!/usr/bin/env bash

## Example of using raw2slc.py
# Prior commands
# create DEM
# dem_import nomecreekdem.tif nomecreekdem.dem nomecreekdem.dem_par

# add 18 seconds

# set common environment variables 
start_pos='[65.3319542166, -146.7114936, 706.8]'
end_pos='[65.32216095, -146.6939016, 767.3]'
mpath='mfunctions_DEMO_GSL_Release-2022-04-25_02-38-55'
geoid_hgt='+11'   # height of the geoid relative to the WGS84 elliposoid for your location
dem_par='DEM/nc.dem_par'
slc_dir='20230802/processed'
nav_data='20230802/nav/ATLANS-20230802-161326_OutC_POSTPROCESSING-replay.xpf.txt_plus18sec'

# The H-pol receive antenna was connected to R1
# The V-pol receive antenna was connected to R2

# -h, -t, and -c options force recalculation of the track, hilbert transform, and focusing of the data

raw='20230802/raw/PF_100MHz_H_6_20230802_205148_H_R1_1.raw'
raw_par='20230802/raw/PF_100MHz_H_6_20230802_205148_H.raw_par'
pol='HH'
gs_proc.py $raw $raw_par $nav_data $dem_par $slc_dir "$start_pos" "$end_pos" $geoid_hgt $pol -m $mpath