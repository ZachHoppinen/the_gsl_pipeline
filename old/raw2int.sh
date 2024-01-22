#!/usr/bin/env bash

# set common environment variables 

# start and end positions
start_pos='[65.3319542166, -146.7114936, 706.8]'
end_pos='[65.32216095, -146.6939016, 767.3]'
# height of the geoid relative to the WGS84 elliposoid for your location
#TODO add programatic geoid height search here
geoid_hgt='+11'
# dem par file (use dem_import to create this file from a geotiff)
dem_par='DEM/nc_hd.dem_par'

# channel to polarization relationship
declare -A polarization_channel
polarization_channel=( ["1"]="H" ["2"]="V")

# matlab functions. should either be in the same directory or in the path
mpath='mfunctions_DEMO_GSL_Release-2022-04-25_02-38-55'

# first positional command line argument
DIRECTORY=$1

# check if directory exists and exit if not
if [ ! -d "$DIRECTORY" ]; then
    echo "$DIRECTORY does not exist"
    exit 1
fi

echo "Running in directory: $DIRECTORY"

if [ ! -f $dem_par ]; then
    echo "DEM par does not exist. $dem_par"
fi
echo "Using DEM: $dem_par with geoid_hgt correction: $geoid_hgt"

# activate correct environment
conda activate gamma

# output directory
slc_dir="$DIRECTORY/processed"
echo "Outputting results into $slc_dir"
if [ ! -d $slc_dir ]; then 
    mkdir $slc_dir
fi

# nav file search
# nav_data='20230801/nav/ATLANS-20230308-120754_OutC_POSTPROCESSING-replay.xpf.txt_plus18sec'
nav_data=$(find $DIRECTORY -type f -name "*.txt_plus18sec")
echo "Using nav_data $nav_data"

if [ -z "$nav_data" ]; then
    echo "No 18 second data found. Trying to generate."
    nav_data=$(find $DIRECTORY -type f -name "ATLANS*.txt")
    plus_18_ml=$(find . -type f -name "plus18.m")

    echo $nav_data
    echo $plus_18_ml

    if [ ! -z "$nav_data" ] && [ ! -z "$plus_18_ml" ]; then

        #TODO add check for more than 1 nav data or plus18
        echo $nav_data
        echo $plus_18_ml
    fi
fi

# recursively search for files that end in .raw
RAWS=$(find $DIRECTORY -type f -name "*.raw")
# converts string with spaces to array
RAWS=($RAWS)

# loop through each raw file
for RAW in ${RAWS[@]}; do
    echo "Processing $RAW"
    
    # get the corresponding raw_par file by slicing off the end underscores
    RAW_PAR="${RAW%%_R*}.raw_par"

    if [ ! -f $RAW_PAR ]; then
        echo "Raw par $RAW_PAR does not exist."
        break
    fi

    #TODO check if output file already exists and skip

    echo "Using raw par file $RAW_PAR"

    pol=${RAW##*_R}
    pol=${pol%%_*}
    pol=${polarization_channel[$pol]}
    echo "POL IS: $pol"

    echo "Final command: "
    echo "gs_proc.py $RAW $RAW_PAR $nav_data $dem_par $slc_dir "$start_pos" "$end_pos" $geoid_hgt $pol -m $mpath"
    gs_proc.py $RAW $RAW_PAR $nav_data $dem_par $slc_dir "$start_pos" "$end_pos" $geoid_hgt $pol -m $mpath

done
