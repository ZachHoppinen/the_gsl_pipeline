#!/usr/bin/env bash

# assign names to DEM, DEM_par, shaded relief intensity image
DEM="DEM/nc_hd.dem"
DEM_par="DEM/nc_hd.dem_par"
shd="DEM/nc_hd.dem.shd"
demWidth=`get_value ${DEM_par} width`
mapshd $DEM $demWidth 3 3 - - $shd 0 0 1 3 1

# calculate VV interferogram, 7x7 rectangular averaging window, no decimation
slc1_VV='20230802/processed/PF_100MHz_H_6_20230802_220138_H_R1_1_HH.slc'
slc2_VV='20230802/processed/PF_100MHz_H_6_20230802_205833_H_R1_1_HH.slc'
intf_VV='20230802/spatialbaseline_ints/HH_20230802_205833_220138.int'
mli1_VV='20230802/spatialbaseline_ints/HH_20230802_205833_220138.mli'
mli2_VV='20230802/spatialbaseline_ints/HH_20230802_205833_220138.mli'
cc_VV='20230802/spatialbaseline_ints/HH_20230802_205833_220138.cc'
DEM_par2='20230802/spatialbaseline_ints/HH_20230802_205833_220138.dem_par'

SLC_intf_geo2 $slc1_VV $slc2_VV $DEM_par $intf_VV $mli1_VV $mli2_VV $cc_VV $DEM_par2 1 1 7 7 0
rasdt_pwr $cc_VV $shd $demWidth 1 0 1 1 .1 1.0 0 cc.cm -
rasmph_pwr $intf_VV $shd $demWidth 1 0 1 1 rmg.cm - 1.5 1. 24
  
mk_kml $DEM_par2 $intf_VV".bmp" $intf_VV".kml"
mk_kml $DEM_par2 $cc_VV".bmp" $cc_VV".kml"

