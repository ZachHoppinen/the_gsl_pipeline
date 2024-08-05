from pathlib import Path

import pandas as pd

import py_gamma as pg

from gs_funcs import get_width

# assign names to DEM, DEM_par, shaded relief intensity image
dem=Path("/data/wcp/dem/wcp.10m.wgs.dem")
dem_par=Path("/data/wcp/dem/wcp.10m.wgs.dem_par")

shade=dem.with_suffix('.dem.shd')
print(shade)
dem_width=get_width(dem_par)
# pg.mapshd(dem, dem_width, 3, 3, '-', '-', shade, 0, 0, 1, 3, 0)

main_dir = Path('/data/wcp/')
SLCs = main_dir.rglob('*.slc')
# print(sorted(list(SLCs)))

from itertools import combinations

int_dir = main_dir.joinpath('ints')
int_dir.mkdir(exist_ok = True)

for slc1, slc2 in combinations(SLCs, 2):

    d1, d2 = [pd.to_datetime(f"{f.stem.split('_')[4]}T{f.stem.split('_')[5]}") for f in [slc1, slc2]]
    if d1 > d2: slc1, slc2 = slc2, slc1

    d1, d2 = [pd.to_datetime(f"{f.stem.split('_')[4]}T{f.stem.split('_')[5]}") for f in [slc1, slc2]]
    mli1, mli2 = slc1.with_suffix('.mli'), slc2.with_suffix('.mli')
    
    d1_d2_str = f"{d1.strftime('%Y-%m-%dT%H-%M')}_{d2.strftime('%Y-%m-%dT%H-%M')}"
    int_fp = int_dir.joinpath(d1_d2_str).with_suffix('.int')
    cc_fp = int_dir.joinpath(d1_d2_str).with_suffix('.cc')
    int_par = int_dir.joinpath(d1_d2_str).with_suffix('.dem_par')

    pg.SLC_intf_geo2(slc1, slc2, dem_par, int_fp, mli1, mli2, cc_fp, int_par, 1, 1, 7, 7, 0)

    print(cc_fp.with_suffix('.cc.bmp'))
    pg.rasdt_pwr(cc_fp, shade, dem_width, 1, 0, 1, 1, .1, 1.0, 0, 'cc.cm' , cc_fp.with_suffix('.cc.bmp'))
    pg.rasmph_pwr(int_fp, shade, dem_width, 1, 0, 1, 1, 'rmg.cm', int_fp.with_suffix('.int.bmp'), 1.5, 1.0, 24)

    pg.mk_kml(int_par, int_fp.with_suffix('.int.bmp'), int_fp.with_suffix('.int.kml'))
    pg.mk_kml(int_par, cc_fp.with_suffix('.cc.bmp'), cc_fp.with_suffix('.cc.kml'))