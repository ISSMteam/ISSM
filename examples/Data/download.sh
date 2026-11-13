#!/bin/bash

# Use wget if available, fall back to curl otherwise
if command -v wget > /dev/null 2>&1; then
	download() { wget "$1"; }
elif command -v curl > /dev/null 2>&1; then
	download() { echo "Downloading $1..." && curl -O -L "$1"; }
else
	echo "Error: neither wget nor curl is available" >&2
	exit 1
fi

download https://ross.ics.uci.edu/tutorials/SquareShelf.nc
download https://ross.ics.uci.edu/tutorials/Antarctica_5km_withshelves_v0.75.nc
download https://ross.ics.uci.edu/tutorials/Greenland_5km_dev1.2.nc
download https://ross.ics.uci.edu/tutorials/Antarctica_ice_velocity.nc
download https://ross.ics.uci.edu/tutorials/greenland_vel_mosaic250_vx_v1.tif
download https://ross.ics.uci.edu/tutorials/greenland_vel_mosaic250_vy_v1.tif
download https://ross.ics.uci.edu/tutorials/CrossOvers2009.mat
download https://ross.ics.uci.edu/tutorials/Box_Greenland_SMB_monthly_1840-2012_5km_cal_ver20141007.nc
download https://data.cresis.ku.edu/data/grids/old_versions/Jakobshavn_2008_2011_Composite.zip
download https://ross.ics.uci.edu/tutorials/GRACE_and_suporting_datasets.zip
