#!/bin/bash
set -eu

#Some cleanup
rm -rf install

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/PolygonClipper.zip' 'PolygonClipper.zip'

#install
mkdir install
cd install
mv ../PolygonClipper.zip .

#uncompress
unzip PolygonClipper.zip

#Make
mex  -compatibleArrayDims gpc.c gpc_mexfile.c -O -output PolygonClip
