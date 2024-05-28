#!/bin/bash
set -eu

#Download from ISSM server
$ISSM_DIR/scripts/DownloadExternalPackage.sh 'https://issm.ess.uci.edu/files/externalpackages/kml_shapefile.zip' 'kml_shapefile.zip'
