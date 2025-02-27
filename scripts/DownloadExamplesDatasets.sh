#!/bin/bash
################################################################################
# This script downloads all datasets needed for running the ISSM tutorials.
#
# The default behavior is to download datasets to the examples/Data directory 
# relative to this script. An alternate output directory can be designated by 
# supplying a command line argument.
#
# NOTE:
# - This script does not clobber existing files, intentionally. To download and 
#	unzip new copies, first remove existing files manually.
################################################################################

## Constants
#
DATASETS_URL="https://issm.jpl.nasa.gov/documentation/tutorials/datasets"
DIRECTORY_PREFIX=$(cd $(dirname "$0"); pwd)"/../examples/Data" # Default behavior is to download datasets to examples/Data directory relative to this script

if [ $# -gt 0 ]; then
	DIRECTORY_PREFIX=$1

	if [ ! -d "${DIRECTORY_PREFIX}" ]; then
		mkdir -p "${DIRECTORY_PREFIX}"
	fi
fi

# Get content of page that hosts datasets, reduce to just datasets list, then
# parse out dataset links
#
# NOTE: Clear DYLD_LIBRARY_PATH in case we have installed our own copy of cURL
#		and $ISSM_DIR/etc/environment.sh has been sourced as there may be a
#		conflict between versions of cURL executable and libcurl
#
dataset_urls=$(\
	DYLD_LIBRARY_PATH=""; \
	/usr/bin/curl -Lks ${DATASETS_URL} |\
	sed '/<!--DATASETS LIST START-->/,/<!--DATASETS LIST END-->/ !d' |\
	sed -n 's/.*<li><a href="\([^"]*\)">.*/\1/p'
)

# Get datasets
#
echo "Downloading examples datasets..."
wget --quiet --no-clobber --directory-prefix="${DIRECTORY_PREFIX}" -i - <<< "${dataset_urls}"

# Expand zip files
unzip -n -d "${DIRECTORY_PREFIX}" "${DIRECTORY_PREFIX}/*.zip"
