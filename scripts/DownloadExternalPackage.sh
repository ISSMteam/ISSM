#!/bin/bash
#
# DownloadExternalPackage.sh
# Generally, used to download a hosted file from a URL
# URL - Location of file to download
# file - File to write to (including path)
# usage: download_external_package.sh URL file

## Constants
#
MSG_ERR_NO_GET_CMD="No supported file download command was found"
MSG_USAGE="usage: $(basename ${0}) [-h] URL file
  URL  : Location of file to download
  file : File to write to (including path)"

## Variables
#
OUT_FILE=""
URL=""

## Check that number of args is 2 (note that this also handles case where user explicitly requests help)
#
if [ $# != 2 ]
then
	echo "$MSG_USAGE"
	exit 0
fi

## Retrieve args
#
URL=$1
OUT_FILE=$2

## Check if OUT_FILE already exists
#
if [ -f ${OUT_FILE} ]
then
	echo "File ${OUT_FILE} already exists and will not be downloaded..."
	exit 0
fi

## Download file
#
if [ ! -z `which curl` ]; then
	curl --silent $URL -o $OUT_FILE

	# Try wget if curl exists but fails
	if [ $? -ne 0 ]; then
		if [ ! -z `which wget` ]; then
			wget --quiet -O $OUT_FILE $URL

			if [ $? -ne 0 ]; then
				echo "Error: both curl and wget failed. To debug, download package manually using curl without '--silent' option and/or wget without '--quiet' option."
				exit 0
			fi
		else
			echo $MSG_ERR_NO_GET_CMD
			exit 0
		fi
	fi
elif [ ! -z `which wget` ]; then
	wget --quiet -O $OUT_FILE $URL
else
	echo $MSG_ERR_NO_GET_CMD
fi
