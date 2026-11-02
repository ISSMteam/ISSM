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
DATASETS_URL="https://issmteam.github.io/ISSM-Documentation/using-issm/tutorials/datasets"
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
dataset_urls=()
while IFS= read -r url; do
	dataset_urls+=("${url}")
done < <(
	curl -Lks "$DATASETS_URL" |
	perl -0777 -ne '
	if (/<div data-marker="datasets-list-start"><\/div>(.*?)<div data-marker="datasets-list-end"><\/div>/s) {
	print $1;
	}
	' |
	grep -Eo 'href="[^"]*"' |
	sed 's/^href="//; s/"$//'
)

# Skip datasets that have already been downloaded
#
urls_to_download=()
for url in "${dataset_urls[@]}"; do
	[ -z "${url}" ] && continue

	file_name=$(basename "${url%%\?*}")
	[ -z "${file_name}" ] && continue

	if [ -f "${DIRECTORY_PREFIX}/${file_name}" ]; then
		echo "File ${file_name} already downloaded, skipping..."
	else
		urls_to_download+=("${url}")
	fi
done

# Get datasets
#
# NOTE: wget is not available by default on macOS, so fall back to cURL
#
if [ ${#urls_to_download[@]} -eq 0 ]; then
	echo "All examples datasets already downloaded"
else
	echo "Downloading examples datasets..."

	if command -v wget > /dev/null 2>&1; then
		printf '%s\n' "${urls_to_download[@]}" |
		wget --quiet --no-clobber --directory-prefix="${DIRECTORY_PREFIX}" -i -
	else
		for url in "${urls_to_download[@]}"; do
			file_name=$(basename "${url%%\?*}")

			curl -Lk --fail --silent --show-error --output "${DIRECTORY_PREFIX}/${file_name}.part" "${url}" &&
			mv "${DIRECTORY_PREFIX}/${file_name}.part" "${DIRECTORY_PREFIX}/${file_name}" ||
			{
				echo "Error: failed to download ${url}" >&2
				rm -f "${DIRECTORY_PREFIX}/${file_name}.part"
			}
		done
	fi
fi

# Expand zip files
unzip -n -d "${DIRECTORY_PREFIX}" "${DIRECTORY_PREFIX}/*.zip"
