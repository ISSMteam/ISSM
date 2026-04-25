#!/bin/bash
################################################################################
# This script determines coverage of Python API compared to MATLAB API, 
# assuming that number of MATLAB scripts will always be >= number of Python
# scripts.
#
# TODO:
# - Add coverage reports for JavaScript, Julia.
# - Generalize by allowing two or more file extension arguments to count and 
#	compare?
################################################################################

## Functions
#
function recursive_process() {
	local dir="${1}"
	for file in "${dir}"/*; do
		if [[ -f "${file}" && "${file}" == */test*.m ]]; then
			((++num_m))
			basename=$(basename "${file}" ".m")
			if [[ -f "${dir}/${basename}.py" ]]; then
				((++num_py))
			else
				missing+=("${file}")
			fi
		elif [ -d "${file}" ]; then
			recursive_process "${file}"
		fi
	done
}

## Variables
#
missing=()
num_m=0
num_py=0

# Process arguments
if [ "$#" -gt 0 ]; then
	TARGET=$1
else
	echo "Error: please supply the path to the directory whose contents you want to analyze."
	exit 1
fi

## Main
#
recursive_process "${TARGET}"

# Print report
coverage=$(printf '%.2f' $(echo "scale=4; $num_py / $num_m * 100" | bc))
echo "In directory ${TARGET} there are ${num_m} MATLAB scripts and ${num_py} Python scripts for a Python API coverage of ${coverage}%"
#echo "The following MATLAB scripts are missing a Python translation..."
#printf "%s\n" "${missing[@]}"
