#!/bin/bash
################################################################################
# This script tests if MATLAB license will expire soon on Linux and macOS.
#
# NOTE:
# - Assumes constant MATLAB_PATH is defined and points to MATLAB executable or 
#	application.
################################################################################

# Check if MATLAB exists
if ! [ -d ${MATLAB_PATH} ]; then
	echo "${MATLAB_PATH} does not point to a MATLAB installation! Please modify MATLAB_PATH variable in Jenkins job and try again."
	exit 1
fi

# Start MATLAB
${MATLAB_PATH}/bin/matlab -nojvm -nosplash -nojvm -r "exit;" &> matlab.log

# Check log for license expiration message
matlabLicenseExpiring=`grep -c -E "Your license will expire in" matlab.log`

if [ ${matlabLicenseExpiring} -ne 0 ]; then
	cat matlab.log
	echo "Login to build node to update MATLAB license."
	exit 1;
fi
