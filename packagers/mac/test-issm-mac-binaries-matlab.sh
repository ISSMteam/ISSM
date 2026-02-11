#!/bin/bash

################################################################################
# This script is intended to test ISSM macOS MATLAB binaries on an end-user 
# machine after successful packaging and signing.
#
# NOTE: Tarball must already exist in INSTALL_DIR
################################################################################

## Constants
#
INSTALL_DIR=${PWD}
MATLAB_NROPTIONS="'benchmark','all','exclude',[125,126,129,234,235,418,420,435,444,445,550,701,702,703,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1201,1202,1203,1204,1205,1206,1207,1208,1301,1302,1303,1304,1401,1402,1601,1602,2002,2003,2004,2006,2007,2008,2010,2011,2012,2013,2020,2021,2051,2052,2053,2084,2085,2090,2091,2092,2101,2424,2425,3001:3300,3480,3481,4001:4100]" # NOTE: Combination of test suites from basic, Dakota, and Solid Earth builds, with tests that require a restart and those that require the JVM excluded
MATLAB_PATH="/Applications/MATLAB_R2022b.app"
PKG="ISSM-macOS-MATLAB"

COMPRESSED_PKG="${PKG}.zip"

## Environment
#
export ISSM_DIR="${INSTALL_DIR}/${PKG}"
export PATH="${PATH}:${ISSM_DIR}/bin:${ISSM_DIR}/scripts"

cd ${INSTALL_DIR}
rm -rf ${PKG}
ditto -xk ${COMPRESSED_PKG} .
cd ${PKG}/test/NightlyRun

# Run tests, redirecting output to logfile and suppressing output to console
echo "Running tests"
rm matlab.log 2> /dev/null
${MATLAB_PATH}/bin/matlab -nosplash -nodesktop -nojvm -r "try, addpath ../../bin; addpath ../../lib; runme(${MATLAB_NROPTIONS}); exit; catch me,fprintf('%s',getReport(me)); exit; end" -logfile matlab.log &> /dev/null

# Check that MATLAB did not exit in error
matlabExitCode=`echo $?`
matlabExitedInError=`grep -c -E "Activation cannot proceed|Error in|Illegal|Invalid MEX-file|Warning: Name is nonexistent or not a directory" matlab.log`

if [[ ${matlabExitCode} -ne 0 || ${matlabExitedInError} -ne 0 ]]; then
	echo "----------MATLAB exited in error!----------"
	cat matlab.log
	echo
	echo "-----------End of matlab.log-----------"
	exit 1
fi

# Check that all tests passed
numTestsFailed=`grep -c -E "FAILED|ERROR" matlab.log`

if [[ ${numTestsFailed} -ne 0 ]]; then
	echo "One or more tests FAILED"
	cat matlab.log
	exit 1
else
	echo "All tests PASSED"
fi
