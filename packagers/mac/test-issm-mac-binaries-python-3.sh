#!/bin/bash

################################################################################
# This script is intended to test ISSM macOS Python 3 binaries on an end-user 
# machine after successful packaging and signing.
#
# NOTE: Tarball must already exist in INSTALL_DIR
################################################################################

## Constants
#
INSTALL_DIR=${PWD}
PKG="ISSM-macOS-Python-3"
PYTHON_NROPTIONS="--benchmark all --exclude 125 126 129 234 235 418 420 435 444 445 550 701 702 703 1101 1102 1103 1104 1105 1106 1107 1108 1109 1110 1201 1202 1203 1204 1205 1206 1207 1208 1301 1302 1303 1304 1401 1402 1601 1602 2002 2003 2004 2005 2006 2007 2008 2010 2011 2012 2013 2020 2021 2051 2052 2053 2084 2085 2090 2091 2092 2101 2424 2425 3001:3300 3480 3481 4001:4100" # NOTE: Combination of test suites from basic, Dakota, and Solid Earth builds, with tests that require a restart and those that require the JVM excluded

COMPRESSED_PKG="${PKG}.zip"

export ISSM_DIR="${INSTALL_DIR}/${PKG}"
export PATH="${PATH}:${ISSM_DIR}/bin:${ISSM_DIR}/scripts"
export PYTHONPATH="${ISSM_DIR}/scripts"
export PYTHONSTARTUP="${PYTHONPATH}/devpath.py"
export PYTHONUNBUFFERED=1 # We don't want Python to buffer output, otherwise issm.exe output is not captured

cd ${INSTALL_DIR}
rm -rf ${PKG}
ditto -xk ${COMPRESSED_PKG} .
cd ${PKG}/test/NightlyRun

# Run tests, redirecting output to logfile and suppressing output to console
echo "Running tests"
rm python.log 2> /dev/null
./runme.py ${PYTHON_NROPTIONS} &> python.log 2>&1

# Check that Python did not exit in error
pythonExitCode=`echo $?`
pythonExitedInError=`grep -c -E "Error|No such file or directory|Permission denied|Standard exception|Traceback|bad interpreter|syntax error|error:" python.log`

if [[ ${pythonExitCode} -ne 0 || ${pythonExitedInError} -ne 0 ]]; then
	echo "----------Python exited in error!----------"
	cat python.log
	echo "-----------End of python.log-----------"
	exit 1
fi

# Check that all tests passed
numTestsFailed=`grep -c -E "FAILED|ERROR" python.log`

if [[ ${numTestsFailed} -ne 0 ]]; then
	echo "One or more tests FAILED"
	cat python.log
	exit 1
else
	echo "All tests PASSED"
fi
