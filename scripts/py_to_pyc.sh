#!/bin/bash

################################################################################
# Compiles all Python source files in a directory (recursively).
#
# Runs quietly if there are no errors. Otherwise, prints errors to console.
################################################################################
COMPILE_LOG='./py_to_pyc.log'
TARGET='.'

if [ "$#" -gt 0 ]; then
	TARGET=$1
fi

echo "Compiling Python source files"
python3 -m compileall -f -q -b ${TARGET}

if [ -s ${COMPILE_LOG} ]; then
	echo "Error(s) occurred while compiling Python scripts!"
	echo "--------------- start: ${COMPILE_LOG} ---------------"
	cat ${COMPILE_LOG}
	echo "---------------- end: ${COMPILE_LOG} ----------------"
	exit 1
fi
