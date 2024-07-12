#!/bin/bash

# Source config file
if [ $# -ne 1 ]; then
	#no config file specified: exit
	echo "No config file specified. Exiting..." >&2 # Error message to stderr.
	exit 1
fi
if [ ! -f "$1" ]; then
	echo "File $1 not found!" >&2 # Error message to stderr.
	exit 1
fi

# Initialize test suite variables (to avoid "-eq: unary operator expected")
NUMCPUS_INSTALL=1
NUMCPUS_RUN=1

# Source configuration script
source $1;

# Number of packages
NUMPACKAGES=$(($(echo ${EXTERNALPACKAGES} | wc -w ) / 2))
source ${ISSM_DIR}/etc/environment.sh

for ((i=1;i<=$NUMPACKAGES;i++)); do
	NUM1=$((2*$i-1))
	NUM2=$((2*$i))
	PACKAGENAME=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM1-$NUM1)
	PACKAGEINST=$(echo $EXTERNALPACKAGES | cut -d " " -f $NUM2-$NUM2)

	# Install if requested or if previous install has not been successful
	cd ${ISSM_DIR}/externalpackages/$PACKAGENAME

	# Do not install package if it already exists
	PREFIX=$(egrep "PREFIX=" ./$PACKAGEINST | cut -d'"' -f 2 | envsubst)
	echo "======================================================";
	echo "       Installing $PACKAGENAME                        ";
	echo "======================================================";

	./$PACKAGEINST $NUMCPUS_INSTALL &> compil.log
	if [ $? -ne 0 ] && [ "${PACKAGENAME}" != "boost" ]; then
		cat compil.log
		echo "======================================================";
		echo "    ERROR: installation of $PACKAGENAME failed        ";
		echo "======================================================";
		exit 1;
	fi
	source ${ISSM_DIR}/etc/environment.sh
done

# ISSM compilation
cd $ISSM_DIR
echo "======================================================";
echo "                    Reconfiguring                     ";
echo "======================================================";
autoreconf -ivf
if [ $? -ne 0 ]; then
	echo "autoreconf failed!"
	exit 1
fi

eval "./configure ${ISSM_CONFIG}"
if [ $? -ne 0 ]; then
	echo "ISSM configuration failed (see options below)"
	echo $ISSM_CONFIG
	echo "ISSM configuration failed!"
	exit 1
fi

# Compile and install ISSM
echo "======================================================"
echo "                    Compiling ISSM                    "
echo "======================================================"
make -j $NUMCPUS_INSTALL install
if [ $? -ne 0 ]; then
	echo "ISSM_COMPILATION failed!"
	exit 1
fi

echo "======================================================"
echo "             Preparing Test scripts                   "
echo "======================================================"

#Prepare MATLAB tests
cat > ${ISSM_DIR}/matlab_ci.m << EOF
% Go to the test directory
cd $ISSM_DIR/test/NightlyRun/

% Add ISSM tools to path
addpath('${ISSM_DIR}/src/m/dev');
devpath;

% Run tests
runme(${MATLAB_NROPTIONS},'quionerror',1);
disp('MATLABEXITEDCORRECTLY');
quit(0);
EOF
