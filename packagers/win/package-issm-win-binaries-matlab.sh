#!/bin/bash

################################################################################
# Packages and tests ISSM distributable package for Windows with MATLAB API.
#
# Options:
# -s/--skiptests		Skip testing during packaging Use if packaging fails 
#						for some reason but build is valid.
#
# NOTE:
# - Assumes that the following constants are defined,
#
#		COMPRESSED_PKG
#		ISSM_DIR
#		PKG
#
# See also:
# - packagers/win/complete-issm-win-binaries-matlab.sh
################################################################################

# Expand aliases within the context of this script
shopt -s expand_aliases

## Override certain aliases
#
alias grep=$(which grep)

## Constants
#
#LIBGFORTRAN="/usr/lib/x86_64-win-gnu/libgfortran.so.5.0.0" # Important that this is the library itself
#LIBGFORTRAN_DIST="${ISSM_DIR}/lib/libgfortran.so.5" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
MATLAB_NROPTIONS="'exclude',[IdFromString('Dakota'),125,126,129,435,550,701,702,703]"
#MATLAB_NROPTIONS="'benchmark','all','exclude',[125,126,129,234,235,418,420,435,444,445,701,702,703,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1201,1202,1203,1204,1205,1206,1207,1208,1301,1302,1303,1304,1401,1402,1601,1602,2002,2003,2004,2006,2007,2008,2010,2011,2012,2013,2020,2021,2051,2052,2053,2084,2085,2090,2091,2092,2101,2424,2425,3001:3300,3480,3481,4001:4100]" # NOTE: Combination of test suites from basic, Dakota, and Solid Earth builds, with tests that require a restart and those that require the JVM excluded
MATLAB_PATH=$(cygpath -u $(cygpath -ms "/c/Program Files/MATLAB/R2023b"))
MSMPI_ROOT="${ISSM_DIR}/externalpackages/msmpi/install"

## Environment
#
export PATH="${ISSM_DIR}/bin:$(getconf PATH)" # Ensure that we pick up binaries from 'bin' directory rather than 'externalpackages'

## Parse options
#
if [ $# -gt 1 ]; then
	echo "Can use only one option at a time"
	exit 1
fi

skip_tests=0

if [ $# -eq 1 ]; then
	case $1 in
		-s|--skiptests) skip_tests=1;					;;
		*) echo "Unknown parameter passed: $1"; exit 1	;;
	esac
fi

# Check if MATLAB exists
if ! [ -d ${MATLAB_PATH} ]; then
	echo "${MATLAB_PATH} does not point to a MATLAB installation! Please modify MATLAB_PATH variable in $(basename $0) and try again."
	exit 1
fi

# Clean up from previous packaging
echo "Cleaning up existing assets"
cd ${ISSM_DIR}
rm -rf ${PKG} ${COMPRESSED_PKG}
mkdir ${PKG}

# Add required binaries and libraries to package and modify them where needed
cd ${ISSM_DIR}/bin

echo "Modify generic"
cat generic_static.m | sed -e "s/generic_static/generic/g" > generic.m

echo "Moving system and MinGW libraries to bin/"
cp /c/msys64/mingw64/bin/libgcc_s_seh-1.dll .
cp /c/msys64/mingw64/bin/libstdc++-6.dll .
cp /c/msys64/mingw64/bin/libwinpthread-1.dll .
cp /c/msys64/mingw64/bin/libgfortran-5.dll .
cp /c/msys64/mingw64/bin/libquadmath-0.dll .

echo "Moving MSMPI binaries and libraries to bin/"
if [ -f ${MSMPI_ROOT}/bin/mpiexec.exe ]; then
	cp ${MSMPI_ROOT}/bin/* .
	cp ${MSMPI_ROOT}/lib/libmsmpi.dll .
	cp ${MSMPI_ROOT}/lib/msmpi.dll .
else
	echo "MSMPI not found"
	exit 1
fi

echo "Moving MEX-files to bin/"
mv ${ISSM_DIR}/lib/*.mexw64 .

# echo "Moving MPICH binaries to bin/"
# if [ -f ${ISSM_DIR}/externalpackages/petsc/install/bin/mpiexec ]; then
# 	cp ${ISSM_DIR}/externalpackages/petsc/install/bin/mpiexec .
# 	cp ${ISSM_DIR}/externalpackages/petsc/install/bin/hydra_pmi_proxy .
# elif [ -f ${ISSM_DIR}/externalpackages/mpich/install/bin/mpiexec ]; then
# 	cp ${ISSM_DIR}/externalpackages/mpich/install/bin/mpiexec .
# 	cp ${ISSM_DIR}/externalpackages/mpich/install/bin/hydra_pmi_proxy .
# else
# 	echo "MPICH not found"
# 	exit 1
# fi

# echo "Moving GDAL binaries to bin/"
# if [ -f ${ISSM_DIR}/externalpackages/gdal/install/bin/gdal-config ]; then
# 	cp ${ISSM_DIR}/externalpackages/gdal/install/bin/gdalsrsinfo .
# 	cp ${ISSM_DIR}/externalpackages/gdal/install/bin/gdaltransform .
# else
# 	echo "GDAL not found"
# 	exit 1
# fi

# echo "Moving GMT binaries to bin/"
# if [ -f ${ISSM_DIR}/externalpackages/gmt/install/bin/gmt-config ]; then
# 	cp ${ISSM_DIR}/externalpackages/gmt/install/bin/gmt .
# 	cp ${ISSM_DIR}/externalpackages/gmt/install/bin/gmtselect .
# else
# 	echo "GMT not found"
# 	exit 1
# fi

# echo "Moving Gmsh binaries to bin/"
# if [ -f ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh ]; then
# 	cp ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh .
# else
# 	echo "Gmsh not found"
# 	exit 1
# fi

cd ${ISSM_DIR}/lib

# echo "Moving libgfortran to lib/"
# cp ${LIBGFORTRAN} ${LIBGFORTRAN_DIST} 2> /dev/null

# echo "Moving GSHHG assets to share/"
# if [ -d ${ISSM_DIR}/externalpackages/gmt/install/share/coast ]; then
# 	mkdir ${ISSM_DIR}/share 2> /dev/null
# 	cp -R ${ISSM_DIR}/externalpackages/gmt/install/share/coast ${ISSM_DIR}/share
# else
# 	echo "GSHHG not found"
# 	exit 1
# fi

# echo "Moving PROJ assets to share/"
# if [ -d ${ISSM_DIR}/externalpackages/proj/install/share/proj ]; then
# 	mkdir ${ISSM_DIR}/share 2> /dev/null
# 	cp -R ${ISSM_DIR}/externalpackages/proj/install/share/proj ${ISSM_DIR}/share
# else
# 	echo "PROJ not found"
# 	exit 1
# fi

# Run tests
if [ ${skip_tests} -eq 0 ]; then
	echo "Running tests"
	cd ${ISSM_DIR}/test/NightlyRun
	rm matlab.log 2> /dev/null

	# Run tests, redirecting output to logfile and suppressing output to console
	export ISSM_DIR_WIN=$(cygpath -w "${ISSM_DIR}")
	${MATLAB_PATH}/bin/matlab -nodesktop -nosplash -nojvm -r "try, addpath ${ISSM_DIR_WIN}/bin ${ISSM_DIR_WIN}/lib; runme(${MATLAB_NROPTIONS}); exit; catch me,fprintf('%s',getReport(me)); exit; end" -logfile matlab.log &

	# Wait for MATLAB to exit
	#
	# TODO:
	# - Replace by adding -wait option to above calls to matlab?
	#
	sleep 5;
	echo "Waiting for MATLAB to exit"
	pid=$(ps -W | grep MATLAB | awk '{print $1}')
	echo '-----------------------------'
	echo "pid: ${pid}"
	echo '-----------------------------'

	# Time out after $max_time seconds because sometimes multiple MATLAB processes get locked in race condition
	timer=0
	max_time=7200
	while [[ $timer -lt $max_time && -n "${pid}" ]]; do
		pid=$(ps -W | grep MATLAB | awk '{print $1}')
		timer=$((timer + 1))
		sleep 1;
	done

	# Check if timer hit $max_time
	if [ $timer -eq $max_time ]; then
		echo "Testing timed out at ${timer} seconds"
		# Kill MATLAB processes
		pid=$(ps -W | grep MATLAB | awk '{print $1}')
		echo "${pid}" | xargs /bin/kill -f
		exit 1
	fi

	# Filter out Windows characters
	cat matlab.log | tr -cd '\11\12\40-\176' > matlab.log2 && mv matlab.log2 matlab.log

	# Check that MATLAB did not exit in error
	matlabExitedInError=`grep -c -E "Activation cannot proceed|Error in|Illegal|Invalid MEX-file|Warning: Name is nonexistent or not a directory" matlab.log`

	if [ ${matlabExitedInError} -ne 0 ]; then
		echo "----------MATLAB exited in error!----------"
		cat matlab.log
		echo "-----------End of matlab.log-----------"

		# Clean up execution directory
		rm -rf ${ISSM_DIR}/execution/*

		exit 1
	fi

	# Check that all tests passed
	sed -i "/FAILED TO establish the default connection to the WindowServer/d" matlab.log # First, need to remove WindowServer error message
	numTestsFailed=`grep -c -E "FAILED|ERROR" matlab.log`

	if [ ${numTestsFailed} -ne 0 ]; then
		echo "One or more tests FAILED"
		cat matlab.log
		exit 1
	else
		echo "All tests PASSED"
	fi
else
	echo "Skipping tests"
fi

# Create package
cd ${ISSM_DIR}
git clean -d -f test # Clean up test directory (before copying to package)
echo "Copying assets to package: ${PKG}"
# NOTE: We do not copy lib directory to package as MATLAB seems to perform differently under Windows and so we package all DDL's and MEX-files in bin directory
cp -rf bin examples scripts share test ${PKG}
mkdir ${PKG}/execution

echo "Cleaning up unneeded/unwanted files"
rm -f ${PKG}/bin/generic_static.* # Remove static versions of generic cluster classes
rm -rf ${PKG}/test/SandBox # Remove testing sandbox from package

# Compress package
echo "Compressing package"
tar -czf ${COMPRESSED_PKG} ${PKG}
