#!/bin/bash

################################################################################
# Packages and tests ISSM distributable package for Linux with MATLAB API.
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
# - packagers/linux/complete-issm-linux-binaries-matlab.sh
#
# TODO:
# - Make sure that all TPL licenses are copied to package
################################################################################

# Expand aliases within the context of this script
shopt -s expand_aliases

## Override certain aliases
#
alias grep=$(which grep)

## Constants
#
LIBGCC="/usr/lib/x86_64-linux-gnu/libgcc_s.so.1" # Important that this is the library itself
LIBGCC_DIST="${ISSM_DIR}/lib/libgcc_s.so.1" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBGFORTRAN="/usr/lib/x86_64-linux-gnu/libgfortran.so.5.0.0" # Important that this is the library itself
LIBGFORTRAN_DIST="${ISSM_DIR}/lib/libgfortran.so.5" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBGMT="${ISSM_DIR}/externalpackages/gmt/install/lib/libgmt.so.6.6.0" # Important that this is the library itself
LIBGMT_DIST="${ISSM_DIR}/lib/libgmt.so.6" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBPSL="${ISSM_DIR}/externalpackages/gmt/install/lib/libpostscriptlight.so.6.6.0" # Important that this is the library itself
LIBPSL_DIST="${ISSM_DIR}/lib/libpostscriptlight.so.6" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBQUADMATH="/usr/lib/x86_64-linux-gnu/libquadmath.so.0.0.0" # Important that this is the library itself
LIBQUADMATH_DIST="${ISSM_DIR}/lib/libquadmath.so.0" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBSTDCXX="/usr/lib/x86_64-linux-gnu/libstdc++.so.6.0.30" # Important that this is the library itself
LIBSTDCXX_DIST="${ISSM_DIR}/lib/libstdc++.so.6.0.30" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it

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

echo "Moving certain shared libraries to lib/"
cp ${LIBGCC} ${LIBGCC_DIST} 2>/dev/null
cp ${LIBGFORTRAN} ${LIBGFORTRAN_DIST} 2>/dev/null
cp ${LIBQUADMATH} ${LIBQUADMATH_DIST} 2>/dev/null
cp ${LIBGMT} ${LIBGMT_DIST} 2>/dev/null
cp ${LIBPSL} ${LIBPSL_DIST} 2>/dev/null
cp ${LIBSTDCXX} ${LIBSTDCXX_DIST} 2>/dev/null

echo "Moving MPICH binaries to bin/"
if [ -f ${ISSM_DIR}/externalpackages/petsc/install/bin/mpiexec ]; then
	cp ${ISSM_DIR}/externalpackages/petsc/install/bin/mpiexec .
	cp ${ISSM_DIR}/externalpackages/petsc/install/bin/hydra_pmi_proxy .
elif [ -f ${ISSM_DIR}/externalpackages/mpich/install/bin/mpiexec ]; then
	cp ${ISSM_DIR}/externalpackages/mpich/install/bin/mpiexec .
	cp ${ISSM_DIR}/externalpackages/mpich/install/bin/hydra_pmi_proxy .
else
	echo "MPICH not found"
	exit 1
fi

echo "Moving GDAL binaries to bin/"
if [ -f ${ISSM_DIR}/externalpackages/gdal/install/bin/gdal-config ]; then
	cp ${ISSM_DIR}/externalpackages/gdal/install/bin/gdalsrsinfo .
	cp ${ISSM_DIR}/externalpackages/gdal/install/bin/gdaltransform .
else
	echo "GDAL not found"
	exit 1
fi

echo "Moving GMT binaries to bin/"
if [ -f ${ISSM_DIR}/externalpackages/gmt/install/bin/gmt ]; then
	cp ${ISSM_DIR}/externalpackages/gmt/install/bin/gmt .
elif [ -f ${ISSM_EXT_STATIC_DIR}/gmt/install/bin/gmtselect ]; then
	cp ${ISSM_DIR}/externalpackages/gmt/install/bin/gmtselect .
else
	echo "GMT not found"
	exit 1
fi

echo "Moving Gmsh binaries to bin/"
if [ -f ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh ]; then
	cp ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh .
else
	echo "Gmsh not found"
	exit 1
fi

echo "Moving GSHHG assets to share/"
if [ -d ${ISSM_DIR}/externalpackages/gshhg/install ]; then
	mkdir ${ISSM_DIR}/share 2>/dev/null
	cp -R ${ISSM_DIR}/externalpackages/gshhg/install/*.nc ${ISSM_DIR}/share
else
	echo "GSHHG not found"
	exit 1
fi

echo "Moving PROJ assets to bin/ and share/"
if [ -d ${ISSM_DIR}/externalpackages/proj/install/share/proj ]; then
	cp ${ISSM_DIR}/externalpackages/proj/install/bin/projinfo .
	mkdir ${ISSM_DIR}/share 2>/dev/null
	cp -R ${ISSM_DIR}/externalpackages/proj/install/share/proj ${ISSM_DIR}/share
else
	echo "PROJ not found"
	exit 1
fi

# Run tests
if [ ${skip_tests} -eq 0 ]; then
	echo "Running tests"
	cd ${ISSM_DIR}/test/NightlyRun
	rm matlab.log 2>/dev/null

	# Run tests, redirecting output to logfile and suppressing output to console
	${MATLAB_PATH}/bin/matlab -nojvm -nosplash -nojvm -r "try, addpath ${ISSM_DIR}/bin ${ISSM_DIR}/lib ${ISSM_DIR}/share; runme(${MATLAB_NROPTIONS}); exit; catch me,fprintf('%s',getReport(me)); exit; end" &> matlab.log

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
cp -rf bin examples lib scripts share test ${PKG}
mkdir ${PKG}/execution
echo "Cleaning up unneeded/unwanted files"
rm -f ${PKG}/bin/generic_static.* # Remove static versions of generic cluster classes
rm -f ${PKG}/lib/*.a # Remove static libraries from package
rm -f ${PKG}/lib/*.la # Remove libtool libraries from package
rm -rf ${PKG}/test/SandBox # Remove testing sandbox from package

# Compress package
echo "Compressing package"
tar -czf ${COMPRESSED_PKG} ${PKG}
