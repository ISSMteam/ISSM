#!/bin/bash

################################################################################
# Packages and tests ISSM distributable package for macOS with Python 3 API.
#
# Options:
# -s/--skiptests		Skip testing during packaging Use if packaging fails 
#						for some reason but build is valid.
#
# NOTE: The following variables *must* be defined in order for this script to 
# work properly,
#
#		COMPRESSED_PKG
#		ISSM_DIR
#		PKG
#		PYTHON_NROPTIONS
#
# See also:
# - packagers/mac/<ARCH>/complete-issm-mac-<ARCH>-binaries-python-3.sh
#
# TODO:
# - Make sure that all TPL licenses are copied to package
################################################################################

# Expand aliases within the context of this script
shopt -s expand_aliases

## Override certain other aliases
#
alias cp=$(which cp)
alias grep=$(which grep)

## Constants
#
LIBGMT="${ISSM_DIR}/externalpackages/gmt/install/lib/libgmt.6.6.0.dylib" # Important that this is the library itself
LIBGMT_DIST="${ISSM_DIR}/lib/libgmt.6.dylib" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it
LIBPSL="${ISSM_DIR}/externalpackages/gmt/install/lib/libpostscriptlight.6.6.0.dylib" # Important that this is the library itself
LIBPSL_DIST="${ISSM_DIR}/lib/libpostscriptlight.6.dylib" # Important the file name matches the SONAME entry in the binaries and other shared libraries which link to it

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
		-s|--skiptests)	skip_tests=1;					;;
		*) echo "Unknown parameter passed: $1"; exit 1	;;
	esac
fi

# Clean up from previous packaging
echo "Cleaning up existing assets"
cd ${ISSM_DIR}
rm -rf ${PKG} ${COMPRESSED_PKG}
mkdir ${PKG}

# Add required binaries and libraries to package and modify them where needed
cd ${ISSM_DIR}/bin

echo "Modify generic"
cat generic_static.py | sed -e "s/generic_static/generic/g" > generic.py

echo "Moving certain shared libraries to lib/"
cp ${LIBGMT} ${LIBGMT_DIST} 2>/dev/null
cp ${LIBPSL} ${LIBPSL_DIST} 2>/dev/null

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

echo "Moving Gmsh binaries to bin/"
if [ -f ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh ]; then
	cp ${ISSM_DIR}/externalpackages/gmsh/install/bin/gmsh .
else
	echo "Gmsh not found"
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
	rm python.log 2>/dev/null

	# Set Python environment
	export PYTHONPATH="${ISSM_DIR}/src/m/dev"
	export PYTHONSTARTUP="${PYTHONPATH}/devpath.py"
	export PYTHONUNBUFFERED=1 # We don't want Python to buffer output, otherwise issm.exe output is not captured

	# Run tests, redirecting output to logfile and suppressing output to console
	./runme.py ${PYTHON_NROPTIONS} &> python.log 2>&1

	# Check that Python did not exit in error
	pythonExitCode=`echo $?`
	pythonExitedInError=`grep -c -E "Error|No such file or directory|Permission denied|Standard exception|Traceback|bad interpreter|syntax error|error:" python.log`

	if [[ ${pythonExitCode} -ne 0 || ${pythonExitedInError} -ne 0 ]]; then
		echo "----------Python exited in error!----------"
		cat python.log
		echo "-----------End of python.log-----------"

		# Clean up execution directory
		rm -rf ${ISSM_DIR}/execution/*

		exit 1
	fi

	# Check that all tests passed
	sed -i '' "/FAILED TO establish the default connection to the WindowServer/d" python.log # First, need to remove WindowServer error message
	numTestsFailed=`grep -c -E "FAILED|ERROR" python.log`

	if [ ${numTestsFailed} -ne 0 ]; then
		echo "One or more tests FAILED"
		cat python.log
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
cp packagers/mac/issm-executable_entitlements.plist ${PKG}/bin/entitlements.plist
# ${ISSM_DIR}/scripts/py_to_pyc.sh ${PKG}/bin # Compile Python source files
echo "Cleaning up unneeded/unwanted files"
# rm -f ${PKG}/bin/*.py # Remove all Python scripts
rm -f ${PKG}/bin/generic_static.* # Remove static versions of generic cluster classes
rm -f ${PKG}/lib/*.a # Remove static libraries from package
rm -f ${PKG}/lib/*.la # Remove libtool libraries from package
rm -rf ${PKG}/test/SandBox # Remove testing sandbox from package

# Compress package
echo "Compressing package"
ditto -ck --sequesterRsrc --keepParent ${PKG} ${COMPRESSED_PKG}
