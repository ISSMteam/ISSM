#!/bin/bash

################################################################################
# Wrapper script to build, package, and transfer to ISSM Web site ISSM 
# distributable package for Windows with MATLAB API.
#
# Normally, we would put this directly into the project configuration under 
# 'Build' -> 'Execute shell', but because it is a bit more involved, it is a 
# good idea to version it.
#
# When no failures/errors occur, performs the following:
# - Builds ISSM according to configuration.
# - Packages executables and libraries.
# - Runs test suite against package.
# - Transmits package to ISSM Web site for distribution.
#
# Options:
# -b/--skipbuild		Skip ISSM compilation.
# -s/--skiptests		Skip ISSM compilation and testing during packaging 
#						step. Use if packaging fails for some reason but build 
#						is valid.
# -t/--transferonly		Transfer package to ISSM Web site only. Use if transfer 
#						fails for some reason to skip building, packaging, and 
#						signing.
#
# NOTE:
# - Use only *one* of the above options at a time, and make sure it is removed 
#	again after a single run.
# - Builds will fail when any of the above options are used on a clean 
#	workspace. For example, if 'Source Code Management' -> 'Check-out Strategy' 
#	select menu is set to "Always check out a fresh copy".
################################################################################

## Constants
#
PKG="ISSM-Windows-MATLAB" # Name of directory to copy distributable files to

COMPRESSED_PKG="${PKG}.tar.gz"

## Environment
#
export COMPRESSED_PKG
export PKG

## Parse options
#
if [ $# -gt 1 ]; then
	echo "Can use only one option at a time"
	exit 1
fi

# NOTE: We could do this with binary switching (i.e. 0011 to sign and transfer, 
#		but the following is self-documenting).
#
build=1
package=1
transfer=1

if [ $# -eq 1 ]; then
	case $1 in
		-b|--skipbuild)		build=0;				shift	;;
		-s|--skiptests)		build=0;						;;
		-t|--transferonly)	build=0;	package=0;			;;
		*) echo "Unknown parameter passed: $1"; exit 1 		;;
	esac
fi

# Build
if [ ${build} -eq 1 ]; then
	./jenkins/jenkins.sh ./jenkins/ross-win-msys2-mingw-msmpi-binaries-matlab

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

# Package
if [ ${package} -eq 1 ]; then
	./packagers/win/package-issm-win-binaries-matlab.sh $1

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

# Transfer distributable package to ISSM Web site
if [ ${transfer} -eq 1 ]; then
	./packagers/win/transfer-issm-win-binaries.sh

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

