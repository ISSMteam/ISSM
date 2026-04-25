#!/bin/bash

################################################################################
# Wrapper script to build, package, and transfer to ISSM Web site ISSM 
# distributable package for Linux with Python 3 API.
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
PKG="ISSM-Linux-Python-3" # Name of directory to copy distributable files to
PYTHON_NROPTIONS="--benchmark all --exclude 125:126 129 234:235 417:418 420 435 444:445 456 550 701:703 1101:1110 1201:1208 1301:1304 1401:1402 1601:1602 2002 2004 2006 2010:2013 2020:2021 2052:2053 2090:2092 2110:2113 2424:2425 3001:3300 3480:3481 4001:4100" # NOTE: Combination of test suites from basic, Dakota, and Solid Earth builds, with tests that require a restart and those that require the JVM excluded

COMPRESSED_PKG="${PKG}.tar.gz"

## Environment
#
export COMPRESSED_PKG
export PKG
export PYTHON_NROPTIONS

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
	./jenkins/jenkins.sh ./jenkins/ross-debian_linux-binaries-python-3

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

# Package
if [ ${package} -eq 1 ]; then
	./packagers/linux/package-issm-linux-binaries-python-3.sh $1

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

# Transfer distributable package to ISSM Web site
if [ ${transfer} -eq 1 ]; then
	./packagers/linux/transfer-issm-linux-binaries.sh

	if [ $? -ne 0 ]; then 
		exit 1
	fi
fi

