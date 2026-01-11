#!/bin/bash

################################################################################
# Wrapper script to build, package, send for signing, and transfer to ISSM 
# website ISSM distributable package with MATLAB API for macOS running on 
# Silicon.
#
# Normally, we would put this directly into the project configuration under 
# 'Build' -> 'Execute shell', but because it is a bit more involved, it is a 
# good idea to version it.
#
# When no failures/errors occur, performs the following:
# - Builds ISSM according to configuration.
# - Packages executables and libraries.
# - Runs test suite against package.
# - Commits compressed package to repository to be signed by JPL Cybersecurity.
# - Retrieves signed package and transmits it to ISSM Web site for 
#	distribution.
#
# Options:
# -b/--skipbuild		Skip ISSM compilation.
# -r/--resign			Skip ISSM compilation and packaging. Use to retrigger 
#						signing/notarization if it fails but build and package 
#						are valid.
# -s/--skiptests		Skip ISSM compilation and testing during packaging 
#						step. Use if packaging fails for some reason but build 
#						is valid.
# -t/--transferonly		Transfer package to ISSM Web site only. Use if transfer 
#						fails for some reason to skip building, packaging, and 
#						signing.
# -u/--unlock			Remove lock file from signed package repository. Use if 
#						build is aborted to allow for subsequent fresh build.
#
# Debugging:
# - Relies on a very tight handshake with project on remote JPL Cybersecurity 
#	Jenkins server. Debugging may be performed locally by running,
#
#		packagers/mac/sign-issm-mac-binaries.sh
#
#	with "AD_IDENTITY", "AD_USERNAME", and "ASC_PROVIDER" hardcoded to Apple 
#	Developer credentials (make sure to also set keychain password in 
#	"ALTOOL_PASSWORD") and "PKG" and "VARIANT_REPO_SUBPATH" properly set for 
#	the package signing to be tested (should match constants "PKG" and 
#	"REPO_BASE_URL" defined in this script).
# - Removing stdout/stderr redirections to null device (> /dev/null 2>&1) can 
#	help debug potential SVN issues.
#
# NOTE:
# - Use only *one* of the above options at a time, and make sure it is removed 
#	again after a single run.
# - Builds will fail when any of the above options are used on a clean 
#	workspace. For example, if 'Source Code Management' -> 'Check-out Strategy' 
#	select menu is set to "Always check out a fresh copy".
# - Assumes that "ISSM_BINARIES_USER" and "ISSM_BINARIES_PASS" are set up in 
#	the 'Bindings' section under a 'Username and password (separated)' binding 
#	(requires 'Credentials Binding Plugin') with 'Credentials' select menu set 
#	to "jenkins/****** (SVN repository for ISSM binaries)".
#
# TODO:
# - Generalize these wrapper scripts even further by defining only constants in 
# config file (the rest of the script is identical).
################################################################################

## Constants
#
PKG="ISSM-macOS-Silicon-MATLAB" # Name of directory to copy distributable files to
VARIANT_REPO_SUBPATH="silicon/matlab"

MATLAB_NROPTIONS="'benchmark','all','exclude',[119,124:126,129,216,234:235,274,362,417:418,420,423,430,433,435,441:442,444:445,448,456,462:464,470:476,508,517,544,546,701:703,808,1101:1110,1201:1208,1301:1304,1401:1402,1601:1602,2002,2004,2006,2010:2013,2020:2021,2052:2053,2085,2090:2092,2110:2113,2424:2425,3001:3300,3480:3481,4001:4100]" # NOTE: Combination of test suites from basic, Dakota, and Solid Earth builds, with tests that require a restart and those that require the JVM excluded
MATLAB_PATH="/Applications/MATLAB_R2023b.app"
SIGNING_REPO_BASE_URL="https://issm.ess.uci.edu/svn/issm-macos-signing"
SIGNING_REPO_URL="${SIGNING_REPO_BASE_URL}/${VARIANT_REPO_SUBPATH}"
SIGNED_REPO_COPY="./signed"
UNSIGNED_REPO_COPY="./unsigned"

COMPRESSED_PKG="${PKG}.zip"
SIGNED_REPO_URL="${SIGNING_REPO_URL}/signed"
UNSIGNED_REPO_URL="${SIGNING_REPO_URL}/unsigned"

## Environment
#
export COMPRESSED_PKG
export MATLAB_NROPTIONS
export MATLAB_PATH
export PKG
export SIGNED_REPO_COPY
export SIGNED_REPO_URL
export UNSIGNED_REPO_COPY
export UNSIGNED_REPO_URL

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
sign=1
transfer=1

if [ $# -eq 1 ]; then
	case $1 in
		-b|--skipbuild)		build=0;							shift	;;
		-r|--resign)		build=0;	package=0;						;;
		-s|--skiptests)		build=0;									;;
		-t|--transferonly)	build=0;	package=0;	sign=0;				;;
		-u|--unlock)		build=0;	package=0;	transfer=0;			;;
		*) echo "Unknown parameter passed: $1"; exit 1 					;;
	esac
fi

# Build
if [ ${build} -eq 1 ]; then
	./jenkins/jenkins.sh ./jenkins/mac-silicon-binaries-matlab

	if [ $? -ne 0 ]; then
		echo "Failure while compiling"
		exit 1
	fi
fi

# Package
if [ ${package} -eq 1 ]; then
	./packagers/mac/package-issm-mac-binaries-matlab.sh $1

	if [ $? -ne 0 ]; then
		echo "Failure during packaging"
		exit 1
	fi

	shift # Clear $1 so that it is not passed to commit_for_signing script
fi

# Commit for signing
if [ ${sign} -eq 1 ]; then
	./packagers/mac/commit_for_signing-issm-mac-binaries.sh $1

	if [ $? -ne 0 ]; then
		echo "Failure while committing package for signing"
		exit 1
	fi
fi

# NOTE: Because Mac build nodes are no longer directly connected to UCI 
#		network and because remote access requires a VPN connection, we can 
#		no longer transfer signed distributables via SSH. For now, there is 
#		a cron job running every five minutes under user jenkins on 
#		ross.ics.uci.edu that runs 
#		/home/jenkins/bin/update-issm-mac-binaries.sh, which checks for 
#		updated, signed distributables in the ISSM Binaries SVN repository
#		and if they are available, copies them to the public directory.
#

# # Transfer distributable package to ISSM Web site
# if [ ${transfer} -eq 1 ]; then
# 	./packagers/mac/transfer-issm-mac-binaries.sh

# 	if [ $? -ne 0 ]; then 
# 		exit 1
# 	fi
# fi

