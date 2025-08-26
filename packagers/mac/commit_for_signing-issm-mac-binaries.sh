#!/bin/bash

################################################################################
# Commits ISSM distributable package with for macOS running on to repository 
# for signing. This repository is polled by a project running on a JPL 
# Cybersecurity Jenkins server and performs the actual signing and 
# notarization.
#
# NOTE: The following variables *must* be defined in order for this script to 
# work properly,
#
#		COMPRESSED_PKG
#		ISSM_BINARIES_REPO_PASS
#		ISSM_BINARIES_REPO_USER
#		SIGNED_REPO_COPY
#		SIGNED_REPO_URL
#		UNSIGNED_REPO_COPY
#		UNSIGNED_REPO_URL
#
# Options:
# -r/--resign			Skip ISSM compilation and packaging. Use to retrigger 
#						signing/notarization if it fails but build and package 
#						are valid.
# -u/--unlock			Remove lock file from signed package repository. Use if 
#						build is aborted to allow for subsequent fresh build.
#
# See also:
# - packagers/mac/<ARCH>/complete-issm-mac-<ARCH>-binaries-<API>.sh
#
# TODO:
# - Generalize checkout_*_repo_copy and validate_*_repo_copy functions (e.g. 
#	pass 'signed' or 'unsigned' as argument)
################################################################################

# Expand aliases within the context of this script
shopt -s expand_aliases

# NOTE: For some reason, calling svn from within the context of this script 
#		gives,
#
#			svn: command not found
#
#		even though it is installed via Homebrew and available at the following 
#		path.
#
alias svn=$(which svn)

## Override certain other aliases
#
alias cp=$(which cp)
alias grep=$(which grep)

## Constants
#
MAX_SIGNING_CHECK_ATTEMPTS=45
NOTARIZATION_LOGFILE="notarization.log"
RETRIGGER_SIGNING_FILE="retrigger.txt"
SIGNING_CHECK_PERIOD=60 # in seconds
SIGNING_LOCK_FILE="signing.lock"

## Functions
#
checkout_signed_repo_copy(){
	echo "Checking out copy of repository for signed packages"

	# NOTE: Get empty copy because we do not want to have to check out package 
	#		from previous signing.
	#
	svn checkout \
		--trust-server-cert \
		--non-interactive \
		--depth empty \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		${SIGNED_REPO_URL} \
		${SIGNED_REPO_COPY} > /dev/null 2>&1
}
checkout_unsigned_repo_copy(){
	echo "Checking out copy of repository for unsigned packages"
	svn checkout \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		${UNSIGNED_REPO_URL} \
		${UNSIGNED_REPO_COPY} > /dev/null 2>&1
}
validate_signed_repo_copy(){
	# Validate copy of repository for signed binaries (e.g. 
	# 'Check-out Strategy' was set to 'Use 'svn update' as much as possible'; 
	# initial checkout failed)
	if [[ ! -d ${SIGNED_REPO_COPY} || ! -d ${SIGNED_REPO_COPY}/.svn ]]; then
		rm -rf ${SIGNED_REPO_COPY}
		checkout_signed_repo_copy
	fi
}
validate_unsigned_repo_copy(){
	# Validate copy of repository for unsigned binaries (e.g. 
	# 'Check-out Strategy' was set to 'Use 'svn update' as much as possible'; 
	# initial checkout failed)
	if [[ ! -d ${UNSIGNED_REPO_COPY} || ! -d ${UNSIGNED_REPO_COPY}/.svn ]]; then
		rm -rf ${UNSIGNED_REPO_COPY}
		checkout_unsigned_repo_copy
	fi
}

## Parse options
#
if [ $# -gt 1 ]; then
	echo "Can use only one option at a time"
	exit 1
fi

retrigger_signing=0
unlock=0

if [ $# -eq 1 ]; then
	case $1 in
		-r|--resign)	retrigger_signing=1;	;;
		-u|--unlock)	unlock=1;				;;
		*) echo "Unknown parameter passed: $1"; exit 1	;;
	esac
fi

validate_signed_repo_copy

if [ ${unlock} -eq 1 ]; then
	# Remove signing lock file from signed package repository so that a new 
	# build can run
	echo "Removing lock file from repository for signed packages"
	svn update \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} > /dev/null 2>&1
	svn delete ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} > /dev/null 2>&1
	svn commit \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		--message "DEL: Removing lock file after failed build" ${SIGNED_REPO_COPY} > /dev/null 2>&1
	svn cleanup ${SIGNED_REPO_COPY} > /dev/null 2>&1

	echo "Remove -u/--unlock option from configuration and run again"
	exit 1
fi

# If lock file exists, a signing build is still in process by JPL Cybersecurity
svn update \
	--trust-server-cert \
	--non-interactive \
	--username ${ISSM_BINARIES_REPO_USER} \
	--password ${ISSM_BINARIES_REPO_PASS} \
	${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} > /dev/null 2>&1

if [ -f ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} ]; then
	echo "Previous signing job still in process by JPL Cybersecurity. Please try again later."
	exit 1
fi

# Commit lock file to repository for signed packages
echo "Committing lock file to repository for signed packages"
touch ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE}
svn add ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} > /dev/null 2>&1
svn commit \
	--trust-server-cert \
	--non-interactive \
	--username ${ISSM_BINARIES_REPO_USER} \
	--password ${ISSM_BINARIES_REPO_PASS} \
	--message "ADD: New lock file" ${SIGNED_REPO_COPY} > /dev/null 2>&1

# Check out copy of repository for unsigned packages
validate_unsigned_repo_copy

if [ ${retrigger_signing} -eq 0 ]; then
	# Commit new compressed package to repository for unsigned binaries
	echo "Committing package to repository for unsigned packages"
	cp ${COMPRESSED_PKG} ${UNSIGNED_REPO_COPY}
	svn add ${UNSIGNED_REPO_COPY}/${COMPRESSED_PKG} > /dev/null 2>&1
	svn commit \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		--message "CHG: New unsigned package" ${UNSIGNED_REPO_COPY} > /dev/null 2>&1
else
	# NOTE: If notarize_only == 1, we commit a dummy file as we do not want to 
	#		have to commit the entire compressed package again simply to 
	#		retrigger the signing build on the remote JPL Cybersecurity Jenkins 
	#		server.
	#
	echo "Attempting to sign existing package again"
	echo $(date +'%Y-%m-%d-%H-%M-%S') > ${UNSIGNED_REPO_COPY}/${RETRIGGER_SIGNING_FILE} # Write datetime stamp to file to ensure modification is made
	svn add ${UNSIGNED_REPO_COPY}/${RETRIGGER_SIGNING_FILE} > /dev/null 2>&1
	svn commit \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		--message "ADD: Retriggering signing with same package (previous attempt failed)" ${UNSIGNED_REPO_COPY} > /dev/null 2>&1
fi

# Check status of signing
echo "Checking progress of signing..."
SIGNING_CHECK_ATTEMPT=0
while [ ${SIGNING_CHECK_ATTEMPT} -lt ${MAX_SIGNING_CHECK_ATTEMPTS} ]; do
	echo "...in progress still; checking again in ${SIGNING_CHECK_PERIOD} seconds"
	sleep ${SIGNING_CHECK_PERIOD}
	svn update \
		--trust-server-cert \
		--non-interactive \
		--username ${ISSM_BINARIES_REPO_USER} \
		--password ${ISSM_BINARIES_REPO_PASS} \
		${SIGNED_REPO_COPY} > /dev/null 2>&1

	if [ ! -f ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} ]; then
		# Retrieve notarization lock file
		svn update \
			--trust-server-cert \
			--non-interactive \
			--username ${ISSM_BINARIES_REPO_USER} \
			--password ${ISSM_BINARIES_REPO_PASS} \
			${SIGNED_REPO_COPY}/${NOTARIZATION_LOGFILE}

		# Check status
		STATUS=$(grep '"status": "Accepted"' ${SIGNED_REPO_COPY}/${NOTARIZATION_LOGFILE} | wc -l)

		if [[ ${STATUS} -gt 0 ]]; then
			echo "Notarization successful!"
			break
		else
			echo "Notarization failed!"
			echo "----------------------- Contents of notarization logfile -----------------------"
			cat ${SIGNED_REPO_COPY}/${NOTARIZATION_LOGFILE}
			echo "--------------------------------------------------------------------------------"

			exit 1
		fi
	else
		((++SIGNING_CHECK_ATTEMPT))
	fi
done

if [ ! -f ${SIGNED_REPO_COPY}/${NOTARIZATION_LOGFILE} ]; then
	echo "Signing timed out!"
	exit 1
fi
