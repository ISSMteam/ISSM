#!/bin/bash

################################################################################
# Intended to be run in the context of a Jenkins project on a JPL 
# Cybersecurity server for signing macOS applications. Polls SCM of the 
# Subversion repository hosted at 
# https://issm.ess.uci.edu/svn/issm-macos-signing/<ARCH>/<API>/unsigned to 
# trigger new builds.
#
# In order to replicate the required Jenkins project configuration,
# - first, navigate to 'Manage Jenkins' -> 'Manage Plugins' and install the 
#	'Credentials Bindings Plugin' if it is not already installed.
# - contact one of the members of the ISSM development team for credentials for 
#	the ISSM binaries repository (mention that the credentials are stored in 
#	ISSM-Infrastructure.pdf).
# - navigate to 'Manage Jenkins' -> 'Manage Credentials' -> <domain> -> 
#	'Add Credentials' and enter the credentials from above.
# - from the 'Dashboard', select 'New Item' -> 'Freestyle project'.
# - under 'Source Code Management', select 'Subversion',
#		- the 'Repository URL' text field should be set to 
#
# 		https://issm.ess.uci.edu/svn/issm-binaries/mac/<ARCH>/<API>/unsigned
#
#		where,
#
#		<ARCH>			'intel' or 'silicon'
#		<API>			'matlab' or 'python'
#
#		- the 'Credentials' select menu should be set to the new credentials 
#		created previously.
#		- the 'Local module directory' text field should be set to the same 
#		value as the constant UNSIGNED_REPO_COPY (set below to './unsigned').
# - under 'Build Triggers', check the box for 'Poll SCM' and set the 
#	'Schedule' text area to "H/5 * * * *".
# - under 'Build Environment', check the box for 'Use secret text(s) or 
#	file(s)', then under 'Bindings' click the 'Add...' button and select 
#	'Username and password (separated)',
#		- set 'Username Variable' to "ISSM_BINARIES_USER".
#		- set 'Password Variable' to "ISSM_BINARIES_PASS".
# - under 'Credentials', select the same, new credentials that created 
#	previously.
# - the contents of this script can be copied/pasted directly into the ‘Build' 
#	-> 'Execute Shell' -> ‘Command' textarea of the project configuration (or 
#	you can simply store the script on disk and call it from there).
# - make sure to click the 'Save' button.
#
# Current point of contact at JPL Cybersecurity:
#	Alex Coward, alexander.g.coward@jpl.nasa.gov
#
# NOTE:
# - Assumes that "ISSM_BINARIES_USER" and "ISSM_BINARIES_PASS" are set up in 
#	the 'Bindings' section under a 'Username and password (separated)' binding 
#	(requires 'Credentials Binding Plugin').
# - For local debugging, the aforementioned credentials can be hardcoded into 
#	the 'USERNAME' and 'PASSWORD' constants below.
################################################################################

# Expand aliases within the context of this script
shopt -s expand_aliases

## Override certain other aliases
#
alias cp=$(which cp)
alias grep=$(which grep)

## Constants
#

## NOTE: The following need to be set with the proper credentials
#
AD_IDENTITY="**********" # Apple Developer identity
AD_USERNAME="**********" # Apple Developer username
ALTOOL_PASSWORD="@keychain:**********" # altool password (assumed to be stored in keychain)
ASC_PROVIDER="**********"

## NOTE: The following need to be set for the particular signing job (see comments for options)
#
PKG="ISSM-macOS-<ARCH>-<API>" # <ARCH>: 'Intel' or 'Silicon'; <API>: 'MATLAB' or 'Python-3'
VARIANT_REPO_SUBPATH="<ARCH>/<API>" # <ARCH>: 'intel' or 'silicon'; <API>: 'matlab' or 'python'

EXE_ENTITLEMENTS_PLIST="${PKG}/bin/entitlements.plist"
MAX_SVN_ATTEMPTS=10
NOTARIZATION_CHECK_ATTEMPTS=20
NOTARIZATION_CHECK_PERIOD=60
NOTARIZATION_LOGFILE="notarization.log"
NOTARIZATION_LOGFILE_PATH="."
PASSWORD=${ISSM_BINARIES_PASS}
SIGNING_REPO_BASE_URL="https://issm.ess.uci.edu/svn/issm-macos-signing"
SIGNING_REPO_URL="${SIGNING_REPO_BASE_URL}/${VARIANT_REPO_SUBPATH}"
SIGNED_REPO_COPY="./signed"
SIGNING_LOCK_FILE="signing.lock"
UNSIGNED_REPO_COPY="./unsigned"
USERNAME=${ISSM_USERNAME}

COMPRESSED_PKG="${PKG}.zip"

SIGNED_REPO_URL="${SIGNING_REPO_URL}/signed"
UNSIGNED_REPO_URL="${SIGNING_REPO_URL}/unsigned"

# NOTE: Uncomment the following for local testing (Jenkins checks out copy of 
#		repository for unsigned packages to working directory)
#

# # Clean up from previous packaging (not necessary for single builds on Jenkins, 
# # but useful when testing packaging locally)
# echo "Cleaning up existing assets"
# rm -rf ${COMPRESSED_PKG} ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE} ${UNSIGNED_REPO_COPY}

# # Check out copy of repository for unsigned packages
# echo "Checking out copy of repository for unsigned packages"
# svn checkout \
# 	--trust-server-cert \
# 	--non-interactive \
# 	--username ${USERNAME} \
# 	--password ${PASSWORD} \
# 	${UNSIGNED_REPO_URL} \
# 	${UNSIGNED_REPO_COPY}

rm -rf ${PKG} ${SIGNED_REPO_COPY}

# Extract package contents
echo "Extracting package contents"
ditto -xk ${UNSIGNED_REPO_COPY}/${COMPRESSED_PKG} .

# Clear extended attributes on all files
xattr -cr ${PKG}

# Build list of ISSM executables, libraries, and packages
ISSM_BINS=$(\
	find ${PKG}/lib -type f -name *.dylib; \
	find ${PKG}/bin -type f -name *.exe; \
	find ${PKG}/lib -type f -name *.mexmaca64; \
	find ${PKG}/lib -type f -name *.mexmaci64; \
	find ${PKG}/test -type f -name *.pkg; \
	find ${PKG}/bin -type f -name *.pyc; \
	find ${PKG}/lib -type f -name *.so; \
)

# Build list of third party executables
THIRD_PARTY_BINS=$(\
	echo ${PKG}/bin/gdalsrsinfo; \
	echo ${PKG}/bin/gdaltransform; \
	echo ${PKG}/bin/gmsh; \
	echo ${PKG}/bin/gmt; \
	echo ${PKG}/bin/hydra_pmi_proxy; \
	echo ${PKG}/bin/mpiexec; \
	echo ${PKG}/bin/projinfo; \
)

# Sign all executables in package
echo "Signing ISSM executables, libraries, and packges in package"
codesign -s ${AD_IDENTITY} -v --timestamp --options=runtime --entitlements ${EXE_ENTITLEMENTS_PLIST} ${ISSM_BINS}
for ISSM_BIN in "${ISSM_BINS[@]}"; do
	codesign -s ${AD_IDENTITY} -v --timestamp --options=runtime ${ISSM_BIN}
done

echo "Signing third-party executables in package"
for THIRD_PARTY_BIN in "${THIRD_PARTY_BINS[@]}"; do
	codesign -s ${AD_IDENTITY} -v --timestamp --options=runtime ${THIRD_PARTY_BIN}
done

# Validate timestamp
echo "Validating timestamp on an ISSM executable"
codesign -dvv ${PKG}/bin/issm.exe
echo "Validating timestamp on a third-party executable"
codesign -dvv ${PKG}/bin/hydra_pmi_proxy

# NOTE: Skipping signature validation because this is not a true package nor app

# Compress signed package
echo "Compressing signed package"
ditto -ck --sequesterRsrc --keepParent ${PKG} ${COMPRESSED_PKG}

# Submit compressed package for notarization
echo "Submitting signed package to Apple for notarization"
xcrun notarytool submit ${COMPRESSED_PKG} --apple-id "$AD_USERNAME" --team-id "$TEAM_ID" --password "$NOTARY_PASSWORD" --wait &> ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE}

echo "Notarization request response received"

# Check if UUID exists in response
HAS_UUID=$(grep 'id: ' ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE})
if [ -z "${HAS_UUID}" ]; then
	echo "Notarization failed!"
	echo "----------------------- Contents of notarization logfile -----------------------"
	cat ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE}
	echo "--------------------------------------------------------------------------------"

	# Clean up
	rm -rf ${PKG} ${COMPRESSED_PKG}

	exit 1
fi

# Get UUID from notarization request response
UUID=$(echo ${HAS_UUID} | sed 's/[[:space:]]*id: //')
echo "UUID: ${UUID}" 

# Check notarization status
#
# NOTE: Currently, this checks if notarization was successful, but we are not 
#		able to staple notarization as this is not a true package nor app and, 
#		at the very least, MATLAB Mex files cannot be stapled. As such, clients 
#		will not be able to clear Gatekeeper if they are offline.
#
echo "Checking notarization status"
SUCCESS=0
xcrun notarytool log ${UUID} --apple-id "$AD_USERNAME" --team-id "$TEAM_ID" --password "$NOTARY_PASSWORD" &> ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE}
STATUS=$(grep '"status": "Accepted"' ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE} | wc -l)

if [[ ${STATUS} -gt 0 ]]; then
	# Staple notarization to all elements of package that were previously signed
	#xcrun stapler staple ${THIRD_PARTY_BINS} # NOTE: Fails with "Stapler is incapable of working with MATLAB Mex files."

	# Validate stapling of notarization
	#xcrun stapler validation ${THIRD_PARTY_BINS} # NOTE: Skipping notarization stapling validation because this is not a true package nor app

	# Compress signed and notarized package
	ditto -ck --sequesterRsrc --keepParent ${PKG} ${COMPRESSED_PKG}

	echo "Notarization successful!"

	# Set flag indicating notarization was successful
	SUCCESS=1
else
	echo "Notarization failed!"
	echo "----------------------- Contents of notarization logfile -----------------------"
	cat ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE}
	echo "--------------------------------------------------------------------------------"
fi

# Check out copy of repository for signed packages
echo "Checking out copy of repository for signed packages"
SVN_ATTEMPT=0
SVN_SUCCESS=0
while [[ ${SVN_ATTEMPT} -lt ${MAX_SVN_ATTEMPTS} && ${SVN_SUCCESS} -eq 0 ]]; do
	rm -rf ${SIGNED_REPO_COPY}
	svn checkout \
		--trust-server-cert \
		--non-interactive \
		--username ${USERNAME} \
		--password ${PASSWORD} \
		${SIGNED_REPO_URL} \
		${SIGNED_REPO_COPY} > /dev/null 2>&1
	if [ $? -eq 0 ]; then
		SVN_SUCCESS=1
		break
	else
		((++SVN_ATTEMPT))
		sleep 5
	fi
done

if [ ${SVN_SUCCESS} -eq 0 ]; then
	echo "Checkout of repository for signed packages failed"
	exit 1
fi

# Copy notarization file to repository for signed packages
cp ${NOTARIZATION_LOGFILE_PATH}/${NOTARIZATION_LOGFILE} ${SIGNED_REPO_COPY}
svn add ${SIGNED_REPO_COPY}/${NOTARIZATION_LOGFILE} > /dev/null 2>&1

# Remove lock file from repository for signed packages
svn delete ${SIGNED_REPO_COPY}/${SIGNING_LOCK_FILE} > /dev/null 2>&1

SVN_ATTEMPT=0
SVN_SUCCESS=0
if [ ${SUCCESS} -eq 1 ]; then
	# Copy signed package to repository for signed packages
	cp ${COMPRESSED_PKG} ${SIGNED_REPO_COPY}
	svn add ${SIGNED_REPO_COPY}/${COMPRESSED_PKG} > /dev/null 2>&1

	# Commit changes
	echo "Committing changes to repository for signed packages"
	while [[ ${SVN_ATTEMPT} -lt ${MAX_SVN_ATTEMPTS} && ${SVN_SUCCESS} -eq 0 ]]; do
		svn commit \
			--trust-server-cert \
			--non-interactive \
			--username ${USERNAME} \
			--password ${PASSWORD} \
			--message "CHG: New signed package (success)" ${SIGNED_REPO_COPY} > /dev/null 2>&1
		if [ $? -eq 0 ]; then
			SVN_SUCCESS=1
			break
		else
			((++SVN_ATTEMPT))
			sleep 5
		fi
	done

	if [ ${SVN_SUCCESS} -eq 0 ]; then
		echo "Commit to repository for signed packages failed"
		exit 1
	fi
else
	# Commit changes
	echo "Committing changes to repository for signed packages"
	while [[ ${SVN_ATTEMPT} -lt ${MAX_SVN_ATTEMPTS} && ${SVN_SUCCESS} -eq 0 ]]; do
		svn commit \
			--trust-server-cert \
			--non-interactive \
			--username ${USERNAME} \
			--password ${PASSWORD} \
			--message "CHG: New signed package (failure)" ${SIGNED_REPO_COPY} > /dev/null 2>&1
		if [ $? -eq 0 ]; then
			SVN_SUCCESS=1
			break
		else
			((++SVN_ATTEMPT))
			sleep 5
		fi
	done

	if [ ${SVN_SUCCESS} -eq 0 ]; then
		echo "Commit to repository for signed packages failed"
		exit 1
	fi

	exit 1
fi
