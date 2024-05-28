#!/bin/bash

################################################################################
# Transfers ISSM distributable package for Windows to ISSM website.
#
# NOTE:
# - Assumes that the following constants are defined,
#
#		COMPRESSED_PKG
#
# See also:
# - packagers/win/complete-issm-win-binaries-matlab.sh
# - packagers/win/complete-issm-win-binaries-python-2.sh
# - packagers/win/complete-issm-win-binaries-python-3.sh
################################################################################

# Transfer package to ISSM Web site
echo "Transferring package to ISSM Web site"
scp -i ~/.ssh/windows_10-vm_to_ross ${COMPRESSED_PKG} jenkins@ross.ics.uci.edu:/var/www/html/${COMPRESSED_PKG}

if [ $? -ne 0 ]; then
	echo "Transfer failed! Verify connection then build this project again (with -t/--transferonly option to skip building and packaging)."
	exit 1
fi
