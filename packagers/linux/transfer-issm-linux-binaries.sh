#!/bin/bash

################################################################################
# Transfers ISSM distributable package for Linux to ISSM website.
#
# NOTE:
# - Assumes that the following constants are defined,
#
#		COMPRESSED_PKG
#
# See also:
# - packagers/linux/complete-issm-linux-binaries-matlab.sh
# - packagers/linux/complete-issm-linux-binaries-python-2.sh
# - packagers/linux/complete-issm-linux-binaries-python-3.sh
################################################################################

# Transfer package to ISSM Web site
echo "Transferring package to ISSM Web site"
scp -P 2222 -i ~/.ssh/debian_12_vm-to-ross ${COMPRESSED_PKG} jenkins@ross.ics.uci.edu:/var/www/html/${COMPRESSED_PKG} # TODO: Fix issue with connecting on port 22

if [ $? -ne 0 ]; then
	echo "Transfer failed! Verify connection then build this project again (with -t/--transferonly option to skip building and packaging)."
	exit 1
fi
