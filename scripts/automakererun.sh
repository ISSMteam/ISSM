#!/bin/sh
#
# automakererun.sh
# Remakes GNU Build System files.
#
# Currently, a wrapper for Autotools' autoreconf script.
#
# Run,
#	autoreconf --help
# for more information.

cd $ISSM_DIR
autoreconf -ivf
