#!/usr/bin/env python

# NOTE: This script is a stripped-down version of
#       $ISSM_DIR/src/m/dev/devpath.py and is intended only for loading ISSM in 
#       order to test our distributable packages. It assumes the following is 
#       set before $ISSM_DIR/test/NightlyRun/runme.py is called, 
#
#           export ISSM_DIR=</path/to/ISSM>
#           export PATH="${PATH}:${ISSM_DIR}/bin:${ISSM_DIR}/scripts"
#           export PYTHONPATH="${ISSM_DIR}/scripts"
#           export PYTHONSTARTUP="${PYTHONPATH}/devpath.py"
#           export PYTHONUNBUFFERED=1
#

import os
import sys

ISSM_DIR = os.getenv('ISSM_DIR')
sys.path.append(ISSM_DIR + '/bin')
sys.path.append(ISSM_DIR + '/lib')

from issmversion import issmversion
