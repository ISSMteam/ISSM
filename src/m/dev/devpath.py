#!/usr/bin/env python3
import os
import sys

# Recover ISSM_DIR and USERNAME
ISSM_DIR = os.getenv('ISSM_DIR')
USERNAME = os.getenv('USER')
JPL_SVN = os.getenv('JPL_SVN')
if ISSM_DIR is None:
    raise NameError('"ISSM_DIR" environment variable is empty! You should define ISSM_DIR in your .cshrc or .bashrc!')

# Go through src/m and append any directory that contains a *.py file to PATH
for root, dirs, files in os.walk(ISSM_DIR + '/src/m'):
    if '.svn' in dirs:
        dirs.remove('.svn')
    for file in files:
        if file.find('.py') != -1:
            if file.find('.pyc') == -1:
                if root not in sys.path:
                    sys.path.append(root)

# Also add the Nightly run directory
if ISSM_DIR + '/test/NightlyRun' not in sys.path:
    sys.path.append(ISSM_DIR + '/test/NightlyRun')
if ISSM_DIR + '/lib' not in sys.path:
    sys.path.append(ISSM_DIR + '/lib')
if ISSM_DIR + '/lib-precompiled' not in sys.path:
    sys.path.append(ISSM_DIR + '/lib-precompiled') # load precompiled MEX files; remove after MEX file compilation is supported on Silicon-based Macs
if ISSM_DIR + '/src/wrappers/python/.libs' not in sys.path:
    sys.path.append(ISSM_DIR + '/src/wrappers/python/.libs')

# If using clusters, we need to have the path to the cluster settings directory
if JPL_SVN is not None:
    jpl_path = JPL_SVN + '/usr/' + USERNAME
    if os.path.exists(jpl_path):
        if jpl_path not in sys.path:
            sys.path.append(jpl_path)
    else:
        print('Warning: devpath.py: cluster settings should be in {}'.format(jpl_path))

try: # Avoid circular import
    from runme import runme  # First, because plotmodel may fail
except:
    pass
from plotmodel import plotmodel

#c = get_ipython().config
#c.InteractiveShellApp.exec_lines = []
#c.InteractiveShellApp.exec_lines.append('%load_ext autoreload')
#c.InteractiveShellApp.exec_lines.append('%autoreload 2')
#c.InteractiveShellApp.exec_lines.append('print "Warning: disable autoreload in startup.py to improve performance." ')

# print("\n  ISSM development path correctly loaded")
# print("Current path is {}\n\n".format(ISSM_DIR))

