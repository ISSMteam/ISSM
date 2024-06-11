#!/usr/bin/env python3
import os
import sys

# Recover ISSM_DIR and USERNAME
ISSM_DIR = os.getenv('ISSM_DIR')
USERNAME = os.getenv('USER')
if ISSM_DIR is None:
    raise NameError('"ISSM_DIR" environment variable is empty. You should define ISSM_DIR in your .zshrc or .bashrc')

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
if ISSM_DIR + '/src/wrappers/python/.libs' not in sys.path:
    sys.path.append(ISSM_DIR + '/src/wrappers/python/.libs')

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

