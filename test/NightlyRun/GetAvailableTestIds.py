#!/usr/bin/env python3

import os
import re


def GetAvailableTestIds():
    """GetAvailableTestIds - Returns a list of all the test files in this 
    directory that match a certain naming scheme

    Usage:
        available_test_ids = GetAvailableTestIds()

    TODO:
    - Add optional parameter that allows user to designate another directory to 
    poll for test files
    """

    test_file_list = [f for f in os.listdir('.') if re.match('test[0-9]+.py$', f)] # File name must follow the format "test<INTEGER>.py"
    return [int(re.search(r'\d+',f.split('.')[0]).group()) for f in test_file_list] # Retrieve test IDs
