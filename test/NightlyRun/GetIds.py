#!/usr/bin/env python3

from IdToName import *
from IdFromString import *
import numpy as np


def GetIds(ids):
    """GetIds - Output IDs from a given list of IDs, strings representing 
    MATLAB-like integer ranges, and test names

    The test names can be any string or substring present in the test's name 
    (first line of corresponding file). Test names are case sensitive.

    Usage:
        ids = GetIds(101)
        ids = GetIds(101:111)
        ids = GetIds('Dakota')
        ids = GetIds([101, 102[, ...])
        ids = GetIds(['Dakota', 'Slc'[, ...]])
        ids = GetIds([101, 'Dakota', 'Slc'[, ...]])
    """

    output_ids = []

    for id in ids:
        if id.isdigit():
            output_ids.append(int(id))
        else:
            if ':' in id: # MATLAB-like range of test ID's
                id_range = id.split(':')
                for i in range(int(id_range[0]), int(id_range[1]) + 1):
                    output_ids.append(i)
            else: # Test name
                for id_from_string in IdFromString(id):
                    output_ids.append(id_from_string)

    return output_ids
