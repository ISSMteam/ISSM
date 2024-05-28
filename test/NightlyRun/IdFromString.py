#!/usr/bin/env python3

from GetAvailableTestIds import *
from IdToName import IdToName
import os
import re

# use verbose = False to print output when this is called by command line


def IdFromString(string, verbose=False):
    """IdFromString - output ids from a given string

    Usage:
        ids = IdFromString(string)
    
    Examples:
        ids = IdFromString('Dakota')
        ids = IdFromString('Slc')
        ids = IdFromString(' * ')
    """

    # Check input
    if not isinstance(string, str):
        raise TypeError('IdFromString error message: input argument is not a string.')
    string = string.replace("'", '')
    string = string.replace('"', '')

    # Get the test ids and names and scan for matches
    ids = []
    idnames = []
    available_test_ids = GetAvailableTestIds()
    for i in available_test_ids:
        name = IdToName(i)
        if (string == ' * ') or (name is not None and string in name):
            ids.append(i)
            idnames.append(name)

    # Return if no test found
    if not ids:
        print("No test matches '%s'." % string)
        return ids

    # Display names
    if verbose:
        idnames = [i for _, i in sorted(zip(ids, idnames), key=lambda pair: pair[0])]

    ids.sort()

    if verbose:
        print("{} tests match '{}':".format(len(ids), string))
        for i in range(len(ids)):
            print("   {} : {}".format(ids[i], idnames[i]))
    return ids
