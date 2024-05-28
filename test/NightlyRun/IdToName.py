#!/usr/bin/env python3


def IdToName(test_id):
    """IdToName - return name of test

    Usage:
        name = IdToName(test_id)
    """
    infile = open('test' + str(test_id) + '.py', 'r')
    file_text = infile.readline()

    string = '#Test Name:'
    name = file_text[len(string) + 1:-1]
    return name
