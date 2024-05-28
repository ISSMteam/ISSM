#!/usr/bin/env python3

import argparse
import os
import re
from sys import float_info
from traceback import format_exc

import numpy as np

# Avoid the following error on Jenkins, 
#
#   "Unable to init server: Could not connect: Connection refused
#
#   (runme.py:28445): Gdk-CRITICAL **: 02:23:15.525: gdk_cursor_new_for_display: assertion 'GDK_IS_DISPLAY (display)' failed"
#
import matplotlib
matplotlib.use('Agg')

try:
    from arch import archread
except: # ISSM_DIR is not on path
    import devpath

from arch import archread
from arch import archwrite
from GetAvailableTestIds import *
from GetIds import *
from IdToName import IdToName
from parallelrange import parallelrange
from loadmodel import loadmodel
from solve import solve
from importlib import import_module


def runme(id=None, exclude=None, benchmark='nightly', procedure='check', output='none', rank=1, numprocs=1):
    """runme - test deck for ISSM nightly runs

    In a test deck directory (for example, test/NightlyRun) the following
    command will launch all existing tests,

        ./runme.py

    To run tests 101 and 102,

        ./runme.py -i 101 102

    Options:
        -i/--id         Followed by the list of ids or (parts of) test names
                        requested
        -e/--exclude    Ids or (parts of) test names to be excluded (same
                        format as id). Does nothing if 'id' is specified with
                        different values.
        -b/--benchmark  'all'           : (all of the tests)
                        'nightly'       : (nightly run/daily run)
                        'validation'    : (validation)
                        'adolc'         : validation of adolc tests
                        'eismint'       : validation of eismint tests
                        'ismip'         : validation of ismip-hom tests
                        'mesh'          : validation of mesh tests
                        'qmu'           : validation of qmu tests
                        'referential'   : validation of referential tests
                        'slc'           : validation of slc tests
                        'thermal'       : validation of thermal tests
                        'tranforcing'   : validation of transient forcing tests
        -p/--procedure  'check'         : run the test (default)
                        'update'        : update the archive
                        'runFromNC'     : run from an existing nc file


    Usage:
        ./runme.py [option [args]]

    Examples:
        ./runme.py
        ./runme.py -i 101
        ./runme.py -i 'SquareShelf'
        ./runme.py -e 2001
        ./runme.py -e 'Dakota' --benchmark 'all'
        ./runme.py -i [[101, 102], ['Dakota', 'Slc']]

    NOTE:
    - Will only run test scripts whose names explicitly follow the convention,

        test<id>.py

    where <id> is any integer.

    TODO:
    - At '#disp test result', make sure precision of output matches that of
    MATLAB.
    - Check for failures that do not raise exceptions (for example, 'Standard
    exception'; see also jenkins/jenkins.sh). These should be counted as
    failures.
    - Add support for 'stoponerror'
    """

    # Get ISSM_DIR variable
    ISSM_DIR = os.environ['ISSM_DIR']

    # Process options
    # Get benchmark {{{
    if benchmark not in ['all', 'nightly', 'ismip', 'eismint', 'thermal', 'mesh', 'validation', 'tranforcing', 'adolc', 'slc', 'qmu']:
        print('runme warning: benchmark \'{}\' not supported, defaulting to test \'nightly\'.'.format(benchmark))
        benchmark = 'nightly'
    # }}}
    # Get procedure {{{
    if procedure not in ['check', 'update', 'runFromNC']:
        print('runme warning: procedure \'{}\' not supported, defaulting to test \'check\'.'.format(procedure))
        procedure = 'check'
    # }}}
    # Get output {{{
    if output not in ['nightly', 'none']:
        print('runme warning: output \'{}\' not supported, defaulting to test \'none\'.'.format(output))
        output = 'none'
    # }}}
    # Get rank and numprocs for multi-threaded runs {{{
    if (numprocs < rank):
        numprocs = 1
    # }}}
    # Get available test IDs {{{
    all_ids = GetAvailableTestIds()
    i1, i2 = parallelrange(rank, numprocs, len(all_ids)) # Get tests for this CPU only
    all_ids = all_ids[i1:i2 + 1]

    # Check value passed to id argument. If is a single integer or string, convert to a list
    requested_ids = GetIds(id)

    if len(requested_ids) > 0:
        ids_to_run = set(requested_ids).intersection(set(all_ids))
        benchmark = None
    else:
        # If no tests are specifically requested, do them all
        ids_to_run = set(all_ids)
    # }}}
    # Get excluded tests {{{
    exclude_ids = GetIds(exclude)
    ids_to_run = ids_to_run.difference(exclude_ids)
    # }}}
    if procedure == 'runFromNC':
        # bamg test
        ids_to_run = ids_to_run.difference([119, 514])
        # smbGEMB format is weird for the test
        ids_to_run = ids_to_run.difference([243, 244, 252, 253])
        # AMR runs where the index is missing from fieldnames
        ids_to_run = ids_to_run.difference([462, 463, 464, 465])
        # test247 solves for thermal and transient which makes it complex to check
        ids_to_run = ids_to_run.difference([247])
        # test 902 is running two models with different stepping
        ids_to_run = ids_to_run.difference([902])
        # size issue in 517 needs investigation
        ids_to_run = ids_to_run.difference([517])

    # Process IDs according to benchmarks {{{
    if benchmark == 'nightly':
        ids_to_run = ids_to_run.intersection(set(range(1, 1000)))
    elif benchmark == 'validation':
        ids_to_run = ids_to_run.intersection(set(range(1001, 2000)))
    elif benchmark == 'ismip':
        ids_to_run = ids_to_run.intersection(set(range(1101, 1200)))
    elif benchmark == 'eismint':
        ids_to_run = ids_to_run.intersection(set(range(1201, 1300)))
    elif benchmark == 'thermal':
        ids_to_run = ids_to_run.intersection(set(range(1301, 1400)))
    elif benchmark == 'mesh':
        ids_to_run = ids_to_run.intersection(set(range(1401, 1500)))
    elif benchmark == 'tranforcing':
        ids_to_run = ids_to_run.intersection(set(range(1501, 1503)))
    elif benchmark == 'referential':
        ids_to_run = ids_to_run.intersection(set(range(1601, 1603)))
    elif benchmark == 'slc':
        ids_to_run = ids_to_run.intersection(set(range(2001, 2500)))
    elif benchmark == 'adolc':
        ids_to_run = ids_to_run.intersection(set(range(3001, 3200)))
    elif benchmark == 'qmu':
        ids_to_run = ids_to_run.intersection(set((218, 234, 235, 417, 418, 420)).union(set(range(412, 414))))
    ids_to_run = list(ids_to_run)
    ids_to_run.sort()

    # }}}

    # Loop over tests and launch sequence
    root = os.getcwd()
    error_count = 0
    resulted_in_error = []
    for id in ids_to_run:
        print(("----------------starting:{}-----------------------".format(id)))
        try:
            # Execute test
            os.chdir(root)
            id_string = IdToName(id)
            print(("----------------running-----------------------"))
            if procedure == 'runFromNC':
                Tmod = import_module('test{}'.format(id))
            else:
                exec(compile(open('test{}.py'.format(id)).read(), 'test{}.py'.format(id), 'exec'), globals())

            # Update archive?
            archive_name = 'Archive' + str(id)
            if procedure == 'update':
                archive_file = os.path.join('..', 'Archives', archive_name + '.arch')
                if os.path.isfile(archive_file):
                    os.remove(archive_file)
                for k, fieldname in enumerate(field_names):
                    field = np.array(field_values[k], dtype=float)
                    if len(field.shape) == 1:
                        if np.size(field):
                            field = field.reshape(np.size(field), 1)
                        else:
                            field = field.reshape(0, 0)
                    elif len(field.shape) == 0:
                        field = field.reshape(1, 1)
                        # Matlab uses base 1, so use base 1 in labels
                    archwrite(archive_file, archive_name + '_field' + str(k + 1), field)
                print(("File {} saved. \n".format(os.path.join('..', 'Archives', archive_name + '.arch'))))
            elif procedure == 'runFromNC':
                print(("----------------loadingNC-----------------------"))
                mdl = loadmodel('test{}ma.nc'.format(id))
                for key in mdl.results.__dict__.keys():
                    if 'Solution' in key:
                        solvetype = re.split('Solution', key)[0]

                # Save the results, scrap them and solve
                loaded_res = mdl.results
                mdl.results = []
                mdl = solve(mdl, solvetype)

                # Loop on the field_names from the nightly test
                for k, fieldname in enumerate(Tmod.field_names):
                    try:
                        # First, look for indexing
                        if re.search(r'\d+$', fieldname):
                            index = int(re.search(r'\d+$', fieldname).group()) - 1
                            fieldname = fieldname[:re.search(r'\d+$', fieldname).start()]
                        elif 'FirstStep' in fieldname:
                            index = 0
                            fieldname = fieldname[:re.search('FirstStep', fieldname).start()]
                        elif 'SecondStep' in fieldname:
                            index = 1
                            fieldname = fieldname[:re.search('SecondStep', fieldname).start()]
                        elif 'ThirdStep' in fieldname:
                            index = 2
                            fieldname = fieldname[:re.search('ThirdStep', fieldname).start()]
                        else:
                            index = 0

                        # Then, check if the key exists in the loaded results
                        try:
                            reskeys = mdl.results.__dict__[solvetype + 'Solution'][index].__dict__.keys()
                        except TypeError:
                            # Most likely a steady state so no subscripting
                            reskeys = mdl.results.__dict__[solvetype + 'Solution'].__dict__.keys()
                        if fieldname not in reskeys:
                            sufixes = ["P1bubble", "P1bubbleCondensed", "LliboutryDuval", "CuffeyTemperate", "SSA", "HO", "FS", "P1xP", "P2xP",
                                       'MINI', 'MINIcondensed', 'TaylorHood', 'XTaylorHood', 'LATaylorHood', 'CrouzeixRaviart', 'LACrouzeixRaviart']
                            namedifs = {'Misfits': 'J',
                                        'D': 'DamageDbar',
                                        'F': 'DamageF',
                                        'MaterialsRheologyB': 'MaterialsRheologyBbar',
                                        'SedimentWaterHead': 'SedimentHead',
                                        'EplWaterHead': 'EplHead',
                                        'SedimentWaterHeadSubstep': 'SedimentHeadSubstep',
                                        'EplWaterHeadSubstep': 'EplHeadSubstep',
                                        'Volume': 'IceVolume',
                                        'Bed': 'Base',
                                        'SMB': 'SmbMassBalance'}

                            if fieldname in namedifs.keys():
                                # Some fields are not consistent
                                fieldname = namedifs[fieldname]
                            elif any([suf in fieldname for suf in sufixes]):
                                # Some tests have loops that mess up naming
                                try:
                                    sufix = sufixes[np.squeeze(np.where([suf in fieldname for suf in sufixes]))]
                                except TypeError:
                                    # Probably several matches; we take the last one which should be the one we're want to run (needs to be controlled in the list above)
                                    sufix = sufixes[np.squeeze(np.where([suf in fieldname for suf in sufixes]))[-1]]
                                fieldname = fieldname[:re.search(sufix, fieldname).start()]
                            elif fieldname.endswith("P") and index == 1:
                                # Looking for P2 but 2 refers to an index, so shift by -1
                                fieldname = fieldname[:-1]
                            else:
                                # Handle case where index selected above is part of the name
                                fieldname = fieldname + str(index + 1)
                        try:
                            field = mdl.results.__dict__[solvetype + 'Solution'][index].__dict__[fieldname]
                            loaded_field = loaded_res.__dict__[solvetype + 'Solution'][index].__dict__[fieldname]
                        except TypeError:
                            # Most likely a steady state so no subscripting
                            try:
                                field = mdl.results.__dict__[solvetype + 'Solution'].__dict__[fieldname]
                                loaded_field = loaded_res.__dict__[solvetype + 'Solution'].__dict__[fieldname]
                            except KeyError:
                                print("WARNING: {}{} does not exist and checking will be skipped".format(fieldname, index + 1))
                                continue
                        except KeyError:
                            print("WARNING: {}{} does not exist and checking will be skipped".format(fieldname, index + 1))
                            continue

                        ref = Tmod.field_values[k]
                        # Get tolerance
                        tolerance = Tmod.field_tolerances[k]
                        # Compute differences for the results computed from the nc file
                        error_diff = np.amax(np.abs(ref - field), axis=0) / (np.amax(np.abs(ref), axis=0) + float_info.epsilon)
                        if not np.isscalar(error_diff):
                            error_diff = error_diff[0]

                        # Compute the differences for the results of the nc file
                        load_diff = np.amax(np.abs(np.squeeze(ref) - loaded_field), axis=0) / (np.amax(np.abs(np.squeeze(ref)), axis=0) + float_info.epsilon)
                        if not np.isscalar(load_diff):
                            load_diff = load_diff[0]

                        # Display test result
                        if (np.any(error_diff > tolerance) or np.isnan(error_diff)) and (np.any(load_diff > tolerance) or np.isnan(load_diff)):
                            if abs(error_diff - load_diff) < tolerance:
                                print(('WARNING difference: {:7.2g} > {:7.2g} test id: {} field: {}{} differs from computation but equal to saved results'.format(error_diff, tolerance, id, fieldname, index + 1)))
                            else:
                                print(('ERROR difference: {:7.2g} > {:7.2g} test id: {} field: {}{} is false in both loaded and computed results'.format(error_diff, tolerance, id, fieldname, index + 1)))
                                error_count += 1
                                resulted_in_error.append(id)
                        elif (np.any(error_diff > tolerance) or np.isnan(error_diff)):
                            print(('ERROR   difference: {:7.2g} > {:7.2g} test id: {} test name: {} field: {}{}'.format(error_diff, tolerance, id, id_string, fieldname, index + 1)))
                            error_count += 1
                            resulted_in_error.append(id)
                        elif (np.any(load_diff > tolerance) or np.isnan(load_diff)):
                            print(('SAVEERROR difference: {:7.2g} > {:7.2g} test id: {} test name: {} saved result : {}{}'.format(load_diff, tolerance, id, id_string, fieldname, index + 1)))
                            error_count += 1
                            resulted_in_error.append(id)
                        else:
                            print(('SUCCESS difference: {:7.2g} < {:7.2g} test id: {} test name: {} field: {}{}'.format(error_diff, tolerance, id, id_string, fieldname, index + 1)))
                        # Display only if there are errors in the results

                    except Exception as message:
                        # Something went wrong; print failure message
                        print((format_exc()))
                        if output == 'nightly':
                            fid = open(os.path.join(ISSM_DIR, 'nightlylog', 'pythonerror.log'), 'a')
                            fid.write('%s' % message)
                            fid.write('\n------------------------------------------------------------------\n')
                            fid.close()
                            print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, fieldname)))
                        else:
                            print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, fieldname)))
                            raise RuntimeError(message)
            # Check test
            else:
                # Load archive
                if os.path.exists(os.path.join('..', 'Archives', archive_name + '.arch')):
                    archive_file = os.path.join('..', 'Archives', archive_name + '.arch')
                else:
                    raise IOError("Archive file '../Archives/{}.arch' does not exist.".format(archive_name))

                for k, fieldname in enumerate(field_names):
                    try:
                        # Get field and tolerance
                        field = np.array(field_values[k])
                        if len(field.shape) == 1:
                            if np.size(field):
                                field = field.reshape(np.size(field), 1)
                            else:
                                field = field.reshape(0, 0)
                        tolerance = field_tolerances[k]

                        # Compare to archive
                        # MATLAB uses base 1, so use base 1 in labels
                        archive = np.array(archread(archive_file, archive_name + '_field' + str(k + 1)))
                        # NOTE: str(np.array(None)) becomes 'None' but np.array(None) is never equal to None: it basically becomes a type of string in an array
                        if str(archive) == 'None':
                            raise NameError("Field name '" + archive_name + '_field' + str(k + 1) + "' does not exist in archive file.")
                        if np.shape(field) != np.shape(archive) and not np.shape(field) in [(1, 1), (0, 0), (1, 0), (0, 1)]:
                            field = field.T
                            if np.shape(field) != np.shape(archive):
                                raise RuntimeError("Field '{}' from test {} is malformed; shape is {}, should be {} or {}".format(fieldname, archive_name[7:], np.shape(field.T), np.shape(archive), np.shape(archive.T)))

                        error_diff = np.amax(np.abs(archive - field), axis=0) / (np.amax(np.abs(archive), axis=0) + float_info.epsilon)
                        if not np.isscalar(error_diff):
                            error_diff = error_diff[0]

                        # Display test result
                        if (np.any(error_diff > tolerance) or np.isnan(error_diff)):
                            print(('ERROR   difference: {:7.2g} > {:7.2g} test id: {} test name: {} field: {}'.format(error_diff, tolerance, id, id_string, fieldname)))
                            error_count += 1
                            resulted_in_error.append(id)
                        else:
                            print(('SUCCESS difference: {:7.2g} < {:7.2g} test id: {} test name: {} field: {}'.format(error_diff, tolerance, id, id_string, fieldname)))

                    except Exception as message:
                        # Something went wrong; print failure message
                        print((format_exc()))
                        if output == 'nightly':
                            fid = open(os.path.join(ISSM_DIR, 'nightlylog', 'pythonerror.log'), 'a')
                            fid.write('%s' % message)
                            fid.write('\n------------------------------------------------------------------\n')
                            fid.close()
                            print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, fieldname)))
                        else:
                            print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, fieldname)))
                            raise RuntimeError(message)

        except Exception as message:
            # Something went wrong; print failure message
            print((format_exc()))
            if output == 'nightly':
                fid = open(os.path.join(ISSM_DIR, 'nightlylog', 'pythonerror.log'), 'a')
                fid.write('%s' % message)
                fid.write('\n------------------------------------------------------------------\n')
                fid.close()
                print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, 'N/A')))
            else:
                print(('FAILURE difference: N/A test id: {} test name: {} field: {}'.format(id, id_string, 'N/A')))
                raise RuntimeError(message)

        print(("----------------finished:{}-----------------------".format(id)))

    if error_count > 0:
        print("{} errors were detected in test {}".format(error_count, np.unique(resulted_in_error)))
    return

if __name__ == '__main__':
    if 'PYTHONSTARTUP' in os.environ:
        PYTHONSTARTUP = os.environ['PYTHONSTARTUP']
        if os.path.exists(PYTHONSTARTUP):
            try:
                exec(compile(open(PYTHONSTARTUP).read(), PYTHONSTARTUP, 'exec'))
            except Exception as e:
                print("PYTHONSTARTUP error: ", e)
        else:
            print(("PYTHONSTARTUP file '{}' does not exist.".format(PYTHONSTARTUP)))

        parser = argparse.ArgumentParser(description='runme - test deck for ISSM nightly runs')
        parser.add_argument('-i', '--id', nargs='*', help='followed by the list of test ID\'s or names to run', default=[])
        parser.add_argument('-e', '--exclude', nargs='*', type=str, help='followed by the list of test ID\'s or names to exclude', default=[])
        parser.add_argument('-b', '--benchmark', help='nightly/ismip/eismint/thermal/mesh/...', default='nightly')
        parser.add_argument('-p', '--procedure', help='check/update', default='check')
        parser.add_argument('-o', '--output', help='nightly/daily/none', default='none')
        parser.add_argument('-r', '--rank', type=int, help='rank', default=1)
        parser.add_argument('-n', '--numprocs', type=int, help='numprocs', default=1)
        args = parser.parse_args()

        md = runme(args.id, args.exclude, args.benchmark, args.procedure, args.output, args.rank, args.numprocs)

        exit(md)
    else:
        print("PYTHONSTARTUP not defined in environment")
        raise RuntimeError()
