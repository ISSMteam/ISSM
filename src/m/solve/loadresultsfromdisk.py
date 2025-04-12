import os
from helpers import fieldnames, struct
from parseresultsfromdisk import parseresultsfromdisk
from postqmu import postqmu
from results import results

def loadresultsfromdisk(md, filename):
    """LOADRESULTSFROMDISK - load results of solution sequence from disk file
    "filename"

    Usage:
        md = loadresultsfromdisk(md=False, filename=False)
    """

    # Check number of inputs/outputs
    if not md or not filename:
        raise ValueError('loadresultsfromdisk: error message.')

    if not md.qmu.isdakota:

        # Check that file exists
        if not os.path.exists(filename):
            err_msg = '==========================================================================\n'
            err_msg += '   Binary file {} not found                                              \n'.format(filename)
            err_msg += '                                                                         \n'
            err_msg += '   This typically results from an error encountered during the simulation\n'
            err_msg += '   Please check for error messages above or in the outlog                \n'
            err_msg += '=========================================================================\n'
            print(err_msg)
            return

        # Initialize md.results if it is not a structure yet
        if not isinstance(md.results, results):
            md.results = results()

        # Load results onto model
        structure = parseresultsfromdisk(md, filename, not md.settings.io_gather)
        if not structure:
            raise RuntimeError('No result found in binary file \'{}\'. Check for solution crash.'.format(filename))
        if not hasattr(structure[0], 'SolutionType'):
            if hasattr(structure[-1], 'SolutionType'):
                structure[0].SolutionType = structure[-1].SolutionType
            else:
                print('Cannot find a solution type in the results! Ascribing one: \'NoneSolution\'.')
                structure[0].SolutionType = 'NoneSolution'
        setattr(md.results, structure[0].SolutionType, structure)

        # Recover solution_type from results
        md.private.solution = structure[0].SolutionType

        # Read log files onto fields
        if os.path.exists(md.miscellaneous.name + '.errlog'):
            with open(md.miscellaneous.name + '.errlog', 'r') as f:
                setattr(getattr(md.results, structure[0].SolutionType)[0], 'errlog', [line[:-1] for line in f])
        else:
            setattr(getattr(md.results, structure[0].SolutionType)[0], 'errlog', '')

        if os.path.exists(md.miscellaneous.name + '.outlog'):
            with open(md.miscellaneous.name + '.outlog', 'r') as f:
                setattr(getattr(md.results, structure[0].SolutionType)[0], 'outlog', [line[:-1] for line in f])
        else:
            setattr(getattr(md.results, structure[0].SolutionType)[0], 'outlog', '')

        if getattr(md.results, structure[0].SolutionType)[0].errlog:
            print('loadresultsfromdisk info message: error during solution. Check your errlog and outlog model fields.')

        # If only one solution, extract it from list for user friendliness
        if len(structure) == 1 and structure[0].SolutionType != 'TransientSolution':
            setattr(md.results, structure[0].SolutionType, structure[0])

    # Postprocess QMU results if necessary
    else:
        md = postqmu(md)

    return md
