import os

from dakota_in_data import *
from expandresponses import *
from expandvariables import *
from helpers import *
from MatlabFuncs import *
from process_qmu_response_data import *
from qmupart2npart import qmupart2npart


def preqmu(md, options):
    """PREQMU - apply Quantification of Margins and Uncertainties techniques
    to a solution sequence (like stressbalance.py, progonstic.py, etc ...),
    using the Dakota software from Sandia.

    Options come from the solve.py routine. They can include Dakota options:

        qmufile : input file for Dakota
        ivap    : <currently unimplemented>
        iresp   : <currently unimplemented>
        imethod : <currently unimplemented>
        iparams : <currently unimplemented>
        ivar    : selection number for variables input (if several are specified in variables)
        iresp   : same thing for response functions
        imethod : same thing for methods
        iparams : same thing for params
    """

    print('preprocessing dakota inputs')
    qmufile = options.getfieldvalue('qmufile', 'qmu')
    # qmufile cannot be changed unless ????script.sh is also changed
    options.addfielddefault('ivar', 0)
    options.addfielddefault('iresp', 0)
    options.addfielddefault('imethod', 0)
    options.addfielddefault('iparams', 0)

    # When running in library mode, the in file needs to be called md.miscellaneous.name.qmu.in
    qmufile = md.miscellaneous.name

    # Retrieve variables and resposnes for this particular analysis.
    #print type(md.qmu.variables)
    #print md.qmu.variables.__dict__
    # Print ivar
    variables = md.qmu.variables  #[ivar]
    responses = md.qmu.responses  #[iresp]

    # Expand variables and responses
    #print variables.__dict__
    #print responses.__dict__
    variables = expandvariables(md, variables)
    responses = expandresponses(md, responses)

    # Go through variables and responses, and check they don't have more than
    # md.qmu.numberofpartitions values. Also determine numvariables and 
    # numresponses
    #{{{
    numvariables = 0
    variable_fieldnames = fieldnames(variables)
    for i in range(len(variable_fieldnames)):
        field_name = variable_fieldnames[i]
        fieldvariables = vars(variables)[field_name]
        for j in range(np.size(fieldvariables)):
            if strncmpi(fieldvariables[j].descriptor, 'scaled_', 7):
                npart = qmupart2npart(fieldvariables[j].partition)
                nt = fieldvariables[j].nsteps
                if nt == 1:
                    if str2int(fieldvariables[j].descriptor, 'last') > npart:
                        raise RuntimeError('preqmu error message: one of the expanded variables has more values than the number of partitions')
        numvariables = numvariables + np.size(vars(variables)[field_name])

    numresponses = 0
    response_fieldnames = fieldnames(responses)
    for i in range(len(response_fieldnames)):
        field_name = response_fieldnames[i]
        fieldresponses = vars(responses)[field_name]
        for j in range(np.size(fieldresponses)):
            if strncmpi(fieldresponses[j].descriptor, 'scaled_', 7):
                npart = partition_npart(fieldresponses[j].partition)
                if str2int(fieldresponses[j].descriptor, 'last') > npart:
                    raise RuntimeError('preqmu error message: one of the expanded responses has more values than the number of partitions')
        numresponses = numresponses + np.size(vars(responses)[field_name])
    # }}}

    # Create in file for Dakota
    #dakota_in_data(md.qmu.method[imethod], variables, responses, md.qmu.params[iparams], qmufile)
    dakota_in_data(md.qmu.method, variables, responses, md.qmu.params, qmufile)

    # Build a list of variables and responses descriptors. the list is not expanded.
    #{{{
    variabledescriptors = []
    # variable_fieldnames = fieldnames(md.qmu.variables[ivar])
    variable_fieldnames = fieldnames(md.qmu.variables)
    for i in range(len(variable_fieldnames)):
        field_name = variable_fieldnames[i]
    #fieldvariables = vars(md.qmu.variables[ivar])[field_name]
        fieldvariables = vars(md.qmu.variables)[field_name]
        if type(fieldvariables) in [list, np.ndarray]:
            for j in range(np.size(fieldvariables)):
                variabledescriptors.append(fieldvariables[j].descriptor)
        else:
            variabledescriptors.append(fieldvariables.descriptor)

    responsedescriptors = []
    # response_fieldnames = fieldnames(md.qmu.responses[iresp])
    response_fieldnames = fieldnames(md.qmu.responses)
    for i in range(len(response_fieldnames)):
        field_name = response_fieldnames[i]
    #fieldresponses = vars(md.qmu.responses[iresp])[field_name]
        fieldresponses = vars(md.qmu.responses)[field_name]
        if type(fieldresponses) in [list, np.ndarray]:
            for j in range(np.size(fieldresponses)):
                responsedescriptors.append(fieldresponses[j].descriptor)
        else:
            responsedescriptors.append(fieldresponses.descriptor)
    # }}}

    # Build a list of variable partitions
    variablepartitions = []
    variablepartitions_npart = []
    variablepartitions_nt = []
    variable_fieldnames = fieldnames(md.qmu.variables)
    for i in range(len(variable_fieldnames)):
        field_name = variable_fieldnames[i]
        fieldvariable = vars(md.qmu.variables)[field_name]
        if type(fieldvariable) in [list, np.ndarray]:
            for j in range(np.size(fieldvariable)):
                if fieldvariable[j].isscaled() or fieldvariable[j].isdistributed():
                    variablepartitions.append(fieldvariable[j].partition)
                    variablepartitions_npart.append(qmupart2npart(fieldvariable[j].partition))
                    if hasattr(fieldvariable[j], 'nsteps'):
                        variablepartitions_nt.append(fieldvariable[j].nsteps)
                    else:
                        variablepartitions_nt.append(1)
                else:
                    variablepartitions.append([])
                    variablepartitions_npart.append(0)
                    variablepartitions_nt.append(1)
        else:
            if fieldvariable.isscaled():
                variablepartitions.append(fieldvariable.partition)
                variablepartitions_npart.append(qmupart2npart(fieldvariable.partition))
                if hasattr(fieldvariable, 'nsteps'):
                    variablepartitions_nt.append(fieldvariable.nsteps)
                else:
                    variablepartitions_nt.append(1)
            else:
                variablepartitions.append([])
                variablepartitions_npart.append(0)
                variablepartitions_nt.append(1)

    # Build a list of response partitions
    responsepartitions = []
    responsepartitions_npart = np.array([])
    response_fieldnames = fieldnames(md.qmu.responses)
    for i in range(len(response_fieldnames)):
        field_name = response_fieldnames[i]
        fieldresponse = vars(md.qmu.responses)[field_name]
        if type(fieldresponses) in [list, np.ndarray]:
            for j in range(np.size(fieldresponses)):
                if fieldresponse[j].isscaled():
                    responsepartitions.append(fieldresponse[j].partition)
                    responsepartitions_npart = np.append(responsepartitions_npart, qmupart2npart(fieldresponse[j].partition))
                else:
                    responsepartitions.append([])
                    responsepartitions_npart = np.append(responsepartitions_npart, 0)
        else:
            if fieldresponse.isscaled():
                responsepartitions.append(fieldresponse.partition)
                responsepartitions_npart = np.append(responsepartitions_npart, qmupart2npart(fieldresponse.partition))
            else:
                responsepartitions.append([])
                responsepartitions_npart = np.append(responsepartitions_npart, 0)

    if responsepartitions_npart.shape[0] != 1:
        responsepartitions_npart = responsepartitions_npart.reshape(1, -1)

    # Register the fields that will be needed by the Qmu model.
    md.qmu.numberofresponses = numresponses
    md.qmu.variabledescriptors = variabledescriptors
    md.qmu.variablepartitions = variablepartitions
    md.qmu.variablepartitions_npart = variablepartitions_npart
    md.qmu.variablepartitions_nt = variablepartitions_nt
    md.qmu.responsedescriptors = responsedescriptors
    md.qmu.responsepartitions = responsepartitions
    md.qmu.responsepartitions_npart = responsepartitions_npart

    # Now, we have to provide all the info necessary for the solutions to 
    # compute the responses. For example, if mass_flux is a response, we need a 
    # profile of points. For a misfit, we need the observed velocity, etc.
    md = process_qmu_response_data(md)

    return md
