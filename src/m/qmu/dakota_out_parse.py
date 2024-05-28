import numpy as np
from os.path import isfile, getsize
import re
from MatlabFuncs import *
from prctile_issm import *
from normfit_issm import *
#move this later
from helpers import *

# NOTE: May be rewritten later to take advantage of Python's file I/O 
# mechanics. As it is written now, it is often difficult to work with, but is 
# analagous to the MATLAB version of dakota_out_parse.


def dakota_out_parse(filei):  # {{{
    """DAKOTA_OUT_PARSE - read a Dakota .out or .dat output file and parse it.

    Usage:
        [method, dresp, scm, pcm, srcm, prcm] = dakota_out_parse(filei)

    where the required input is,
        filei         (character, name of .out file)

    the required output is,
        method        (character, Dakota method name)
        dresp         (structure array, responses)

    and the optional output is,
        scm           (double array, simple correlation matrix)
        pcm           (double array, partial correlation matrix)
        srcm          (double array, simple rank correlation matrix)
        prcm          (double array, partial rank correlation matrix)

    The filei will be prompted for if empty. The fields of dresp are particular 
    to the data contained within the file. The scm, pcm, srcm, and prcm are 
    output by Dakota only for the sampling methods.

    This function reads a Dakota .out output file and parses it into the Python 
    runtime. It operates in a content-driven fashion, where it skips the 
    intermediate data and then parses whatever output data it encounters in the 
    order in which it exists in the file, rather than searching for data based 
    on the particular method (this makes it independent of method). It also can 
    read and parse the .dat tabular_output file.

    This data would typically be used for plotting and other postprocessing 
    within MATLAB or Excel.

    TODO:
    - Figure out why output from Dakota is different under MATLAB and Python 
    (is it the input file that we write?)

    "Copyright 2009, by the California Institute of Technology. ALL RIGHTS 
    RESERVED. United States Government Sponsorship acknowledged. Any commercial 
    use must be negotiated with the Office of Technology Transfer at the 
    California Institute of Technology. (NTR 47078)

    This software may be subject to U.S. export control laws. By accepting this 
    software, the user agrees to comply with all applicable U.S. export laws 
    and regulations. User has the responsibility to obtain export licenses, or 
    other export authority as may be required before exporting such information 
    to foreign countries or providing access to foreign persons."
    """

    if filei is None:
        help(dakota_out_parse)
        return

    if not isfile(filei) or getsize(filei) == 0:
        filei = str(eval(input('Input file? ')))

    #fidi = fopen(sprintf('%s', filei), 'r')
    #try:
    with open(filei, 'r') as fidi:
        #  check the first line for the Dakota tabular output file
        method = []
        fline = fidi.readline()
        if getsize(filei) == 0 or fline == '':
            raise RuntimeError('File ' + filei + ' is empty')

        dresp = [] # of struct()
        scm = struct()
        pcm = struct()
        srcm = struct()
        prcm = struct()

        if '%eval_id' in fline:
            method = 'unknown'
            dresp = dak_tab_out(fidi, fline)
            return [method, dresp, scm, pcm, srcm, prcm]
        else:
            fidi.seek(0, 0)

        # loop through the file to find the Dakota method name
        fline = findline(fidi, 'method', True)
        if fline is None:
            #do nothing
            pass
        else:
            if fline[6] == ',':
                fline = fidi.readline()
                [ntokens, tokens] = fltokens(fline)
                method = tokens[0].strip()
                print('Dakota method = \'' + method + '\'')
            elif fline[6] in ['N', 'n']:
                fline = findline(fidi, 'methodName = ')
                [ntokens, tokens] = fltokens(fline)
                method = tokens[2].strip()
                print('Dakota methodName = \'' + method + '\'')

        # loop through the file to find the function evaluation summary
        counter = 0
        fline = ''
        nfeval = nfeval_read(fidi, fline)

        # process each results section based on content of the file
        while counter < 10:
            # because python makes file I/O difficult
            # if we see 10 + blank lines in a row then we have reached EOF
            # (tests show actual maximum number of blank lines is around 5)
            if fline == '' or fline.isspace():
                counter += 1
            else:
                counter = 0
            # ipos = ftell(fidi)
            fline = fidi.readline()
            if fline == '' or fline.isspace():
                pass
            elif '<<<<< Function evaluation summary' in fline:
                nfeval = nfeval_read(fidi, fline)
            elif 'Statistics based on ' in fline:
                nsamp = nsamp_read(fidi, fline)
            elif 'Moments for each response function' in fline:
                dresp = moments_read(fidi, dresp, fline)
            elif 'Moment-based statistics for each response function' in fline:
                dresp = mbstats_read(fidi, dresp, fline)
            elif '95% confidence intervals for each response function' in fline:
                dresp = cis_read(fidi, dresp, fline)
            elif 'Probabilities for each response function' in fline or 'Level mappings for each response function' in fline:
                dresp = cdfs_read(fidi, dresp, fline)
            elif 'Probability Density Function (PDF) histograms for each response function' in fline:
                dresp = pdfs_read(fidi, dresp, fline)
            elif 'Simple Correlation Matrix' in fline:
                scm = corrmat_read(fidi, 'Simple Correlation Matrix', fline)
            elif 'Partial Correlation Matrix' in fline:
                pcm = corrmat_read(fidi, 'Partial Correlation Matrix', fline)
            elif 'Simple Rank Correlation Matrix' in fline:
                srcm = corrmat_read(fidi, 'Simple Rank Correlatio:n Matrix', fline)
            elif 'Partial Rank Correlation Matrix' in fline:
                prcm = corrmat_read(fidi, 'Partial Rank Correlation Matrix', fline)
            elif 'MV Statistics for ' in fline:
                dresp = mvstats_read(fidi, dresp, fline)
            elif '<<<<< Best ' in fline:
                dresp = best_read(fidi, dresp, fline)
            elif 'The following lists volumetric uniformity measures' in fline:
                dresp = vum_read(fidi, dresp, fline)
            elif '<<<<< Iterator ' in fline and (len(fline) > 26) and (' completed.' in fline[15:]):
                method = itcomp_read(fidi, fline)
            elif '-----' in fline:
                pass
            else:
                'Unexpected line: ' + str(fline)

            # fidi.seek(ipos, 0)

    # loop through the file to verify the end

    # fline = findline(fidi, '<<<<< Single Method Strategy completed')
    # if not ischar(fline)
    #     return
    #
        print('End of file successfully reached')
    #close(fidi)
    #except Exception as err:
    #print "ERROR in dakota_out_parse: " + err
    #raise err
    #raise RuntimeError(filei + ' could not be opened')

    return [method, dresp, scm, pcm, srcm, prcm]
    # }}}


def dak_tab_out(fidi, fline):  # {{{
    """DAK_TAB_OUT - function to parse the Dakota tabular output file
    """

    print('Reading Dakota tabular output file')

    # Process column headings of matrix (skipping eval_id)
    [ntokens, tokens] = fltokens(fline)

    if strncmpi(fline, '%eval_id interface', 18): # Dakota versions >= 6
        offset = 2
    else:  # Dakota versions < 6
        offset = 1

    desc = ['' for i in range(ntokens - offset)]
    data = np.zeros((1, ntokens - offset))

    for i in range(ntokens - offset):
        desc[i] = str(tokens[i + offset])

    print('Number of columns (Dakota V + R) = {}'.format(ntokens - 2))

    # Process rows of matrix
    nrow = 0
    while True:
        fline = fidi.readline()

        if fline == '' or fline.isspace():
            break

        if nrow > 0:
            data = np.concatenate((data, [np.zeros(ntokens - offset)]))

        [ntokens, tokens] = fltokens(fline)

        # Add row values to matrix (skipping eval_id)
        for i in range(ntokens - offset):
            data[nrow, i] = tokens[i + offset]

        nrow = nrow + 1

    print('Number of rows (Dakota func evals) = ' + str(nrow))

    # Calculate statistics
    if (np.size(data, 0) > 1):
        #dmean = mean(data)
        #dstddev = std(data, 0)
        [dmean, dstddev, dmeanci, dstddevci] = normfit_issm(data, 0.05)
    else:
        dmean = np.zeros((1, np.size(data, 1)))
        dstddev = np.zeros((1, np.size(data, 1)))
        dmeanci = np.zeros((2, np.size(data, 1)))
        dstddevci = np.zeros((2, np.size(data, 1)))
        for i in range(np.size(data, 1)):
            [dmean[0, i], dstddev[0, i], dmeanci[:, i], dstddevci[:, i]] = normfit_issm(data[:, i], 0.05)

    dmin = data.min(axis=0)
    dquart1 = prctile_issm(data, 25, 0)
    dmedian = np.median(data, axis=0)
    dquart3 = prctile_issm(data, 75, 0)
    dmax = data.max(axis=0)
    dmin95 = prctile_issm(data, 5, 0)
    dmax95 = prctile_issm(data, 95, 0)

    # NOTE: The following line may cause the following warning (should not 
    # crash or invalidate results) when one of the inputs does not change with 
    # respect to the other(s), causing an internal divide-by-zero error,
    #
    #       /usr/local/lib/python2.7/dist-packages/numpy/lib/function_base.py:3163:
    #       RuntimeWarning: invalid value encountered in true_divide
    #       c /= stddev[:, None]
    #
    # (and/or the same but with "c /= stddev[None, :]")

    # Equivalent to Dakota scm, MATLAB corrcoef, and Excel correl
    dcorrel = np.corrcoef(data.T)

    # Divide the data into structures for consistency
    dresp = []
    for i in range(len(desc)):
        dresp.append(struct())
        dresp[i].descriptor = desc[i]
        dresp[i].sample = data[:, i]
        dresp[i].mean = dmean[i]
        dresp[i].stddev = dstddev[i]
        dresp[i].meanci = dmeanci[:, i]
        dresp[i].stddevci = dstddevci[:, i]
        dresp[i].min = dmin[i]
        dresp[i].quart1 = dquart1[i]
        dresp[i].median = dmedian[i]
        dresp[i].quart3 = dquart3[i]
        dresp[i].max = dmax[i]
        dresp[i].dmin95 = dmin95[i]
        dresp[i].dmax95 = dmax95[i]

    #  draw box plot

    # figure
    # subplot(2, 1, 1)
    # plot_boxplot(dresp)

    #  draw normal probability plot

    # subplot(2, 1, 2)
    # plot_normplot(dresp)

    return dresp
    # }}}


def nfeval_read(fidi, fline):  # {{{
    #  function to find and read the number of function evaluations

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, '<<<<< Function evaluation summary')
        nfeval = 0
        return

    [ntokens, tokens] = fltokens(fline)
    nfeval = tokens[4]
    print('  Dakota function evaluations = ' + str(int(nfeval)))

    return nfeval
    # }}}


def nsamp_read(fidi, fline):  # {{{
    #  function to find and read the number of samples

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'Statistics based on ')
        return

    [ntokens, tokens] = fltokens(fline)
    nsamp = tokens[3]
    print('  Dakota samples = ' + str(int(nsamp)))

    return nsamp
    # }}}


def moments_read(fidi, dresp, fline):  # {{{
    #  function to find and read the moments

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'Moments for each response function')
        return

    print('Reading moments for response functions:')

    while True:
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        [ntokens, tokens] = fltokens(fline)

    #  add new response function and moments

        dresp.append(struct())
        dresp[-1].descriptor = tokens[0]
        print('  ' + str(dresp[-1].descriptor))
        dresp[-1].mean = tokens[3]
        dresp[-1].stddev = tokens[6]
        dresp[-1].coefvar = tokens[12]

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def mbstats_read(fidi, dresp, fline):  # {{{
    #  function to find and read the moment - based statistics

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'Moment-based statistics for each response function')
        return

    print('Reading moment-based statistics for response functions:')

    #  skip column headings of moment - based statistics

    fline = fidi.readline()

    while True:
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        [ntokens, tokens] = fltokens(fline)

    #  add new response function and moment - based statistics

        dresp.append(struct())
        dresp[-1].descriptor = tokens[0]
        print('  ' + str(dresp[-1].descriptor))
        dresp[-1].mean = tokens[1]
        dresp[-1].stddev = tokens[2]
        dresp[-1].skewness = tokens[3]
        dresp[-1].kurtosis = tokens[4]

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def cis_read(fidi, dresp, fline):  # {{{
    #  function to find and read the confidence intervals

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, '95% confidence intervals for each response function')
        return

    print('Reading 95% confidence intervals for response functions:')

    while True:
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        [ntokens, tokens] = fltokens(fline)
        #  check for column headings in Dakota 5.2
        if (ntokens == 4):
            fline = fidi.readline()
            if fline == '' or fline.isspace():
                break

            [ntokens, tokens] = fltokens(fline)

        #  find response function associated with confidence intervals
        idresp = -1
        for i in range(len(dresp)):
            if strcmpi(tokens[0], dresp[i].descriptor):
                idresp = i
                break

        if idresp < 0:
            idresp = len(dresp)
            dresp.append(struct())
            dresp[idresp].descriptor = tokens[0]
            print('  ' + str(dresp[idresp].descriptor))

        #  add confidence intervals to response functions
        dresp[i].meanci = np.array([[np.nan], [np.nan]])
        dresp[i].stddevci = np.array([[np.nan], [np.nan]])

        if (ntokens == 14):
            dresp[i].meanci[0, 0] = tokens[4]
            dresp[i].meanci[1, 0] = tokens[5]
            dresp[i].stddevci[0, 0] = tokens[11]
            dresp[i].stddevci[1, 0] = tokens[12]
        else:
            dresp[i].meanci[0, 0] = tokens[1]
            dresp[i].meanci[1, 0] = tokens[2]
            dresp[i].stddevci[0, 0] = tokens[3]
            dresp[i].stddevci[1, 0] = tokens[4]

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def cdfs_read(fidi, dresp, fline):  # {{{
    #  function to find and read the cdf's

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'Probabilities for each response function')
        if fline is None:
            fline = findline(fidi, 'Level mappings for each response function')
            if fline is None:
                return

    print('Reading CDF''s for response functions:')

    while fline == '' or fline.isspace():
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

    #  process header line of cdf

        while (fline != '' and not fline.isspace()):
            [ntokens, tokens] = fltokens(fline)

    #  find response function associated with cdf
    # idresp is an index, so it can be 0, default to - 1
            idresp = -1
            for i in range(len(dresp)):
                if strcmpi(tokens[5], dresp[i].descriptor):
                    idresp = i
                    break
            if idresp < 0:
                idresp = len(dresp)
                dresp.append(struct())
                dresp[idresp].descriptor = tokens[5]
                print('  ' + str(dresp(idresp).descriptor))

            #  skip column headings of cdf
            fline = fidi.readline()
            fline = fidi.readline()

            #  read and add cdf table to response function
            fline = fidi.readline()
            icdf = 0
            while (fline != '' and not fline.isspace()) and not strncmpi(fline, 'Cumulative Distribution Function', 32):
                [ntokens, tokens] = fltokens(fline)
                icdf = icdf + 1
                dresp[idresp].cdf = np.zeros((icdf, 4))
                dresp[idresp].cdf[icdf - 1, 0:4] = np.nan
                #  in later versions of Dakota, uncalculated columns are now blank
                itoken = 0
                for i in range(len(fline) / 19):
                    if not isempty(fline[(i - 1) * 19:i * 19]):
                        itoken = itoken + 1
                        dresp[idresp].cdf[icdf - 1, i] = tokens[itoken]

                fline = fidi.readline()

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def pdfs_read(fidi, dresp, fline):  # {{{
    #  function to find and read the pdf's

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'Probability Density Function (PDF) histograms for each response function')
        return

    print('Reading PDF''s for response functions:')

    while (fline != '' and not fline.isspace()):
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        #  process header line of pdf
        while (fline != '' and not fline.isspace()):
            [ntokens, tokens] = fltokens(fline)

            #  find response function associated with pdf
            # idresp is an index, so it can be 0, default to - 1
            idresp = -1
            for i in range(len(dresp)):
                if strcmpi(tokens[2], dresp[i].descriptor):
                    idresp = i
                    break

            if idresp < 0:
                idresp = len(dresp)
                dresp.append(struct)
                dresp[idresp].descriptor = tokens[2]
                print('  ' + str(dresp[idresp].descriptor))

            #  skip column headings of pdf
            fline = fidi.readline()
            fline = fidi.readline()

            #  read and add pdf table to response function
            fline = fidi.readline()
            ipdf = 0
            while (fline != '' and not fline.isspace()) and not strncmpi(fline, 'PDF for', 7):
                [ntokens, tokens] = fltokens(fline)
                ipdf = ipdf + 1
                dresp[idresp].pdf = np.zeros((ipdf, 4))
                dresp[idresp].pdf[ipdf - 1, 0:3] = np.nan
                for i in range(3):
                    dresp[idresp].pdf[ipdf - 1, i] = tokens[i]

                fline = fidi.readline()

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def corrmat_read(fidi, cmstr, fline):  # {{{
    #  function to find and read a correlation matrix

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, cmstr)
        if fline == '' or fline.isspace():
            cmat = struct()
            return

    print('Reading ' + fline)

    cmat.title = fline

    while (fline != '' and not fline.isspace()):
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        #  process column headings of matrix
        [ntokens, tokens] = fltokens(fline)
        cmat.column = np.empty((1, ntokens))
        cmat.column.fill(0.0)
        cmat.row = np.empty((1, 1))
        cmat.row.fill(0.0)
        cmat.matrix = np.zeros((1, ntokens))

        for i in range(ntokens):
            cmat.column[1, i] = str(tokens[i])

        #  process rows of matrix, reading until blank line
        nrow = 0
        while True:
            fline = fidi.readline()
            if fline == '' or fline.isspace():
                break

            [ntokens, tokens] = fltokens(fline)

            #  add row heading to matrix
            nrow = nrow + 1
            cmat.row[nrow - 1, 0] = str(tokens[0])

            #  add row values to matrix
            for i in range(1, ntokens):
                cmat.matrix[nrow - 1, i - 1] = tokens[i]

    return cmat
    # }}}


def mvstats_read(fidi, dresp, fline):  # {{{
    #  function to find and read the MV statistics

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'MV Statistics for ')
        if fline is None:
            return

    print('Reading MV statistics for response functions:')

    ndresp = 0

    while (fline != '' and not fline.isspace()) and strncmpi(fline, 'MV Statistics for ', 18):

        #  add new response function and moments
        [ntokens, tokens] = fltokens(fline)
        dresp.append(struct())
        dresp[-1].descriptor = tokens[3]
        print('  ' + str(dresp[-1].descriptor))
        fline = fidi.readline()
        [ntokens, tokens] = fltokens(fline)
        dresp[-1].mean = tokens[4]
        fline = fidi.readline()
        [ntokens, tokens] = fltokens(fline)
        dresp[-1].stddev = tokens[6]

        #  read and add importance factors to response function
        idvar = 0
        fline = fidi.readline()
        if fline == '' or fline.isspace():
            break

        # shape: [[0], [0], [0]...]
        dresp[-1].var = []
        dresp[-1].impfac = []
        dresp[-1].sens = []

        while (fline != '' and not fline.isspace()) and strncmpi(fline, '  Importance Factor for variable ', 33):
            [ntokens, tokens] = fltokens(fline)
            idvar = idvar + 1
            dresp[-1].var.append(str(tokens[4]))
            dresp[-1].impfac.append(tokens[6])
            if (ntokens >= 10):
                dresp[-1].sens.append(tokens[9])
            else:
                dresp[-1].sens.append(np.nan)

            fline = fidi.readline()

        #  if importance factors missing, skip to cdf
        if not idvar:
            print('    Importance Factors not available')
            dresp[-1].var = []
            dresp[-1].impfac = []
            dresp[-1].sens = []
            while type(fline) == str and (fline != '' and not fline.isspace()) and not strncmpi(fline, 'Cumulative Distribution Function', 32) and not strncmpi(fline, 'MV Statistics for ', 18) and not strncmp(fline, ' - ', 1):
                fline = fidi.readline()

        #  process header line of cdf
        icdf = 0

        # If there is a warning it MAY involve a lot of spaces; skip over them
        if fline == '' or fline.isspace():
            fline = fidi.readline()
            # Usually: "Warning: negligible standard deviation renders CDF results suspect."
            if strncmpi(fline, 'Warn', 4):
                fline = fidi.readline()
                if fline == '' or fline.isspace():
                    fline = fidi.readline()

        while (fline != '' and not fline.isspace()) and strncmpi(fline, 'Cumulative Distribution Function', 32):
            [ntokens, tokens] = fltokens(fline)

            #  find response function associated with cdf
            # idresp is an index, so it can be 0, default to - 1
            idresp = -1
            for i in range(len(dresp)):
                if strcmpi(tokens[5], dresp[i].descriptor):
                    idresp = i
                    break

            if idresp < 0:
                idresp = len(dresp)
                dresp.append(struct())
                dresp[idresp].descriptor = tokens[5]
                print('  ' + str(dresp[idresp].descriptor))

            #  skip column headings of cdf
            fline = fidi.readline()
            fline = fidi.readline()

            #  read and add cdf table to response function
            fline = fidi.readline()
            while (fline != '' and not fline.isspace()) and not strncmpi(fline, 'MV Statistics for ', 18) and not strncmp(fline, ' - ', 1):
                [ntokens, tokens] = fltokens(fline)
                icdf = icdf + 1
                dresp[idresp].cdf = np.zeros((icdf, 4))
                dresp[idresp].cdf[icdf - 1, 0] = tokens[0]
                dresp[idresp].cdf[icdf - 1, 1] = tokens[1]
                if (ntokens == 4):
                    dresp[idresp].cdf[icdf - 1, 2] = tokens[2]
                    dresp[idresp].cdf[icdf - 1, 3] = tokens[3]
                else:
                    dresp[idresp].cdf[icdf - 1, 2] = np.nan
                    dresp[idresp].cdf[icdf - 1, 3] = np.nan

                fline = fidi.readline()

        #  if cdf missing, skip to end of response function
        if not icdf:
            print('    Cumulative Distribution Function not available')
            dresp[ndresp].cdf = []
            while (fline != '' and not fline.isspace()) and not strncmpi(fline, 'MV Statistics for ', 18) and not strncmp(fline, ' - ', 1):
                fline = fidi.readline()

    print('  Number of Dakota response functions = ' + str(len(dresp)))

    return dresp
    # }}}


def best_read(fidi, dresp, fline):  # {{{
    #  function to find and read the best evaluation

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, ' < < < < < Best ')
        if fline is None:
            return

    if isempty(dresp):
        dresp.append(struct())
        dresp[-1].best = struct()

    print('Reading values for best function evaluation:')

    while (fline != '' and not fline.isspace()) and strncmpi(fline, ' < < < < < Best ', 11):
        [ntokens, tokens] = fltokens(fline)

    #  read and add best parameter(s)

        if strncmpi(str(tokens[2]), 'parameter', 9):
            print('  ' + fline)

            fline = fidi.readline()
            dresp.best.param = []
            dresp.best.descriptor = ''

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, ' < < < < < Best ', 11):
                [ntokens, tokens] = fltokens(fline)
                dresp.best.param.append([0])
                dresp.best.param[-1] = tokens[0]
                dresp.best.descriptor = str(tokens[1])
                fline = fidi.readline()

        #  read and add best objective function(s)
        elif strncmpi(str(tokens[2]), 'objective', 9) and strncmpi(str(tokens[3]), 'function', 8):
            print('  ' + fline)

            fline = fidi.readline()
            dresp.best.of = []

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, ' < < < < < Best ', 11):
                [ntokens, tokens] = fltokens(fline)
                dresp.best.of.append(0)
                dresp.best.of[-1] = tokens[0]
                fline = fidi.readline()

        #  read and add best residual norms
        elif strncmpi(str(tokens[2]), 'residual', 8) and strncmpi(str(tokens[3]), 'norm', 4):
            print('  ' + fline)
            dresp.best.norm = tokens[5]
            dresp.best.hnormsq = tokens[10]

            fline = fidi.readline()

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, ' < < < < < Best ', 11):
                fline = fidi.readline()

        #  read and add best residual term(s)
        elif strncmpi(str(tokens[2]), 'residual', 8) and strncmpi(str(tokens[3]), 'term', 4):
            print('  ' + fline)

            fline = fidi.readline()
            dresp.best.res = []

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, '<<<<<Best ', 11):
                [ntokens, tokens] = fltokens(fline)
                dresp.best.res.append(0)
                dresp.best.res[-1] = tokens[0]
                fline = fidi.readline()

        #  read and add best constraint value(s)
        elif strncmpi(str(tokens[2]), 'constraint', 10) and strncmpi(str(tokens[3]), 'value', 5):
            print('  ' + fline)

            fline = fidi.readline()
            dresp.best.nc = []

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, '<<<<<Best ', 11):
                [ntokens, tokens] = fltokens(fline)
                dresp.best.nc.append(0)
                dresp.best.nc[-1] = tokens[0]
                fline = fidi.readline()

        #  read and add best data captured
        elif strncmpi(str(tokens[2]), 'data', 4) and strncmpi(str(tokens[3]), 'captured', 8):
            print('  ' + fline)
            dresp.best.eval = tokens[7]

            fline = fidi.readline()

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, ' < < < < < Best ', 11):
                fline = fidi.readline()

        #  read until next best or blank or end
        else:
            print('  ' + fline + '  (ignored)')

            fline = fidi.readline()

            while (fline != '' and not fline.isspace()) and not strncmpi(fline, ' < < < < < Best ', 11):
                fline = fidi.readline()

    return dresp
    # }}}


def vum_read(fidi, dresp, fline):  # {{{
    #  function to find and read the volumetric uniformity measures

    if fline is None or fline == '' or fline.isspace():
        fline = findline(fidi, 'The following lists volumetric uniformity measures')
        if fline is None:
            return

    if isempty(dresp):
        dresp.append(struct())
        dresp[-1].vum = []

    print('Reading measures for volumetric uniformity')
    fline = fidi.readline()
    fline = fidi.readline()

    while (fline != '' and not fline.isspace()):
        [ntokens, tokens] = fltokens(fline)
        check = tokens[0].lower()
        if check == 'chi':
            dresp.vum.chi = tokens[3]
        elif check == 'd':
            dresp.vum.d = tokens[3]
        elif check == 'h':
            dresp.vum.h = tokens[3]
        elif check == 'tau':
            dresp.vum.tau = tokens[3]

        fline = fidi.readline()

    return dresp
    # }}}


def itcomp_read(fidi, fline):  # {{{
    #  function to find and read the iterator completion

    if fline is None or fline == '' or fline.isspace():
        while True:
            fline = findline(fidi, '<<<<< Iterator ')
            if fline is None:
                return

            if (len(fline) > 26) and not (' completed.' in fline[15:]):
                break

    [ntokens, tokens] = fltokens(fline)
    method = tokens[2]
    print('Dakota iterator \'' + str(method) + '\' completed')

    return method
    # }}}


def findline(fidi, string, goto_line=False):  # {{{
    #  function to find a file line starting with a specified string
    #  by default, return to previous position, before search
    #  if final argument is True, return such that fidi.readline() will read the line
    #    immediately after the searched line (or the first line if the search failed)

    ipos = fidi.tell()
    npos = 0
    for fline in fidi:
        npos += len(fline)
        if (strncmpi(fline, string, len(string))):
            if goto_line:
                fidi.seek(npos, 0)
            else:
                fidi.seek(ipos, 0)
            return fline

    #  issue warning and reset file position
    print('Warning: findline:str_not_found: String ' + str(string) + ' not found in file')
    fidi.seek(ipos, 0)
    return None
    # }}}


def fltokens(fline):  # {{{
    #  function to parse a file line into tokens
    if fline is None:
        ntokens = -1
        tokens = []
        return [None, None]

    if fline == '' or fline.isspace():
        ntokens = 0
        tokens = []
        return [None, None]

    # split wherever ' ' (space) or ':' occur
    strings = re.split(':| ', fline)
    # remove blank strings
    strings = [a for a in strings if (a != '' and not a.isspace())]

    ntokens = 0
    tokens = ['' for i in range(len(strings))]

    # try to format substrings to float where possible and count tokens and ignore invalid values
    for i in range(len(strings)):
        if isempty(strings[i]):
            continue

        # if the string is a number, make it a float, otherwise leave it alone
        try:
            tokens[ntokens] = float(strings[i])
        except ValueError:
            tokens[ntokens] = strings[i]

        ntokens = ntokens + 1

    return [ntokens, tokens]
    # }}}
