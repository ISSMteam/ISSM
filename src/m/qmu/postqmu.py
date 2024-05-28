from copy import deepcopy
from os import getpid, stat
from os.path import isfile
from subprocess import call

from dakota_out_parse import *
from helpers import *
import loadresultsfromdisk as loadresultsfromdisk # There is a name conflict somewhere
from results import results, resultsdakota, solution


def postqmu(md):
    """POSTQMU - Deal with Dakota output results in files
    
    Usage:
        md = postqmu(md)

    TODO:
    - Add checks to Popen
    """

    # check to see if dakota returned errors in the err file
    qmuerrfile = str(md.miscellaneous.name) + '.qmu.err'

    if isfile(qmuerrfile) and stat(qmuerrfile).st_size > 0:
        with open(qmuerrfile, 'r') as fide:
            fline = fide.read()
            print(fline)

        raise RuntimeError('Dakota returned error in ' + str(qmuerrfile) + ' file.')

    # parse inputs and results from dakota
    qmuinfile = str(md.miscellaneous.name) + '.qmu.in'
    qmuoutfile = str(md.miscellaneous.name) + '.qmu.out'
    [method, dresp_out, scm, pcm, srcm, prcm] = dakota_out_parse(qmuoutfile)
    dakotaresults = resultsdakota()
    dakotaresults.dresp_out = dresp_out
    dakotaresults.scm       = scm
    dakotaresults.pcm       = pcm
    dakotaresults.srcm      = srcm
    dakotaresults.prcm      = prcm

    if isfile('dakota_tabular.dat'):
        # only need a subset of the outputs; dakota_out_parse handles .dat seperately
        [method, dresp_dat, _, _, _, _] = dakota_out_parse('dakota_tabular.dat')
        dakotaresults.dresp_dat = dresp_dat

    if md.qmu.output and md.qmu.statistics.method[0]['name'] == 'None':
        if md.qmu.method.method == 'nond_sampling':
            dakotaresults.modelresults = []
            md2 = deepcopy(md)
            md2.results = results()
            md2.qmu.isdakota = 0
            for i in range(md2.qmu.method.params.samples):
                outbin_name = '{}.outbin.{}'.format(md2.miscellaneous.name, (i + 1))
                print('Reading qmu file {}'.format(outbin_name))
                md2 = loadresultsfromdisk.loadresultsfromdisk(md2, outbin_name)
                dakotaresults.modelresults.append(deepcopy(md2.results))
            del md2

    if md.qmu.statistics.method[0]['name'] != 'None':
        md.qmu.isdakota = 0
        md = loadresultsfromdisk.loadresultsfromdisk(md, '{}.stats'.format(md.miscellaneous.name))
        md.qmu.isdakota = 1

    # put dakotaresults in their right location.
    md.results.dakota = deepcopy(dakotaresults)

    # move all the individual function evalutations into zip files
    if not md.qmu.isdakota:
        subproc_args = 'zip -mq params.in.zip params.in.[1-9]*'
        call(subproc_args, shell=True)
        subproc_args = 'zip -mq results.out.zip results.out.[1-9]*'
        call(subproc_args, shell=True)
        subproc_args = 'zip -mq matlab.out.zip matlab*.out.[1-9]*'
        call(subproc_args, shell=True)

    return md
