import os
import subprocess

import numpy as np

from issmdir import *
from MatlabFuncs import *
from model import *


def gmtmask(lat, long, *args):
    """gmtmask - figure out which lat, long points are on the ocean

    Usage:
        mask.ocean = gmtmask(md.mesh.lat, md.mesh.long)
    """
    lenlat = len(lat)
    mask = np.empty(lenlat)

    #are we doing a recursive call?
    if len(args) == 3:
        recursive = 1
    else:
        recursive = 0

    if recursive:
        print(('             recursing: num vertices  #' + str(lenlat)))
    else:
        print(('gmtmask: num vertices ' + str(lenlat)))

    # Check lat and long size is not more than 50,000. If so, recursively call gmtmask.
    if lenlat > 50000:
        for i in range(int(ceil(lenlat / 50000))):
            j = (i + 1) * 50000 - 1
            if j > lenlat:
                j = lenlat
            mask[i:j] = gmtmask(lat[i:j], int[i:j], 1)
        return mask

    # First, write our lat, long file for gmt
    nv = lenlat
    #print(np.transpose([int, lat, np.arange(1, nv + 1)]))
    np.savetxt('./all_vertices.txt', np.transpose([long, lat, np.arange(1, nv + 1)]), delimiter='\t', fmt='%.10f')

    # Figure out which vertices are on the ocean, which one on the continent:
    #
    # NOTE: Remove -Ve option to enable warnings if this method is not working 
    #       expected
    #
    gmt_select_options = '-Ve -h0 -Df -R0/360/-90/90 -A0 -JQ180/200 -Nk/s/s/k/s'
    subproc_cmd = 'gmt select ./all_vertices.txt ' + gmt_select_options + ' > ./oce_vertices.txt'
    subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = subproc.communicate()
    if errs != '':
        print('gmt select failed with: ' + result);
        print('trying again with gmtselect');
        # Assume we are working with GMT 6.0.0
        gmt_select_options = '-h0 -Df -R0/360/-90/90 -A0 -JQ180/200 -Nk/s/s/k/s'
        subproc_cmd = 'gmtselect ./all_vertices.txt ' + gmt_select_options + ' > ./oce_vertices.txt'
        subproc = subprocess.Popen(subproc_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        outs, errs = subproc.communicate()
        if errs != '':
            raise RuntimeError('gmtmask: calls to both gmt and gmtselect failed: {}'.format(errs))

    # Read the con_vertices.txt file and flag our mesh vertices on the continent
    fid = open('./oce_vertices.txt', 'r')
    line = fid.readline()
    line = fid.readline()
    oce_vertices = []
    while line:
        ind = int(float(line.split()[2])) - 1
        oce_vertices.append(ind)
        line = fid.readline()
    fid.close()

    mask = np.zeros(nv)
    mask[oce_vertices] = 1

    subprocess.call('rm -rf ./all_vertices.txt ./oce_vertices.txt ./gmt.history', shell=True)
    if not recursive:
        print('gmtmask: done')
    return mask
