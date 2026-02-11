#!/usr/bin/env python3
import matplotlib.pyplot as plt

def plot_none(md, options, fig, axgrid, gridindex):
    '''
    PLOT_NONE - plot nothing, just apply options

    Usage:
        plot_none(md,options,fig,axgrid,gridindex)

    See also: ??
    '''
    return

    options=options.addfielddefault('colorbar','none')
    options=options.addfielddefault('map','none')
    #FIXME: What do you mean 'axis" in matlab? Check 'axis' option in python.
    #options=options.addfielddefault('axis','equal')

    #TODO: overlay option in "plot_none.m"

    #apply options
    applyoptions(md,[],options);
