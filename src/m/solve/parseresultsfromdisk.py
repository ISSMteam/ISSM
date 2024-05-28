from collections import OrderedDict
import struct
import numpy as np
from results import solution


def parseresultsfromdisk(md, filename, iosplit):  #{{{
    if iosplit:
        saveres = parseresultsfromdiskiosplit(md, filename)
    else:
        saveres = parseresultsfromdiskioserial(md, filename)
        #saveres = parseresultsfromdiskioserialsequential(md, filename)
    return saveres
# }}}


def parseresultsfromdiskiosplit(md, filename):  # {{{
    #Open file
    try:
        fid = open(filename, 'rb')
    except IOError:
        raise IOError("parseresultsfromdisk error message: could not open '{}' for binary reading.".format(filename))

    saveres = []

    #if we have done split I/O, ie, we have results that are fragmented across patches,
    #do a first pass, and figure out the structure of results
    loadres = ReadDataDimensions(fid)
    while loadres:
        # Get time and step
        if loadres['step'] > len(saveres):
            for i in range(len(saveres), loadres['step'] - 1):
                saveres.append(None)
            saveres.append(resultsclass.results())
        setattr(saveres[loadres['step'] - 1], 'step', loadres['step'])
        setattr(saveres[loadres['step'] - 1], 'time', loadres['time'])

        # Add result
        setattr(saveres[loadres['step'] - 1], loadres['fieldname'], float('NaN'))

        # Read next result
        loadres = ReadDataDimensions(fid)

    #do a second pass, and figure out the size of the patches
    fid.seek(0)  #rewind
    loadres = ReadDataDimensions(fid)
    while loadres:

        #read next result
        loadres = ReadDataDimensions(fid)

    #third pass, this time to read the real information
    fid.seek(0)  #rewind
    loadres = ReadData(fid, md)
    while loadres:

        #Get time and step
        if loadres['step'] > len(saveres):
            for i in range(len(saveres), loadres['step'] - 1):
                saveres.append(None)
            saveres.append(saveresclass.saveres())
        setattr(saveres[loadres['step'] - 1], 'step', loadres['step'])
        setattr(saveres[loadres['step'] - 1], 'time', loadres['time'])

        # Add result
        setattr(saveres[loadres['step'] - 1], loadres['fieldname'], loadres['field'])

        # Read next result
        loadres = ReadData(fid, md)

    # Close file
    fid.close()

    return saveres
# }}}

def parseresultsfromdiskioserial(md, filename):  # {{{
    # Open file
    try:
        fid = open(filename, 'rb')
    except IOError as e:
        raise IOError('parseresultsfromdisk error message: could not open {} for binary reading'.format(filename))

    # Collect all results in a list
    allresults = []
    while True:
        # Read next result
        result = ReadData(fid, md)

        if result is None:
            if allresults == []:
                raise Exception('no results found in binary file ' + filename)
            else:
                break

        allresults.append(result)
    fid.close()

    # Now, process all results and find out how many steps we have
    numresults = len(allresults)
    allsteps = np.zeros((numresults, 1))
    for i in range(numresults):
        allsteps[i] = allresults[i]['step']
    pos = np.where(allsteps != -9999)
    allsteps = np.sort(np.unique(allsteps[pos]))

    # Ok, now construct structure
    results = solution()

    for i in range(numresults):
        result = allresults[i]
        index = 0
        if result['step'] != -9999:
            index = np.where(result['step'] == allsteps)[0][0]
            setattr(results[index], 'step', result['step'])
        if result['time'] != -9999:
            setattr(results[index], 'time', result['time'])
        setattr(results[index], result['fieldname'], result['field'])
    return results
# }}}

def ReadData(fid, md):  # {{{
    """READDATA

    Usage:
        field = ReadData(fid, md)
    """

    # Read field
    try:
        length = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        fieldname = struct.unpack('{}s'.format(length), fid.read(length))[0][:-1]
        fieldname = fieldname.decode() # strings are binaries when stored so need to be converted back
        time = struct.unpack('d', fid.read(struct.calcsize('d')))[0]
        step = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        datatype = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        M = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        if datatype == 1:
            field = np.array(struct.unpack('{}d'.format(M), fid.read(M * struct.calcsize('d'))), dtype=float)

        elif datatype == 2:
            field = struct.unpack('{}s'.format(M), fid.read(M))[0][:-1]
            field = field.decode()

        elif datatype == 3:
            N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            #field = transpose(fread(fid, [N M], 'double'))
            field = np.zeros(shape=(M, N), dtype=float)
            for i in range(M):
                field[i, :] = struct.unpack('{}d'.format(N), fid.read(N * struct.calcsize('d')))

        elif datatype == 4:
            N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            #field = transpose(fread(fid, [N M], 'int'))
            field = np.zeros(shape=(M, N), dtype=int)
            for i in range(M):
                field[i, :] = struct.unpack('{}i'.format(N), fid.read(N * struct.calcsize('i')))

        else:
            raise TypeError('cannot read data of datatype {}'.format(datatype))

        # Process units here FIXME: this should not be done here!
        yts = md.constants.yts
        if fieldname == 'BalancethicknessThickeningRate':
            field = field * yts
        elif fieldname == 'HydrologyWaterVx':
            field = field * yts
        elif fieldname == 'HydrologyWaterVy':
            field = field * yts
        elif fieldname == 'Vx':
            field = field * yts
        elif fieldname == 'Vy':
            field = field * yts
        elif fieldname == 'Vz':
            field = field * yts
        elif fieldname == 'Vel':
            field = field * yts
        elif fieldname == 'VxShear':
            field = field * yts
        elif fieldname == 'VyShear':
            field = field * yts
        elif fieldname == 'VxBase':
            field = field * yts
        elif fieldname == 'VyBase':
            field = field * yts
        elif fieldname == 'VxSurface':
            field = field * yts
        elif fieldname == 'VySurface':
            field = field * yts
        elif fieldname == 'VxAverage':
            field = field * yts
        elif fieldname == 'VyAverage':
            field = field * yts
        elif fieldname == 'VxDebris':
            field = field * yts
        elif fieldname == 'VyDebris':
            field = field * yts
        elif fieldname == 'BasalforcingsGroundediceMeltingRate':
            field = field * yts
        elif fieldname == 'BasalforcingsFloatingiceMeltingRate':
            field = field * yts
        elif fieldname == 'BasalforcingsSpatialDeepwaterMeltingRate':
            field = field * yts
        elif fieldname == 'BasalforcingsSpatialUpperwaterMeltingRate':
            field = field * yts
        elif fieldname == 'TotalFloatingBmb':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalFloatingBmbScaled':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalGroundedBmb':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalGroundedBmbScaled':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalSmb':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalSmbScaled':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalSmbMelt':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'TotalSmbRefreeze':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'GroundinglineMassFlux':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'IcefrontMassFlux':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'IcefrontMassFluxLevelset':
            field = field / pow(10.0, 12) * yts # (GigaTon/year)
        elif fieldname == 'SmbMassBalance':
            field = field * yts
        elif fieldname == 'SmbPrecipitation':
            field = field * yts
        elif fieldname == 'SmbRain':
            field = field * yts
        elif fieldname == 'SmbRunoff':
            field = field * yts
        elif fieldname == 'SmbRunoffSubstep':
            field = field * yts
        elif fieldname == 'SmbEvaporation':
            field = field * yts
        elif fieldname == 'SmbRefreeze':
            field = field * yts
        elif fieldname == 'SmbEC':
            field = field * yts
        elif fieldname == 'SmbAccumulation':
            field = field * yts
        elif fieldname == 'SmbMelt':
            field = field * yts
        elif fieldname == 'SmbMAdd':
            field = field * yts
        elif fieldname == 'SmbWAdd':
            field = field * yts
        elif fieldname == 'CalvingCalvingrate':
            field = field * yts
        elif fieldname == 'Calvingratex':
            field = field * yts
        elif fieldname == 'Calvingratey':
            field = field * yts
        elif fieldname == 'CalvingMeltingrate':
            field = field * yts
        elif fieldname == 'LoveKernelsReal' or fieldname == 'LoveKernelsImag':
            nlayer = md.materials.numlayers
            degmax = md.love.sh_nmax
            nfreq = md.love.nfreq
            r0 = md.love.r0
            g0 = md.love.g0
            mu0 = md.love.mu0
            rr = md.materials.radius
            rho = md.materials.density
            rho_avg_partial = np.diff(np.power(rr, 3), n=1, axis=0)
            rho_avg = ((rho * rho_avg_partial) / rho_avg_partial.sum()).sum()
            temp_field = np.empty((degmax + 1, nfreq, nlayer + 1, 6))
            temp_field.fill(0.0)
            for ii in range(degmax + 1):
                for jj in range(nfreq):
                    for kk in range(nlayer + 1):
                        if kk < nlayer: # NOTE: Upper bound of range is non-inclusive (compare to src/m/solve/parseresultsfromdisk.m)
                            ll = ii * (nlayer + 1) * 6 + (kk * 6 + 1) + 3
                            temp_field[ii, jj, kk, 0] = field[ll + (0 - 1), jj] * r0        # mm = 4
                            temp_field[ii, jj, kk, 1] = field[ll + (1 - 1), jj] * mu0       # mm = 5
                            temp_field[ii, jj, kk, 2] = field[ll + (2 - 1), jj] * r0        # mm = 6
                            temp_field[ii, jj, kk, 3] = field[ll + (3 - 1), jj] * mu0       # mm = 1
                            temp_field[ii, jj, kk, 4] = field[ll + (4 - 1), jj] * r0 * g0   # mm = 2
                            temp_field[ii, jj, kk, 5] = field[ll + (5 - 1), jj] * g0        # mm = 3
                        else: # surface
                            ll = (ii + 1) * (nlayer + 1) * 6 - 2
                            temp_field[ii, jj, kk, 0] = field[ll + (0 - 1), jj] * r0
                            temp_field[ii, jj, kk, 2] = field[ll + (1 - 1), jj] * r0
                            temp_field[ii, jj, kk, 4] = field[ll + (2 - 1), jj] * r0 * g0
                            # surface BC
                            temp_field[ii, jj, kk, 3] = 0
                            if md.love.forcing_type == 9:
                                temp_field[ii, jj, kk, 1] = 0
                                temp_field[ii, jj, kk, 5] = (2 * (ii + 1) - 1) / r0 - (ii + 1) * field[ll + (2 - 1), jj] * g0
                            elif md.love.forcing_type == 11:
                                temp_field[ii, jj, kk, 1] = -(2 * ii + 1) * rho_avg / 3
                                temp_field[ii, jj, kk, 5] = (2 * (ii + 1) - 1) / r0 - (ii + 1) * field[ll + (2 - 1), jj] * g0
            field = temp_field

        if time != -9999:
            time = time / yts

        saveres = OrderedDict()
        saveres['fieldname'] = fieldname
        saveres['time'] = time
        saveres['step'] = step
        saveres['field'] = field

    except struct.error as e:
        saveres = None

    return saveres
# }}}


def ReadDataDimensions(fid):  # {{{
    """READDATADIMENSIONS - read data dimensions, step and time, but not the data itself.

    Usage:
        field = ReadDataDimensions(fid)
    """

    # Read field
    try:
        length = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        fieldname = struct.unpack('{}s'.format(length), fid.read(length))[0][:-1]
        time = struct.unpack('d', fid.read(struct.calcsize('d')))[0]
        step = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        datatype = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        M = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
        N = 1 # default
        if datatype == 1:
            fid.seek(M * 8, 1)
        elif datatype == 2:
            fid.seek(M, 1)
        elif datatype == 3:
            N = struct.unpack('i', fid.read(struct.calcsize('i')))[0]
            fid.seek(N * M * 8, 1)
        else:
            raise TypeError("cannot read data of datatype {}".format(datatype))

        saveres = OrderedDict()
        saveres['fieldname'] = fieldname
        saveres['time'] = time
        saveres['step'] = step
        saveres['M'] = M
        saveres['N'] = N

    except struct.error as Err:
        print(Err)
        saveres = None

    return saveres
# }}}

def addfieldtorecord(a, descr):  # {{{
    if a.dtype.fields is None:
        raise ValueError('\'a\' must be a structured numpy array')
    b = np.empty(a.shape, dtype=a.dtype.descr + descr)
    for name in a.dtype.names:
        b[name] = a[name]

    return b
# }}}
