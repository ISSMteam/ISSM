from copy import deepcopy

import numpy as np

from fielddisplay import fielddisplay
from generic import generic
from getsubattr import getsubattr
from issmsettings import issmsettings
from meshintersect3d import meshintersect3d
from miscellaneous import miscellaneous
from model import model
from modelmerge3d import modelmerge3d
from pairoptions import pairoptions
from planetradius import planetradius
from plotmodel import plotmodel
from private import private
from setsubattr import setsubattr
from TwoDToThreeD import TwoDToThreeD


class sealevelmodel(object):
    """SEALEVELMODEL class definition

    Usage:
        slm = sealevelmodel(*args)

        where args is a variable list of options

    Example:
        slm = sealevelmodel(
            'icecap', md_greenland,
            'icecap', md_antarctica,
            'earth', md_earth
        )
    """

    def __init__(self, *args):  # {{{
        self.icecaps = [] # list of land/ice models; name should be changed later
        self.earth = 0 # model for the whole earth
        self.basins = []  # list of basins, matching icecaps, where shapefile info is held
        self.cluster = 0
        self.miscellaneous = 0
        self.settings = 0
        self.private = 0
        self.mergedcaps = 0
        self.transitions = []
        self.eltransitions = []
        self.planet = ''

        self.setdefaultparameters()

        if len(args):
            options = pairoptions(*args)

            # Recover all the icecap models
            self.icecaps = options.getfieldvalue('ice_cap', [])

            # Recover the earth models
            self.earth = options.getfieldvalue('earth', 0)

            # Set planet type
            self.planet = options.getfieldvalue('planet', 'earth')
    # }}}

    def __repr__(self):  # {{{
        s = '{}\n'.format(fielddisplay(self, 'icecaps', 'ice caps'))
        s += '{}\n'.format(fielddisplay(self, 'earth', 'earth'))
        s += '{}\n'.format(fielddisplay(self, 'settings', 'settings properties'))
        s += '{}\n'.format(fielddisplay(self, 'cluster', 'cluster parameters (number of cpus...'))
        s += '{}\n'.format(fielddisplay(self, 'miscellaneous', 'miscellaneous fields'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.icecaps = []
        self.earth = []
        self.cluster = generic()
        self.miscellaneous = miscellaneous()
        self.settings = issmsettings()
        self.private = private()
        self.transitions = []
        self.eltransitions = []
        self.planet = 'earth'
    # }}}

    @staticmethod
    def checkconsistency(slm, solutiontype):  # {{{
        # Is the coupler turned on?
        #for i in range(len(slm.icecaps)):
        #    if not slm.icecaps[i].transient.iscoupler:
        #        print('Warning: sealevelmodel.py::checkconsistency: icecap model {} should have the transient coupler option turned on!'.format(slm.icecaps[i].miscellaneous.name))

        #if not slm.earth.transient.iscoupler:
        #    print('Warning: sealevelmodel.py::checkconsistency: earth model should have the transient coupler option turned on!')

        # Check that the transition vectors have the right size
        if slm.earth.mesh.numberofvertices != len(slm.earth.solidearth.transfercount):
            raise Exception('sealevelmodel.py::checkconsistency: earth.solidearth.transfercount should be of size earth.mesh.numberofvertices')

        # Check that run frequency is the same everywhere
        for i in range(len(slm.icecaps)):
            if slm.icecaps[i].solidearth.settings.runfrequency != slm.earth.solidearth.settings.runfrequency:
                raise Exception('sealevelmodel.py::checkconsistency: icecap model {} should have the same run frequency as earth!'.format(slm.icecaps[i].miscellaneous.name))

        # Make sure steric_rate is the same everywhere
        for i in range(len(slm.icecaps)):
            md = slm.icecaps[i]
            if np.nonzero(md.dsl.sea_surface_height_above_geoid - slm.earth.dsl.sea_surface_height_above_geoid[slm.transitions[i]]) != []:
                raise Exception('sealevelmodel.py::checkconsistency: steric rate on ice cap {} is not the same as for the earth'.format(md.miscellaneous.name))

        # Make sure grd is the same everywhere
        for i in range(len(slm.icecaps)):
            md = slm.icecaps[i]
            if md.solidearth.settings.isgrd != slm.earth.solidearth.settings.isgrd:
                raise RuntimeError('sealevelmodel.py::checkconsistency: isgrd on ice cap {} is not the same as for the earth\n'.format(md.miscellaneous.name))

        # Make sure that there is no solid earth external forcing on the basins
        for i in range(len(slm.icecaps)):
            md = slm.icecaps[i]
            if md.solidearth.external:
                raise RuntimeError('sealevelmodel.py::checkconsistency: cannot run external forcings on an ice sheet when running a coupling earth/ice sheet model')
        # Make sure that we have the right grd model for computing our sealevel patterns
        for i in range(len(slm.icecaps)):
            md = slm.icecaps[i]
            if md.solidearth.settings.grdmodel != 0:
                raise RuntimeError('sealevelmodel.py::checkconsistency: ice sheets do not run GRD module, specify solidearth.settings.grdmodel=0 on ice cap {}'.format(i))
    # }}}

    def mergeresults(self):  # {{{
        champs = fieldnames(self.icecaps[0].results.TransientSolution)
        for i in range(len(self.mergedcaps / 2)):
            md = self.mergedcaps[2 * i]
            trans = self.mergedcaps[2 * i + 1]
            #icecaps = self.icecaps[self.range[2 * i + 2]]
            for j in range(len(self.icecaps[0].results.TransientSolution)):
                for k in range(len(champs)):
                    if isinstance(getattr(icecaps[0].results.TransientSolution[j], champs[k]), float):
                        # Vertex or element?
                        if len(getattr(icecaps[0].results.TransientSolution[j], champs[k]) == icecaps[0].mesh.numberofvertices):
                            setattr(md.results.TransientSolution[j], champs[k], np.zeros(md.mesh.numberofvertices))
                            for l in range(len(trans)):
                                resultcap = getattr(icecaps[l].results.TransientSolution[j], champs[k])
                                setattr(getattr(md.results.TransientSolution[j], champs[k]), trans[l], resultcap)
                        else:
                            if champs[k] == 'IceVolume' or champs[k] == 'IceVolumeAboveFlotation':
                                setattr(md.results.TransientSolution, champs[k], 0)
                                for l in range(len(trans)):
                                    resultcap = getattr(icecaps[l].results.TransientSolution[j], champs[k])
                                    setattr(md.results.TransientSolution[j], champs[k], getattr(md.results.TransientSolution[j], champs[k]) + resultcap)
                            elif champs[k] == 'time':
                                setattr(md.results.TransientSolution[j], champs[k], getattr(icecaps[0].results.TransientSolution[j], champs[k]))
                            else:
                                continue
                    else:
                        continue
            self.mergedcaps[2 * i] = md
    # }}}

    def listcaps(self):  # {{{
        for i in range(len(self.icecaps)):
            print('{}: {}'.format(i, self.icecaps[i].miscellaneous.name))
    # }}}

    def ncaps(self):  # {{{
        return len(self.icecaps)
    # }}}

    def continents(self):  # {{{
        list = []
        for i in range(len(self.basins)):
            list.append = self.basins[i].continent
        return np.unique(list)
    # }}}

    def basinsfromcontinent(self, continent):  # {{{
        list = []
        for i in range(len(self.icecaps)):
            if self.basins[i].continent == continent:
                list.append = self.basins[i].name
        return np.unique(list)
    # }}}

    def addbasin(self, bas):  # {{{
        if bas.__class__.__name__ != 'basin':
            raise Exception('addbasin method only takes a \'basin\' class object as input')
        self.basins.append(bas)
    # }}}

    def intersections2d(self, *args):  # {{{
        options = pairoptions(*args)
        force = options.getfieldvalue('force', 0)

        # Initialize, to avoid issues of having more transitions than meshes
        self.transitions = []
        self.eltransitions = []

        # For elements
        xe = np.mean(self.earth.mesh.x[self.earth.mesh.elements - 1], axis=1)
        ye = np.mean(self.earth.mesh.y[self.earth.mesh.elements - 1], axis=1)

        for i in range(len(self.icecaps)):
            mdi = self.icecaps[i]

            # For elements
            xei = np.mean(mdi.mesh.x[mdi.mesh.elements - 1], axis=1)
            yei = np.mean(mdi.mesh.y[mdi.mesh.elements - 1], axis=1)

            print('Computing vertex intersections for basin {}'.format(self.basins[i].name))

            self.transitions.append(meshintersect2d(self.earth.mesh.x, self.earth.mesh.y, mdi.mesh.x, mdi.mesh.y, 'force', force))
            self.eltransitions.append(meshintersect2d(xe, ye, xei, yei, 'force', force))
    # }}}

    def intersections(self, *args):  # {{{
        options = pairoptions(*args)
        force = options.getfieldvalue('force', 0)

        # Initialize, to avoid issues of having more transitions than meshes
        self.transitions = []
        self.eltransitions = []
        self.earth.solidearth.transfercount = np.zeros(self.earth.mesh.numberofvertices)

        # For elements
        xe = np.mean(self.earth.mesh.x[self.earth.mesh.elements - 1], axis=1)
        ye = np.mean(self.earth.mesh.y[self.earth.mesh.elements - 1], axis=1)
        ze = np.mean(self.earth.mesh.z[self.earth.mesh.elements - 1], axis=1)

        for i in range(len(self.icecaps)):
            mdi = self.icecaps[i]
            mdi = TwoDToThreeD(mdi, self.planet)

            # For elements
            xei = np.mean(mdi.mesh.x[mdi.mesh.elements - 1], axis=1)
            yei = np.mean(mdi.mesh.y[mdi.mesh.elements - 1], axis=1)
            zei = np.mean(mdi.mesh.z[mdi.mesh.elements - 1], axis=1)

            print('Computing vertex intersections for basin {}'.format(self.basins[i].name))

            self.transitions.append(meshintersect3d(self.earth.mesh.x, self.earth.mesh.y, self.earth.mesh.z, mdi.mesh.x, mdi.mesh.y, mdi.mesh.z, 'force', force))
            self.eltransitions.append(meshintersect3d(xe, ye, ze, xei, yei, zei, 'force', force))

            self.earth.solidearth.transfercount[self.transitions[i]] = self.earth.solidearth.transfercount[self.transitions[i]] + 1

        for i in range(len(self.icecaps)):
            self.icecaps[i].solidearth.transfercount = self.earth.solidearth.transfercount[self.transitions[i]]
    # }}}

    def checkintersections(self):  # {{{
        flags = np.zeros(self.earth.mesh.numberofvertices, 1)
        for i in range(len(self.basins)):
            flags[self.transitions[i]] = i
        plotmodel(self.earth, 'data', flags, 'coastlines', 'on')
    # }}}

    def checkbasinconsistency(self):  # {{{
        for i in range(len(self.basins)):
            self.basins[i].checkconsistency()
    # }}}

    def basinindx(self, *args):  # {{{
        options = pairoptions(*args)
        continent = options.getfieldvalue('continent', 'all')
        bas = options.getfieldvalue('basin', 'all')

        # Expand continent list #{{{
        if type(continent) == np.ndarray:
            if continent.shape[1] == 1:
                if continent[0] == 'all':
                    # Need to transform this into a list of continents
                    continent = []
                    for i in range(len(self.basins)):
                        continent.append(self.basins[i].continent)
                    continent = np.unique(continent)
                else:
                    pass  # Nothing to do: assume we have a list of continents
            else:
                pass  # Nothing to do: assume we have a list of continents
        else:
            if continent == 'all':
                # Need to transform this into a list of continents
                continent = []
                for i in range(len(self.basins)):
                    continent.append(self.basins[i].continent)
                continent = np.unique(continent)
            else:
                pass  # Nothing to do: assume we have a list of continents
        # }}}

        # Expand basins list using the continent list above and the extra bas discriminator #{{{
        if type(bas) == np.ndarray:
            if bas.shape[1] == 1:
                if bas[0] == 'all':
                    # Need to transform this into a list of basins
                    baslist = []
                    for i in range(len(self.basins)):
                        if self.basins[i].iscontinentany(continent):
                            baslist.append(i)
                    baslist = np.unique(baslist)
                else:
                    bas = bas[0]
                    baslist = []
                    for i in range(len(self.basins)):
                        if self.basins[i].iscontinentany(continent):
                            if self.basins[i].isnameany(bas):
                                baslist.append(i)
            else:
                # We have a list of basin names
                baslist = []
                for i in range(len(bas)):
                    basname = bas[i]
                    for j in range(len(self.basins)):
                        if self.basins[j].iscontinentany(continent):
                            if self.basins[j].isnameany(basname):
                                baslist.append(j)
                    baslist = np.unique(baslist)
        else:
            if bas == 'all':
                baslist = []
                for i in range(len(self.basins)):
                    if self.basins[i].iscontinentany(continent):
                        baslist.append(i)
                baslist = np.unique(baslist)
            else:
                baslist = []
                for i in range(len(self.basins)):
                    if self.basins[i].iscontinentany(continent):
                        if self.basins[i].isnameany(bas):
                            baslist.append(i)
                baslist = np.unique(baslist)

        return baslist
        # }}}
    # }}}

    def addicecap(self, md):  # {{{
        if not type(md) == model:
            raise Exception("addicecap method only takes a 'model' class object as input")

        self.icecaps.append(md)
    # }}}

    def basinsplot3d(self, *args):  # {{{
        for i in range(len(self.basins)):
            self.basins[i].plot3d(*args)
    # }}}

    def caticecaps(self, *args):  # {{{
        # Recover options
        options = pairoptions(*args)
        tolerance = options.getfieldvalue('tolerance', .65)
        loneedgesdetect = options.getfieldvalue('loneedgesdetect', 0)

        # Make 3D model
        models = deepcopy(self.icecaps)
        for i in range(len(models)):
            models[i] = TwoDToThreeD(models[i], self.planet)

        # Plug all models together
        md = models[0]
        for i in range(1, len(models)):
            md = modelmerge3d(md, models[i], 'tolerance', tolerance)
            md.private.bamg.landmask = np.hstack((md.private.bamg.landmask, models[i].private.bamg.landmask))

        # Look for lone edges if asked for it
        if loneedgesdetect:
            edges = loneedges(md)
            # TODO: Reconfigure the following in the process of bringing plotting online
            plotmodel(md, 'data', md.mask.land_levelset)
            for i in range(len(edges)):
                ind1 = edges(i, 1)
                ind2 = edges(i, 2)
                plot3([md.mesh.x[ind1], md.mesh.x[ind2]], [md.mesh.y[ind1], md.mesh.y[ind2]], [md.mesh.z[ind1], md.mesh.z[ind2]], 'g*-')

        # Plug into earth
        self.earth = md

        # Create mesh radius
        self.earth.mesh.r = planetradius('earth') * np.ones((md.mesh.numberofvertices, ))
    # }}}

    def caticecaps2d(self, *args):  # {{{
        # Recover options
        options = pairoptions(*args)
        tolerance = options.getfieldvalue('tolerance', 1e-5)
        loneedgesdetect = options.getfieldvalue('loneedgesdetect', 0)
        models = self.icecaps

        # Plug all models together
        md = models[0]
        for i in range(1, len(models)):
            md = modelmerge2d(md, models[i], 'tolerance', tolerance)

        # Look for lone edges if asked for it
        if loneedgesdetect:
            edges = loneedges(md)
            # TODO: Reconfigure the following in the process of bringing plotting online
            plotmodel(md, 'data', md.mask.land_levelset)
            for i in range(len(edges)):
                ind1 = edges(i, 1)
                ind2 = edges(i, 2)
                plot([md.mesh.x[ind1], md.mesh.x[ind2]], [md.mesh.y[ind1], md.mesh.y[ind2]], 'g*-')

        # Plug into earth
        self.earth = md
    # }}}

    def viscousiterations(self):  # {{{
        for i in range(len(self.icecaps)):
            ic = self.icecaps[i]
            mvi = ic.results.TransientSolution[0].StressbalanceConvergenceNumSteps
            for j in range(1, len(ic.results.TransientSolution) - 1):
                mvi = np.max(mvi, ic.results.TransientSolution[j].StressbalanceConvergenceNumSteps)
            print("{}, {}: {}".format(i, self.icecaps[i].miscellaneous.name, mvi))
    # }}}

    def maxtimestep(self):  # {{{
        for i in range(len(self.icecaps)):
            ic = self.icecaps[i]
            mvi = len(ic.results.TransientSolution)
            timei = ic.results.TransientSolution[-1].time
            print("{}, {}: {}/{}".format(i, self.icecaps[i].miscellaneous.name, mvi, timei))

        mvi = len(self.earth.results.TransientSolution)
        timei = self.earth.results.TransientSolution[-1].time
        print("Earth: {}/{}", mvi, timei)
    # }}}

    def transfer(self, string):  # {{{
        # Recover field size in one icecap
        n = getsubattr(self.icecaps[0], string).shape[0]

        if n == self.icecaps[0].mesh.numberofvertices:
            setsubattr(self.earth, string, np.zeros((self.earth.mesh.numberofvertices, ))) # Assign array of zeros to target attribute
            earth_attr = getsubattr(self.earth, string) # Retrieve reference to target attribute
            for i in range(len(self.icecaps)):
                earth_attr[self.transitions[i]] = getsubattr(self.icecaps[i], string)
        elif n == (self.icecaps[0].mesh.numberofvertices + 1):
            # Dealing with transient dataset
            # Check that all timetags are similar between all icecaps #{{{
            for i in range(len(self.icecaps)):
                capfieldi = getsubattr(self.icecaps[i], string)
                for j in range(i + 1, len(self.icecaps)):
                    capfieldj = getsubattr(self.icecaps[j], string)
                    if capfieldi[-1, :] != capfieldj[-1, :]:
                        raise Exception("Time stamps for {} field are different between icecaps {} and {}".format(string, i, j))
            capfield1 = getsubattr(self.icecaps[0], string)
            times = capfield1[-1, :]
            nsteps = len(times)
            # }}}
            # Initialize #{{{
            field = np.zeros((self.earth.mesh.numberofvertices + 1, nsteps))
            field[-1, :] = times  # Transfer the times only, not the values
            # }}}
            # Transfer all the time fields #{{{
            for i in range(len(self.icecaps)):
                capfieldi = getsubattr(self.icecaps[i], string)
                for j in range(nsteps):
                    field[self.transitions[i], j] = capfieldi[0:-1, j]  # Transfer only the values, not the time
            setsubattr(self.earth, string, field)  # Do not forget to plug the field variable into its final location
            # }}}
        elif n == (self.icecaps[0].mesh.numberofelements):
            setsubattr(self.earth, string, np.zeros((self.earth.mesh.numberofelements, ))) # Assign array of zeros to target attribute
            earth_attr = getsubattr(self.earth, string) # Retrieve reference to target attribute
            for i in range(len(self.icecaps)):
                earth_attr[self.eltransitions[i]] = getsubattr(self.icecaps[i], string)
        else:
            raise Exception('not supported yet')
    # }}}

    def homogenize(self, noearth=0):  # {{{
        mintimestep = np.inf

        for i in range(len(self.icecaps)):
            ic = self.icecaps[i]
            mintimestep = np.min(mintimestep, len(ic.results.TransientSolution))

        if not noearth:
            mintimestep = np.min(mintimestep, len(self.earth.results.TransientSolution))

        for i in range(len(self.icecaps)):
            ic = self.icecaps[i]
            ic.resuts.TransientSolution = ic.results.TransientSolution[:mintimestep]
            self.icecaps[i] = ic

        ic = self.earth

        if not noearth:
            ic.results.TransientSolution = ic.results.TransientSolution[:mintimestep]

        self.earth = ic

        return self
    # }}}

    def initializemodels(self):  # {{{
        for i in range(len(self.basins)):
            md = model()
            md.miscellaneous.name = self.basins[i].name
            self.addicecap(md)
    # }}}
