import numpy as np

from checkfield import *
from fielddisplay import fielddisplay
from GetAreas import *
from project3d import *
from WriteData import *

class SMBarma(object):
    """SMBarma class definition

    Usage:
        SMBarma = SMBarma()
    """

    def __init__(self, *args):  # {{{
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.ar_order = 0
        self.ma_order = 0
        self.arlag_coefs = np.nan
        self.ma_order = 0
        self.malag_coefs = np.nan
        self.polynomialparams = np.nan
        self.datebreaks = np.nan
        self.basin_id = np.nan
        self.lapserates = np.nan
        self.elevationbins = np.nan
        self.refelevation = np.nan
        self.datebreaks = np.nan
        self.steps_per_step = 1
        self.averaging = 0
        self.requested_outputs = []

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'lapserates', 'basin-specific SMB lapse rates applied in each elevation bin, 1 row per basin, 1 column per bin, dimension 3 can be of size 12 to prescribe monthly varying values [m ice eq yr^-1 m^-1] (default: no lapse rate)'))
        s += '{}\n'.format(fielddisplay(self, 'elevationbins', 'basin-specific separations between elevation bins, 1 row per basin, 1 column per limit between bins, dimension 3 can be of size 12 to prescribe monthly varying values [m] (default: no basin separation)'))
        s += '{}\n'.format(fielddisplay(self, 'refelevation', 'basin-specific reference elevations at which SMB is calculated, and from which SMB is downscaled using lapserates (default: basin mean elevation) [m]'))
        s += '{}\n'.format(fielddisplay(self, 'steps_per_step', 'number of smb steps per time step'))
        s += '{}\n'.format(fielddisplay(self, 'averaging', 'averaging methods from short to long steps'))
        s += '\t\t{}\n'.format('0: Arithmetic (default)')
        s += '\t\t{}\n'.format('1: Geometric')
        s += '\t\t{}\n'.format('2: Harmonic')
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.ar_order = 0.0 # Autoregression model of order 0
        self.ma_order = 0.0 # Moving-average model of order 0
    # }}}

    def extrude(self, md):  # {{{
        return self # Nothing for now
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['SmbMassBalance']
    # }}}

    def initialize(self, md):  # {{{
        if self.ar_order == 0:
            self.ar_order = 1 # Dummy 1 value for autoregression
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression coefficients all set to 0
            print('      smb.ar_order (order of autoregressive model) not specified: order of autoregressive model set to 0')
        if self.ma_order == 0:
            self.ma_order = 1 # Dummy 1 value for moving-average
            self.malag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average coefficients all set to 0
            print('      smb.ma_order (order of moving-average model) not specified: order of moving-average model set to 0')
        if self.arma_timestep == 0:
            self.arma_timestep = md.timestepping.time_step # ARMA model has no prescribed time step
            print('      smb.arma_timestep (timestep of ARMA model) not specified: set to md.timestepping.time_step')
        if np.all(np.isnan(self.arlag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression model of order 0
            print('      smb.arlag_coefs (AR lag coefficients) not specified: order of autoregressive model set to 0')
        if np.all(np.isnan(self.malag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average model of order 0
            print('      smb.malag_coefs (MA lag coefficients) not specified: order of moving-average model set to 0')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        """
        TODO:
        - Ensure that checks on shape of self.lapserates are same as those under MATLAB as matrix addressing is quite different here
        """
        if 'MasstransportAnalysis' in analyses:
            nbas = md.smb.num_basins
            nprm = md.smb.num_params
            nbrk = md.smb.num_breaks
            md = checkfield(md, 'fieldname', 'smb.num_basins', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'smb.num_params', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'smb.num_breaks', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'smb.basin_id', 'Inf', 1, '>=', 0, '<=', md.smb.num_basins, 'size', [md.mesh.numberofelements])
            # if len(np.shape(self.polynomialparams)) == 1:
            #     self.polynomialparams = np.array([[self.polynomialparams]])
            if nbas > 1 and nbrk >= 1 and nprm > 1:
                md = checkfield(md, 'fieldname', 'smb.polynomialparams', 'NaN', 1, 'Inf', 1, 'size', [nbas, nbrk + 1, nprm], 'numel', nbas * (nbrk + 1) * nprm)
            elif nbas == 1:
                md = checkfield(md, 'fieldname', 'smb.polynomialparams', 'NaN', 1, 'Inf', 1, 'size', [nprm, nbrk + 1], 'numel', nbas * (nbrk + 1) * nprm)
            elif nbrk == 0:
                md = checkfield(md, 'fieldname', 'smb.polynomialparams', 'NaN', 1, 'Inf', 1, 'size', [nbas, nprm], 'numel', nbas * (nbrk + 1) * nprm)
            elif nprm == 1:
                md = checkfield(md, 'fieldname', 'smb.polynomialparams', 'NaN', 1, 'Inf', 1, 'size', [nbas, nbrk], 'numel', nbas * (nbrk + 1) * nprm)
            md = checkfield(md, 'fieldname', 'smb.ar_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'smb.ma_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'smb.arma_timestep', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', md.timestepping.time_step) # Autoregression time step cannot be finer than ISSM timestep
            md = checkfield(md, 'fieldname', 'smb.arlag_coefs', 'NaN', 1, 'Inf', 1, 'size', [nbas, md.smb.ar_order])
            md = checkfield(md, 'fieldname', 'smb.malag_coefs', 'NaN', 1, 'Inf', 1, 'size', [nbas, md.smb.ma_order])
            if nbrk > 0:
                md = checkfield(md, 'fieldname', 'smb.datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,nbrk])
            elif np.size(md.smb.datebreaks) == 0 or np.all(np.isnan(md.smb.datebreaks)):
                pass
            else:
                raise RuntimeError('md.smb.num_breaks is 0 but md.smb.datebreaks is not empty')

            if np.any(np.isnan(self.refelevation) is False) or np.size(self.refelevation) > 1:
                if len(np.shape(self.refelevation)) == 1:
                    self.refelevation = np.array([self.refelevation])
                md = checkfield(md, 'fieldname', 'smb.refelevation', 'NaN', 1, 'Inf', 1, '>=', 0, 'size', [1, nbas], 'numel', nbas)

            if (np.any(np.isnan(self.lapserates) is False) or np.size(self.lapserates) > 1):
                nbas = md.smb.num_basins
                if len(np.shape(self.lapserates)) == 1:
                    nbins = 1
                    self.lapserates = np.reshape(self.lapserates,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 2):
                    nbins = np.shape(self.lapserates)[1]
                    self.lapserates = np.reshape(self.lapserates,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 3):
                    nbins = np.shape(self.lapserates)[1]
                ntmlapse = np.shape(self.lapserates)[2]
                if len(np.shape(self.elevationbins)) < 3:
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,max(1,nbins-1),ntmlapse])
                md = checkfield(md, 'fieldname', 'smb.lapserates', 'NaN', 1, 'Inf', 1, 'size', [nbas,nbins,ntmlapse], 'numel', md.smb.num_basins*nbins*ntmlapse)
                md = checkfield(md, 'fieldname', 'smb.elevationbins', 'NaN', 1, 'Inf', 1, 'size', [nbas,max(1,nbins-1),ntmlapse], 'numel', nbas*max(1,nbins-1)*ntmlapse)
                for rr in range(nbas):
                    if(np.all(self.elevationbins[rr,0:-1]<=self.elevationbins[rr,1:])==False):
                        raise TypeError('md.smb.elevationbins should have rows in order of increasing elevation')
            elif (np.any(np.isnan(self.elevationbins) is False) or np.size(self.elevationbins) > 1):
                # Elevationbins specified but not lapserates: this will inevitably lead to inconsistencies
                nbas = md.smb.num_basins
                if len(np.shape(self.elevationbins)) == 1:
                    nbins = 1
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 2):
                    nbins = np.shape(self.elevationbins)[1]
                    self.elevationbins = np.reshape(self.elevationbins,[nbas,nbins,1])
                elif(len(np.shape(self.lapserates)) == 3):
                    nbins = np.shape(self.lapserates)[1]
                nbins = nbins - 1
                ntmlapse = np.shape(self.lapserates)[2]
                md = checkfield(md, 'fieldname', 'smb.lapserates', 'NaN', 1, 'Inf', 1, 'size', [nbas, nbins * ntmlapse], 'numel', nbas * nbins * ntmlapse)
                md = checkfield(md, 'fieldname', 'smb.elevationbins', 'NaN', 1, 'Inf', 1, 'size', [nbas, max(1, nbins - 1) * ntmlapse], 'numel', nbas * max(1, nbins - 1) * ntmlapse)

        md = checkfield(md, 'fieldname', 'smb.steps_per_step', '>=', 1, 'numel', [1])
        md = checkfield(md, 'fieldname', 'smb.averaging', 'numel', [1], 'values', [0, 1, 2])
        md = checkfield(md, 'fieldname', 'smb.requested_outputs', 'stringrow', 1)
        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        nbas = md.smb.num_basins
        nprm = md.smb.num_params
        nper = md.smb.num_breaks + 1
        if(np.any(np.isnan(md.smb.lapserates))):
            templapserates = np.zeros((nbas, 2, 12))
            print('      smb.lapserates not specified: set to 0')
            tempelevationbins = np.zeros((nbas, 1, 12)) # Dummy elevation bins
            nbins    = 2
            ntmlapse = 12
        else:
            if len(np.shape(md.smb.lapserates)) == 1:
                nbins    = 1
                ntmlapse = 1
            elif len(np.shape(md.smb.lapserates)) == 2:
                nbins    = np.shape(md.smb.lapserates)[1]
                ntmlapse = 1
            elif len(np.shape(md.smb.lapserates)) == 3:
                nbins    = np.shape(md.smb.lapserates)[1]
                ntmlapse = np.shape(md.smb.lapserates)[2]
            templapserates    = np.reshape(md.smb.lapserates,[nbas, nbins, ntmlapse])
            tempelevationbins = np.reshape(md.smb.elevationbins, [nbas, max(1, nbins - 1), ntmlapse])
        temprefelevation  = np.copy(md.smb.refelevation)
        # Scale the parameters
        polyparamsScaled   = np.copy(md.smb.polynomialparams)
        polyparams2dScaled = np.zeros((nbas, nper * nprm))
        if nprm > 1:
            # Case 3D
            if nbas > 1 and nper > 1:
                for ii in range(nprm):
                    polyparamsScaled[:, :, ii] = polyparamsScaled[:, :, ii] * (1 / yts) ** (ii + 1)
                # Fit in 2D array
                for ii in range(nprm):
                    polyparams2dScaled[:, ii * nper:(ii + 1) * nper] = 1 * polyparamsScaled[:, :, ii]
            # Case 2D and higher-order params at increasing row index
            elif nbas == 1:
                for ii in range(nprm):
                    polyparamsScaled[ii, :] = polyparamsScaled[ii, :] * (1 / yts) ** (ii + 1)
                # Fit in row array
                for ii in range(nprm):
                    polyparams2dScaled[0, ii * nper:(ii + 1) * nper] = 1 * polyparamsScaled[ii, :]
            # Case 2D and higher-order params at increasing column index
            elif nper == 1:
                for ii in range(nprm):
                    polyparamsScaled[:, ii] = polyparamsScaled[:, ii] * (1 / yts) ** (ii + 1)
                # 2D array is already in correct format
                polyparams2dScaled = np.copy(polyparamsScaled)
        else:
            polyparamsScaled   = polyparamsScaled * (1 / yts)
            # 2D array is already in correct format
            polyparams2dScaled = np.copy(polyparamsScaled)

        if nper == 1:
            dbreaks = np.zeros((nbas, 1))
        else:
            dbreaks = np.copy(md.smb.datebreaks)

        if ntmlapse == 1:
            templapserates    = np.repeat(templapserates, 12, axis = 2)
            tempelevationbins = np.repeat(tempelevationbins, 12, axis = 2)
        if np.any(np.isnan(md.smb.refelevation)):
            temprefelevation = np.zeros((nbas)).reshape(1, nbas)
            areas = GetAreas(md.mesh.elements, md.mesh.x, md.mesh.y)
            for ii, bid in enumerate(np.unique(md.smb.basin_id)):
                indices = np.where(md.smb.basin_id == bid)[0]
                elemsh = np.zeros((len(indices)))
                for jj in range(len(indices)):
                    elemsh[jj] = np.mean(md.geometry.surface[md.mesh.elements[indices[jj], :] - 1])
                temprefelevation[0, ii] = np.sum(areas[indices] * elemsh) / np.sum(areas[indices])
            if(np.any(templapserates != 0)):
                print('      smb.refelevation not specified: Reference elevations set to mean surface elevation of basins')
        nbins = np.shape(templapserates)[1]
        temp2dlapserates    = np.zeros((nbas, nbins * 12))
        temp2delevationbins = np.zeros((nbas, max(12, (nbins - 1) * 12)))
        for ii in range(12):
            temp2dlapserates[:, ii * nbins:(ii + 1) * nbins] = templapserates[:, :, ii]
            temp2delevationbins[:, ii * (nbins - 1):(ii + 1) * (nbins - 1)] = tempelevationbins[:, :, ii]

        WriteData(fid, prefix, 'name', 'md.smb.model', 'data', 13, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'num_basins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'num_breaks', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'num_params', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ar_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'ma_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'arma_timestep', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'basin_id', 'data', self.basin_id - 1, 'name', 'md.smb.basin_id', 'format', 'IntMat', 'mattype', 2)  # 0-indexed
        WriteData(fid, prefix, 'data', polyparams2dScaled, 'name', 'md.smb.polynomialparams', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'arlag_coefs', 'format', 'DoubleMat', 'name', 'md.smb.arlag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'smb', 'fieldname', 'malag_coefs', 'format', 'DoubleMat', 'name', 'md.smb.malag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'data', dbreaks, 'name', 'md.smb.datebreaks', 'format', 'DoubleMat','scale',yts)
        WriteData(fid, prefix, 'data', temp2dlapserates, 'name', 'md.smb.lapserates', 'format', 'DoubleMat', 'scale', 1 / yts, 'yts', yts)
        WriteData(fid, prefix, 'data', temp2delevationbins, 'name', 'md.smb.elevationbins', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'data', temprefelevation, 'name', 'md.smb.refelevation', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'data', nbins, 'name', 'md.smb.num_bins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'steps_per_step', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'averaging', 'format', 'Integer')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.smb.requested_outputs', 'format', 'StringArray')

    # }}}
