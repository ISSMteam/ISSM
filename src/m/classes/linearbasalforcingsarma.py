import numpy as np

from checkfield import *
from fielddisplay import fielddisplay
from project3d import *
from WriteData import *

class linearbasalforcingsarma(object):
    """linearbasalforcingsarma class definition

    Usage:
        linearbasalforcingsarma = linearbasalforcingsarma()
    """

    def __init__(self, *args):  # {{{
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.datebreaks       = np.nan
        self.ar_order = 0
        self.ma_order = 0
        self.arma_timestep = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.basin_id = np.nan
        self.groundedice_melting_rate = np.nan
        self.deepwater_elevation = np.nan
        self.upperwater_melting_rate = np.nan
        self.upperwater_elevation = np.nan
        self.geothermalflux = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   surface forcings parameters:\n'
        s += '   autoregressive model is applied for deepwater_melting_rate\n'
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
        s += '{}\n'.format(fielddisplay(self, 'deepwater_elevation', 'basin-specific elevation of ocean deepwater [m]'))
        s += '{}\n'.format(fielddisplay(self, 'upperwater_melting_rate', 'basin-specic basal melting rate (positive if melting applied for floating ice whith base >= upperwater_elevation) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'upperwater_elevation', 'basin-specific elevation of ocean upperwater [m]'))
        s += '{}\n'.format(fielddisplay(self, 'groundedice_melting_rate','node-specific basal melting rate (positive if melting) [m/yr]'))
        s += '{}\n'.format(fielddisplay(self, 'geothermalflux','node-specific geothermal heat flux [W/m^2]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.ar_order = 0.0 # Autoregression model of order 0
        self.ma_order = 0.0 # Moving-average model of order 0
        return self
    # }}}

    def extrude(self, md):  # {{{
        return self # Nothing for now
    # }}}

    def initialize(self, md):  # {{{
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices))
            print("      no basalforcings.groundedice_melting_rate specified: values set as zero")
        if np.all(np.isnan(self.trend)):
            self.trend = np.zeros((1, self.num_basins)) # No trend in SMB
            print('      basalforcings.trend (trend) not specified: value set to 0')
        if self.ar_order == 0:
            self.ar_order = 1 # Dummy 1 value for autoregression
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression coefficients all set to 0
            print('      basalforcings.ar_order (order of autoregressive model) not specified: order of autoregressive model set to 0')
        if self.arma_timestep == 0:
            self.arma_timestep = md.timestepping.time_step # ARMA model has no prescribed time step
            print('      basalforcings.arma_timestep (timestep of ARMA model) not specified: set to md.timestepping.time_step')
        if np.all(np.isnan(self.arlag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ar_order)) # Autoregression model of order 0
            print('      basalforcings.arlag_coefs (AR lag coefficients) not specified: order of autoregressive model set to 0')
        if np.all(np.isnan(self.malag_coefs)):
            self.arlag_coefs = np.zeros((self.num_basins, self.ma_order)) # Moving-average model of order 0
            print('      basalforcings.malag_coefs (MA lag coefficients) not specified: order of moving-average model set to 0')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if 'MasstransportAnalysis' in analyses:
            nbas = md.basalforcings.num_basins;
            nprm = md.basalforcings.num_params;
            nbrk = md.basalforcings.num_breaks;
            md = checkfield(md, 'fieldname', 'basalforcings.num_basins', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.groundedice_melting_rate', 'NaN', 1, 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'basalforcings.num_params', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.num_breaks', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)

            if len(np.shape(self.deepwater_elevation)) == 1:
                self.deepwater_elevation = np.array([self.deepwater_elevation])
                self.upperwater_elevation = np.array([self.upperwater_elevation])
                self.upperwater_melting_rate = np.array([self.upperwater_melting_rate])
            if len(np.shape(self.polynomialparams)) == 1:
                self.polynomialparams = np.array([[self.polynomialparams]])
            if(nbas>1 and nbrk>=1 and nprm>1):
                md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm)
            elif(nbas==1):
                md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm)
            elif(nbrk==0):
                md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm)
            elif(nprm==1):
                md = checkfield(md,'fieldname','basalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk],'numel',nbas*(nbrk+1)*nprm)
            md = checkfield(md, 'fieldname', 'basalforcings.deepwater_elevation', 'NaN', 1, 'Inf', 1, 'size', [1, md.basalforcings.num_basins], 'numel', md.basalforcings.num_basins)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_elevation', 'NaN', 1, 'Inf', 1, '<=', 0, 'size', [1, md.basalforcings.num_basins], 'numel', md.basalforcings.num_basins)
            md = checkfield(md, 'fieldname', 'basalforcings.upperwater_melting_rate', 'NaN', 1, 'Inf', 1,'>=', 0, 'size', [1, md.basalforcings.num_basins], 'numel', md.basalforcings.num_basins)
            md = checkfield(md, 'fieldname', 'basalforcings.basin_id', 'Inf', 1, '>=', 0, '<=', md.basalforcings.num_basins, 'size', [md.mesh.numberofelements])


            md = checkfield(md, 'fieldname', 'basalforcings.ar_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.ma_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'basalforcings.arma_timestep', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', md.timestepping.time_step) # Autoregression time step cannot be finer than ISSM timestep
            md = checkfield(md, 'fieldname', 'basalforcings.arlag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.basalforcings.num_basins, md.basalforcings.ar_order])
            md = checkfield(md, 'fieldname', 'basalforcings.malag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.basalforcings.num_basins, md.basalforcings.ma_order])
            if(nbrk>0):
                md = checkfield(md, 'fieldname', 'basalforcings.datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,nbrk])
            elif(np.size(md.basalforcings.datebreaks)==0 or np.all(np.isnan(md.basalforcings.datebreaks))):
                pass
            else:
                raise RuntimeError('md.basalforcings.num_breaks is 0 but md.basalforcings.datebreaks is not empty')

        if 'BalancethicknessAnalysis' in analyses:
            raise Exception('not implemented yet!')
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            raise Exception('not implemented yet!')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        nbas = md.basalforcings.num_basins;
        nprm = md.basalforcings.num_params;
        nper = md.basalforcings.num_breaks+1;
        # Scale the parameters #
        polyparamsScaled   = np.copy(md.basalforcings.polynomialparams)
        polyparams2dScaled = np.zeros((nbas,nper*nprm))
        if(nprm>1):
            # Case 3D #
            if(nbas>1 and nper>1):
                for ii in range(nprm):
                    polyparamsScaled[:,:,ii] = polyparamsScaled[:,:,ii]*(1/yts)**(ii+1)
                # Fit in 2D array #
                for ii in range(nprm):
                    polyparams2dScaled[:,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[:,:,ii]
            # Case 2D and higher-order params at increasing row index #
            elif(nbas==1):
                for ii in range(nprm):
                    polyparamsScaled[ii,:] = polyparamsScaled[ii,:]*(1/yts)**(ii+1)
                # Fit in row array #
                for ii in range(nprm):
                    polyparams2dScaled[0,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[ii,:]
            # Case 2D and higher-order params at incrasing column index #
            elif(nper==1):
                for ii in range(nprm):
                    polyparamsScaled[:,ii] = polyparamsScaled[:,ii]*(1/yts)**(ii+1)
                # 2D array is already in correct format #
                polyparams2dScaled = np.copy(polyparamsScaled)
        else:
            polyparamsScaled   = polyparamsScaled*(1/yts)
            # 2D array is already in correct format #
            polyparams2dScaled = np.copy(polyparamsScaled)
        if(nper==1):
            dbreaks = np.zeros((nbas,1))
        else:
            dbreaks = np.copy(md.basalforcings.datebreaks)

        WriteData(fid, prefix, 'name', 'md.basalforcings.model', 'data', 9, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'groundedice_melting_rate', 'name', 'md.basalforcings.groundedice_melting_rate', 'format', 'DoubleMat', 'mattype', 1, 'scale', 1. / yts, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'geothermalflux', 'name', 'md.basalforcings.geothermalflux', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'num_basins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'num_params', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'num_breaks', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ar_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'ma_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'arma_timestep', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'basin_id', 'data', self.basin_id - 1, 'name', 'md.basalforcings.basin_id', 'format', 'IntMat', 'mattype', 2)  # 0-indexed
        WriteData(fid, prefix, 'data', polyparams2dScaled, 'name', 'md.basalforcings.polynomialparams', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'arlag_coefs', 'format', 'DoubleMat', 'name', 'md.basalforcings.arlag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'malag_coefs', 'format', 'DoubleMat', 'name', 'md.basalforcings.malag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'data', dbreaks, 'name', 'md.basalforcings.datebreaks', 'format', 'DoubleMat','scale',yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'deepwater_elevation', 'format', 'DoubleMat', 'name', 'md.basalforcings.deepwater_elevation')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_melting_rate', 'format', 'DoubleMat', 'name', 'md.basalforcings.upperwater_melting_rate', 'scale', 1 / yts, 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'upperwater_elevation', 'format', 'DoubleMat', 'name', 'md.basalforcings.upperwater_elevation')

    # }}}
