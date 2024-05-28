import numpy as np

from checkfield import *
from fielddisplay import fielddisplay
from project3d import *
from WriteData import *
from GetAreas import *

class hydrologyarmapw(object):
    """HYDROLOGYARMAPW class definition

    Usage:
        hydrologyarmapw = hydrologyarmapw()
    """

    def __init__(self, *args):  # {{{
        self.num_basins = 0
        self.num_params = 0
        self.num_breaks = 0
        self.polynomialparams = np.nan
        self.arma_timestep = 0
        self.ar_order = 0
        self.ma_order = 0
        self.arlag_coefs = np.nan
        self.malag_coefs = np.nan
        self.datebreaks = np.nan
        self.basin_id = np.nan
        self.monthlyfactors = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   hydrologyarmapw\n'
        s += 'subglacial water pressure is calculated as Pw=monthlyfactor[month]*(rho_water*g*bed+Pw_arma) where Pw_arma is the perturbation calculated as an ARMA process\n'
        s += 'polynomialparams includes the constant, linear trend, quadratic trend, etc. of the ARMA process\n'
        s += 'arlag_coefs and malag_coefs include the coefficients of the ARMA process\n'
        s += '{}\n'.format(fielddisplay(self, 'num_basins', 'number of different basins [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'basin_id', 'basin number assigned to each element [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'num_breaks', 'number of different breakpoints in the piecewise-polynomial (separating num_breaks+1 periods)'))
        s += '{}\n'.format(fielddisplay(self, 'num_params', 'number of different parameters in the piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(fielddisplay(self, 'monthlyfactors', 'monthly multiplicative factor on the subglacial water pressure, specified per basin (size:[num_basins,12])'))
        s += '{}\n'.format(fielddisplay(self, 'polynomialparams', 'coefficients for the polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders, ex: polyparams=cat(num_params,intercepts,trendlinearcoefs,trendquadraticcoefs)'))
        s += '{}\n'.format(fielddisplay(self, 'datebreaks', 'dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'ar_order', 'order of the autoregressive model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'ma_order', 'order of the moving-average model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'arma_timestep', 'time resolution of the ARMA model [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'arlag_coefs', 'basin-specific vectors of AR lag coefficients [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'malag_coefs', 'basin-specific vectors of MA lag coefficients [unitless]'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # No default parameters
        return self # Nothing for now
    # }}}

    def extrude(self, md):  # {{{
        self.basin_id = project3d(md,'vector',self.basin_id,'type','element')
        return self # Nothing for now
    # }}}

    def defaultoutputs(self, md):  # {{{
        return ['FrictionWaterPressure']
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        #Early return
        if 'HydrologyArmapwAnalysis' not in analyses:
            return md

        nbas = md.hydrology.num_basins
        nprm = md.hydrology.num_params
        nbrk = md.hydrology.num_breaks
        
        md = checkfield(md, 'fieldname', 'hydrology.num_basins', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'hydrology.num_params', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'hydrology.num_breaks', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.basin_id', 'Inf', 1, '>=', 0, '<=', md.hydrology.num_basins, 'size', [md.mesh.numberofelements])

        # Check if monthlyfactors are provided
        if(np.size(md.hydrology.monthlyfactors)>1 or np.all(np.isnan(md.hydrology.monthlyfactors))==False):
            md = checkfield(md,'fieldname','hydrology.monthlyfactors','NaN',1,'Inf',1,'size',[md.hydrology.num_basins,12])
            if(np.any(md.hydrology.monthlyfactors!=1) and md.timestepping.time_step>=1):
                raise RuntimeError('md.timestepping.time_step is too large to use hydrologyarmapw() with monthlyfactors')

        if len(np.shape(self.polynomialparams)) == 1:
            self.polynomialparams = np.array([[self.polynomialparams]])
        if(nbas>1 and nbrk>=1 and nprm>1):
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm) 
        elif(nbas==1):
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm) 
        elif(nbrk==0):
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm)
        elif(nprm==1):
            md = checkfield(md,'fieldname','hydrology.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1],'numel',nbas*(nbrk+1)*nprm)
        md = checkfield(md, 'fieldname', 'hydrology.ar_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.ma_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'hydrology.arma_timestep', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', md.timestepping.time_step) # Autoregression time step cannot be finer than ISSM timestep
        md = checkfield(md, 'fieldname', 'hydrology.arlag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.hydrology.num_basins, md.hydrology.ar_order])
        md = checkfield(md, 'fieldname', 'hydrology.malag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.hydrology.num_basins, md.hydrology.ma_order])
        if(nbrk>0):
            md = checkfield(md, 'fieldname', 'hydrology.datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,nbrk])
        elif(np.size(md.hydrology.datebreaks)==0 or np.all(np.isnan(md.hydrology.datebreaks))):
            pass
        else:
            raise RuntimeError('md.hydrology.num_breaks is 0 but md.hydrology.datebreaks is not empty')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        nbas = md.hydrology.num_basins;
        nprm = md.hydrology.num_params;
        nper = md.hydrology.num_breaks+1;
        # Scale the parameters #
        polyparamsScaled   = np.copy(md.hydrology.polynomialparams)
        polyparams2dScaled = np.zeros((nbas,nper*nprm))
        if(nprm>1):
            # Case 3D #
            if(nbas>1 and nper>1):
                for ii in range(nprm):
                    polyparamsScaled[:,:,ii] = polyparamsScaled[:,:,ii]*(1/yts)**(ii)
                # Fit in 2D array #
                for ii in range(nprm):
                    polyparams2dScaled[:,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[:,:,ii]
            # Case 2D and higher-order params at increasing row index #
            elif(nbas==1):
                for ii in range(nprm):
                    polyparamsScaled[ii,:] = polyparamsScaled[ii,:]*(1/yts)**(ii)
                # Fit in row array #
                for ii in range(nprm):
                    polyparams2dScaled[0,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[ii,:]
            # Case 2D and higher-order params at incrasing column index #
            elif(nper==1):
                for ii in range(nprm):
                    polyparamsScaled[:,ii] = polyparamsScaled[:,ii]*(1/yts)**(ii)
                # 2D array is already in correct format #
                polyparams2dScaled = np.copy(polyparamsScaled)
        else:
            # 2D array is already in correct format and no need for scaling#
            polyparams2dScaled = np.copy(polyparamsScaled)
        
        if(nper==1):
            dbreaks = np.zeros((nbas,1))
        else:
            dbreaks = np.copy(md.hydrology.datebreaks)

        # If no monthlyfactors provided: set them all to 1 #
        if(np.size(md.hydrology.monthlyfactors)==1):
            tempmonthlyfactors = np.ones((nbas,12))
        else:
            tempmonthlyfactors = np.copy(md.hydrology.monthlyfactors)

        WriteData(fid, prefix, 'name', 'md.hydrology.model', 'data', 7, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'num_basins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'num_breaks', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'num_params', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'ar_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'ma_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'arma_timestep', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'basin_id', 'data', self.basin_id - 1, 'name', 'md.hydrology.basin_id', 'format', 'IntMat', 'mattype', 2)  # 0-indexed
        WriteData(fid, prefix, 'data', polyparams2dScaled, 'name', 'md.hydrology.polynomialparams', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'arlag_coefs', 'format', 'DoubleMat', 'name', 'md.hydrology.arlag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'hydrology', 'fieldname', 'malag_coefs', 'format', 'DoubleMat', 'name', 'md.hydrology.malag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'data', dbreaks, 'name', 'md.hydrology.datebreaks', 'format', 'DoubleMat','scale',yts)
        WriteData(fid,prefix,'data',tempmonthlyfactors,'name','md.hydrology.monthlyfactors','format','DoubleMat')
        WriteData(fid,prefix,'data',{'FrictionWaterPressure'},'name','md.hydrology.requested_outputs','format','StringArray')

    # }}}
