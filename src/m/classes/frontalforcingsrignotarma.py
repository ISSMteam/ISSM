# -*- coding: utf-8 -*-
import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from MatlabFuncs import *
from WriteData import WriteData


class frontalforcingsrignotarma(object):
    """FRONTALFORCINGSRIGNOTARMA class definition

    Usage:
        frontalforcingsrignotarma = frontalforcingsrignotarma()
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
        self.monthlyvals_intercepts = np.nan
        self.monthlyvals_trends = np.nan
        self.monthlyvals_numbreaks = 0
        self.monthlyvals_datebreaks = np.nan
        self.basin_id = np.nan
        self.subglacial_discharge = np.nan
        self.isdischargearma = 0
        self.sd_ar_order = 0
        self.sd_ma_order = 0
        self.sd_arma_timestep = 0
        self.sd_arlag_coefs = np.nan
        self.sd_malag_coefs = np.nan
        self.sd_monthlyfrac = np.nan
        self.sd_num_breaks  = 0
        self.sd_num_params  = 0
        self.sd_polynomialparams = np.nan
        self.sd_datebreaks = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            error('constructor not supported')

    def __repr__(self):  # {{{
        s = '   Frontalforcings parameters:\n'
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
        s += '{}\n'.format(fielddisplay(self, 'isdischargearma','whether an ARMA model is also used for the subglacial discharge (if 0: subglacial_discharge is used, if 1: sd_ parameters are used)'))
        s += '{}\n'.format(fielddisplay(self, 'subglacial_discharge', 'sum of subglacial discharge for each basin [m/d]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_ar_order','order of the subglacial discharge autoregressive model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_ma_order','order of the subglacial discharge moving-average model [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_arma_timestep','time resolution of the subglacial discharge autoregressive model [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_arlag_coefs','basin-specific vectors of AR lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_malag_coefs','basin-specific vectors of MA lag coefficients for subglacial discharge [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_monthlyfrac','basin-specific vectors of 12 values with fraction of the annual discharge occuring every month [unitless]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_num_params','number of different parameters in the subglacial discharge piecewise-polynomial (1:intercept only, 2:with linear trend, 3:with quadratic trend, etc.)'))
        s += '{}\n'.format(fielddisplay(self, 'sd_num_breaks','number of different breakpoints in the subglacial discharge piecewise-polynomial (separating sd_num_breaks+1 periods)'))
        s += '{}\n'.format(fielddisplay(self, 'sd_datebreaks','dates at which the breakpoints in the piecewise polynomial occur (1 row per basin) [yr]'))
        s += '{}\n'.format(fielddisplay(self, 'sd_polynomialparams','coefficients for the sd_polynomial (const,trend,quadratic,etc.),dim1 for basins,dim2 for periods,dim3 for orders'))
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        self.basin_id = np.nan
        self.num_basins = 0
        self.subglacial_discharge = np.nan
        self.ar_order = 0.0  # Autoregression model of order 0
        self.ma_order = 0.0  # Moving-average model of order 0
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not (solution == 'TransientSolution') or not md.transient.ismovingfront:
            return md

        nbas  = md.frontalforcings.num_basins;
        nprm  = md.frontalforcings.num_params;
        nbrk  = md.frontalforcings.num_breaks;
        nMbrk = md.frontalforcings.monthlyvals_numbreaks;
        md = checkfield(md, 'fieldname', 'frontalforcings.num_basins', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.num_params', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.num_breaks', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.basin_id', 'Inf', 1, '>=', 0, '<=', md.frontalforcings.num_basins, 'size', [md.mesh.numberofelements])
        md = checkfield(md, 'fieldname', 'frontalforcings.subglacial_discharge', '>=', 0, 'NaN', 1, 'Inf', 1, 'timeseries', 1)
        if len(np.shape(self.polynomialparams)) == 1:
            self.polynomialparams = np.array([[self.polynomialparams]])
        if(nbas>1 and nbrk>=1 and nprm>1):
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk+1,nprm],'numel',nbas*(nbrk+1)*nprm)
        elif(nbas==1):
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(nbrk+1)*nprm)
        elif(nbrk==0):
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nprm],'numel',nbas*(nbrk+1)*nprm)
        elif(nprm==1):
            md = checkfield(md,'fieldname','frontalforcings.polynomialparams','NaN',1,'Inf',1,'size',[nbas,nbrk],'numel',nbas*(nbrk+1)*nprm)
        md = checkfield(md, 'fieldname', 'frontalforcings.ar_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.ma_order', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
        md = checkfield(md, 'fieldname', 'frontalforcings.arma_timestep', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', md.timestepping.time_step) # ARMA time step cannot be finer than ISSM timestep
        md = checkfield(md, 'fieldname', 'frontalforcings.arlag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.frontalforcings.num_basins, md.frontalforcings.ar_order])
        md = checkfield(md, 'fieldname', 'frontalforcings.malag_coefs', 'NaN', 1, 'Inf', 1, 'size', [md.frontalforcings.num_basins, md.frontalforcings.ma_order])
        if(nbrk>0):
            md = checkfield(md, 'fieldname', 'frontalforcings.datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,nbrk])
        elif(np.size(md.frontalforcings.datebreaks)==0 or np.all(np.isnan(md.frontalforcings.datebreaks))):
            pass
        else:
            raise RuntimeError('md.frontalforcings.num_breaks is 0 but md.frontalforcings.datebreaks is not empty')
        

        ### Check if some monthly forcings are provided ###
        if(np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts))==False or np.all(np.isnan(md.frontalforcings.monthlyvals_trends))==False or np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))==False):
            isMonthly = True
        else:
            isMonthly = False
        if(np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))):
            isMonthlyTrend = True
        else:
            isMonthlyTrend = False
        if(isMonthly):
            md = checkfield(md, 'fieldname', 'frontalforcings.monthlyvals_numbreaks', 'numel', 1, 'NaN', 1, 'Inf', 1, '>=', 0)
            if(nbas>1 and nMbrk>=1):
                md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nbas,12,nMbrk+1],'numel',nbas*(nMbrk+1)*12)
                if(isMonthlyTrend):
                    md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nbas,12,nMbrk+1],'numel',nbas*(nMbrk+1)*12)
            elif(nbas==1):
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nMbrk+1,12],'numel',nbas*(nMbrk+1)*12)
               if(isMonthlyTrend):
                  md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nMbrk+1,12],'numel',nbas*(nMbrk+1)*12)
            elif(nMbrk==0):
               md = checkfield(md,'fieldname','frontalforcings.monthlyvals_intercepts','NaN',1,'Inf',1,'size',[nbas,12],'numel',nbas*(nMbrk+1)*12)
               if(isMonthlyTrend):
                  md = checkfield(md,'fieldname','frontalforcings.monthlyvals_trends','NaN',1,'Inf',1,'size',[nbas,12],'numel',nbas*(nMbrk+1)*12)
        if(nMbrk>0):
            md = checkfield(md, 'fieldname', 'frontalforcings.monthlyvals_datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,nMbrk])
        elif(np.size(md.frontalforcings.monthlyvals_datebreaks)==0 or np.all(np.isnan(md.frontalforcings.monthlyvals_datebreaks))):
            pass
        else:
            raise RuntimeError('md.frontalforcings.monthlyvals_numbreaks is 0 but md.frontalforcings.monthlyvals_datebreaks is not empty')

        ### Chacking subglacial discharge ###
        md = checkfield(md, 'fieldname', 'frontalforcings.isdischargearma', 'values', [0, 1])
        if(self.isdischargearma==0):
            md = checkfield(md,'fieldname','frontalforcings.subglacial_discharge','>=',0,'NaN',1,'Inf',1,'timeseries',1)
        else:
            sdnbrk  = md.frontalforcings.sd_num_breaks
            sdnprm  = md.frontalforcings.sd_num_params
            md = checkfield(md,'fieldname','frontalforcings.sd_ar_order','numel',1,'NaN',1,'Inf',1,'>=',0)
            md = checkfield(md,'fieldname','frontalforcings.sd_ma_order','numel',1,'NaN',1,'Inf',1,'>=',0)
            md = checkfield(md,'fieldname','frontalforcings.sd_arma_timestep','numel',1,'NaN',1,'Inf',1,'>=',max(1,md.timestepping.time_step)) #ARMA time step cannot be finer than ISSM timestep and annual timestep
            md = checkfield(md,'fieldname','frontalforcings.sd_arlag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.sd_ar_order])
            md = checkfield(md,'fieldname','frontalforcings.sd_malag_coefs','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,md.frontalforcings.sd_ma_order])
            md = checkfield(md,'fieldname','frontalforcings.sd_monthlyfrac','NaN',1,'Inf',1,'size',[md.frontalforcings.num_basins,12])
            if(np.any(abs(np.sum(self.sd_monthlyfrac,axis=1)-1)>1e-3)):
                raise RuntimeError('the 12 entries for each basin of md.frontalforcings.sd_monthlyfrac should add up to 1')
            md = checkfield(md,'fieldname','frontalforcings.sd_num_params','numel',1,'NaN',1,'Inf',1,'>',0)
            md = checkfield(md,'fieldname','frontalforcings.sd_num_breaks','numel',1,'NaN',1,'Inf',1,'>=',0)
            if len(np.shape(self.sd_polynomialparams)) == 1:
                self.sd_polynomialparams = np.array([[self.sd_polynomialparams]])
            if(nbas>1 and sdnbrk>=1 and sdnprm>1):
                md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnbrk+1,sdnprm],'numel',nbas*(sdnbrk+1)*sdnprm)
            elif(nbas==1):
                md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nprm,nbrk+1],'numel',nbas*(sdnbrk+1)*sdnprm)
            elif(sdnbrk==0):
                md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnprm],'numel',nbas*(sdnbrk+1)*sdnprm)
            elif(sdnprm==1):
                md = checkfield(md,'fieldname','frontalforcings.sd_polynomialparams','NaN',1,'Inf',1,'size',[nbas,sdnbrk],'numel',nbas*(sdnbrk+1)*sdnprm)
            if(sdnbrk>0):
                md = checkfield(md, 'fieldname', 'frontalforcings.sd_datebreaks', 'NaN', 1, 'Inf', 1, 'size', [nbas,sdnbrk])
            elif(np.size(md.frontalforcings.sd_datebreaks)==0 or np.all(np.isnan(md.frontalforcings.sd_datebreaks))):
                pass
            else:
                raise RuntimeError('md.frontalforcings.sd_num_breaks is 0 but md.frontalforcings.sd_datebreaks is not empty')

        return md
    # }}}

    def extrude(self, md):  # {{{
        # Nothing for now
        return self
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts
        nbas = md.frontalforcings.num_basins;
        nprm = md.frontalforcings.num_params;
        nper = md.frontalforcings.num_breaks+1;
        # Scale the parameters #
        polyparamsScaled   = np.copy(md.frontalforcings.polynomialparams)
        polyparams2dScaled = np.zeros((nbas,nper*nprm))
        if(nprm>1):
            # Case 3D #
            if(nbas>1 and nper>1):
                for ii in range(nprm):
                    polyparamsScaled[:,:,ii] = polyparamsScaled[:,:,ii]*(1/yts)**ii
                # Fit in 2D array #
                for ii in range(nprm):
                    polyparams2dScaled[:,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[:,:,ii]
            # Case 2D and higher-order params at increasing row index #
            elif(nbas==1):
                for ii in range(nprm):
                    polyparamsScaled[ii,:] = polyparamsScaled[ii,:]*(1/yts)**ii
                # Fit in row array #
                for ii in range(nprm):
                    polyparams2dScaled[0,ii*nper:(ii+1)*nper] = 1*polyparamsScaled[ii,:]
            # Case 2D and higher-order params at incrasing column index #
            elif(nper==1):
                for ii in range(nprm):
                    polyparamsScaled[:,ii] = polyparamsScaled[:,ii]*(1/yts)**ii
                # 2D array is already in correct format #
                polyparams2dScaled = np.copy(polyparamsScaled)
        else:
            # 2D array is already in correct format and no need for scaling #
            polyparams2dScaled = np.copy(polyparamsScaled)
        if(nper==1):
            dbreaks = np.zeros((nbas,1))
        else:
            dbreaks = np.copy(md.frontalforcings.datebreaks)

        ### Deal with montly effects ###
        nMper = md.frontalforcings.monthlyvals_numbreaks+1
        if(np.any(np.isnan(md.frontalforcings.monthlyvals_intercepts))):
            interceptsM = np.zeros((nbas,12)) #monthly intercepts not provided, set to 0
            trendsM     = np.zeros((nbas,12)) #set monthly trends also to 0
        else:
            interceptsM3d = md.frontalforcings.monthlyvals_intercepts
            if(np.any(np.isnan(md.frontalforcings.monthlyvals_trends))):
                trendsM3d = 0*interceptsM3d #monthly trends not provided, set to 0
            else:
                trendsM3d = md.frontalforcings.monthlyvals_trends
        # Create 2D arrays from 3D arrays if needed #
        if(nMper>1 and np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts))==False):
            interceptsM = np.zeros((nbas,12*nMper)) 
            trendsM     = np.zeros((nbas,12*nMper))
            for ii in range(nMper):
                interceptsM[:,ii*12:(ii+1)*12] = 1*interceptsM3d[:,:,ii]
                trendsM[:,ii*12:(ii+1)*12] = 1*trendsM3d[:,:,ii]
        elif(nMper==1 and np.all(np.isnan(md.frontalforcings.monthlyvals_intercepts))==False):
            interceptsM = 1*interceptsM3d
            trendsM     = 1*trendsM3d
        if(nMper==1):
            dMbreaks = np.zeros((nbas,1))
        else:
            dMbreaks = np.copy(md.frontalforcings.monthlyvals_datebreaks)

        ### Deal with the subglacial discharge polynomial ###
        if(self.isdischargearma):
            sdnprm  = md.frontalforcings.sd_num_params
            sdnper  = md.frontalforcings.sd_num_breaks+1
            sdpolyparamsScaled   = np.copy(md.frontalforcings.sd_polynomialparams)
            sdpolyparams2dScaled = np.zeros((nbas,sdnper*sdnprm))
            if(sdnprm>1):
                # Case 3D #
                if(nbas>1 and sdnper>1):
                    for ii in range(sdnprm):
                        sdpolyparamsScaled[:,:,ii] = sdpolyparamsScaled[:,:,ii]*(1/yts)**ii
                    # Fit in 2D array #
                    for ii in range(sdnprm):
                        sdpolyparams2dScaled[:,ii*sdnper:(ii+1)*sdnper] = 1*sdpolyparamsScaled[:,:,ii]
                # Case 2D and higher-order params at increasing row index #
                elif(nbas==1):
                    for ii in range(sdnprm):
                        sdpolyparamsScaled[ii,:] = sdpolyparamsScaled[ii,:]*(1/yts)**ii
                    # Fit in row array #
                    for ii in range(nprm):
                        sdpolyparams2dScaled[0,ii*sdnper:(ii+1)*sdnper] = 1*sdpolyparamsScaled[ii,:]
                # Case 2D and higher-order params at incrasing column index #
                elif(sdnper==1):
                    for ii in range(sdnprm):
                        sdpolyparamsScaled[:,ii] = sdpolyparamsScaled[:,ii]*(1/yts)**ii
                    # 2D array is already in correct format #
                    sdpolyparams2dScaled = np.copy(sdpolyparamsScaled)
            else:
                # 2D array is already in correct format and no need for scaling #
                sdpolyparams2dScaled = np.copy(sdpolyparamsScaled)
            if(sdnper==1):
                sd_dbreaks = np.zeros((nbas,1))
            else:
                sd_dbreaks = np.copy(md.frontalforcings.sd_datebreaks)




        WriteData(fid, prefix, 'name', 'md.frontalforcings.parameterization', 'data', 3, 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'num_basins', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'num_breaks', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'num_params', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'ar_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'ma_order', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'arma_timestep', 'format', 'Double', 'scale', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'basin_id', 'data', self.basin_id - 1, 'name', 'md.frontalforcings.basin_id', 'format', 'IntMat', 'mattype', 2)  # 0-indexed
        WriteData(fid, prefix, 'data', polyparams2dScaled, 'name', 'md.frontalforcings.polynomialparams', 'format', 'DoubleMat')
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'arlag_coefs', 'format', 'DoubleMat', 'name', 'md.frontalforcings.arlag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'malag_coefs', 'format', 'DoubleMat', 'name', 'md.frontalforcings.malag_coefs', 'yts', yts)
        WriteData(fid, prefix, 'data', dbreaks, 'name', 'md.frontalforcings.datebreaks', 'format', 'DoubleMat','scale',yts)
        WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','monthlyvals_numbreaks','format','Integer')
        WriteData(fid,prefix,'data',dMbreaks,'name','md.frontalforcings.monthlyvals_datebreaks','format','DoubleMat','scale',yts)
        WriteData(fid,prefix,'data',interceptsM,'name','md.frontalforcings.monthlyvals_intercepts','format','DoubleMat')
        WriteData(fid,prefix,'data',trendsM,'name','md.frontalforcings.monthlyvals_trends','format','DoubleMat','scale',1/yts)
        WriteData(fid,prefix,'object',self,'fieldname','isdischargearma','format','Boolean')
        if(self.isdischargearma==0):
            WriteData(fid, prefix, 'object', self, 'class', 'frontalforcings', 'fieldname', 'subglacial_discharge', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', yts)
        else:
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_num_breaks','format','Integer')
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_num_params','format','Integer')
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_ar_order','format','Integer')
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_ma_order','format','Integer')
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_arma_timestep','format','Double','scale',yts)
            WriteData(fid,prefix,'data',sdpolyparams2dScaled,'name','md.frontalforcings.sd_polynomialparams','format','DoubleMat')
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_arlag_coefs','format','DoubleMat','name','md.frontalforcings.sd_arlag_coefs','yts',yts)
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_malag_coefs','format','DoubleMat','name','md.frontalforcings.sd_malag_coefs','yts',yts)
            WriteData(fid,prefix,'data',sd_dbreaks,'name','md.frontalforcings.sd_datebreaks','format','DoubleMat','scale',yts)
            WriteData(fid,prefix,'object',self,'class','frontalforcings','fieldname','sd_monthlyfrac','format','DoubleMat','name','md.frontalforcings.sd_monthlyfrac','yts',yts)
    # }}}


