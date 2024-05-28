from math import log

import numpy as np

from structtoobj import structtoobj
from checkfield import checkfield
from project3d import project3d
from WriteData import WriteData

class basalforcingspico(object):
    """PICO BASAL FORCINGS class definition

       Usage:
          basalforcingspico=basalforcingspico();
    """
    def __init__(self, *args): # {{{
        self.num_basins                = 0
        self.basin_id                  = np.nan
        self.maxboxcount               = 0
        self.overturning_coeff         = np.nan
        self.gamma_T                   = 0.
        self.farocean_temperature      = np.nan
        self.farocean_salinity         = np.nan
        self.isplume                   = 0
        self.geothermalflux            = np.nan
        self.groundedice_melting_rate  = np.nan

        if len(args) == 0:
            self = self.setdefaultparameters()
        elif len(args) == 1:
            self.setdefaultparameters()
            self=structtoobj(self,args[0]);
        else:
            raise Exception('constructor not supported.')
    # }}}

    def __repr__(self): # {{{
        s = '   PICO basal melt rate parameterization:\n';
        s += '{}\n'.format(fielddisplay(self,'num_basins','number of basins the model domain is partitioned into [unitless]'))
        s += '{}\n'.format(fielddisplay(self,'basin_id','basin number assigned to each node [unitless]'))
        s += '{}\n'.format(fielddisplay(self,'maxboxcount','maximum number of boxes initialized under all ice shelves'))
        s += '{}\n'.format(fielddisplay(self,'overturning_coeff','overturning strength [m^3/s]'))
        s += '{}\n'.format(fielddisplay(self,'gamma_T','turbulent temperature exchange velocity [m/s]'))
        s += '{}\n'.format(fielddisplay(self,'farocean_temperature','depth averaged ocean temperature in front of the ice shelf for basin i [K]'))
        s += '{}\n'.format(fielddisplay(self,'farocean_salinity','depth averaged ocean salinity in front of the ice shelf for basin i [psu]'))
        s += '{}\n'.format(fielddisplay(self,'isplume','boolean to use buoyant plume melt rate parameterization from Lazeroms et al., 2018 (default false)'))
        s += '{}\n'.format(fielddisplay(self,'geothermalflux','geothermal heat flux [W/m^2]'))
        s += '{}\n'.format(fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) [m/yr]'))

        return s
    # }}}

    def extrude(self,md): # {{{
        self.basin_id=project3d(md,'vector',self.basin_id,'type','element','layer',1)
        self.geothermalflux=project3d(md,'vector',self.geothermalflux,'type','element','layer',1) #bedrock only gets geothermal flux
        self.groundedice_melting_rate=project3d(md,'vector',self.groundedice_melting_rate,'type','node','layer',1)

        return self
    # }}}

    def initialize(self,md): # {{{
        if np.isnan(self.maxboxcount):
            self.maxboxcount = 5
            print('      no maximum number of boxes set, setting value to 5')
        if np.isnan(self.overturning_coeff):
            self.overturning_coeff = 1e6*np.ones((md.mesh.numberofvertices,1)) #m^3/s
            print('      no overturning strength set, setting value to 1e6')
        if np.isnan(self.gamma_T):
            self.gamma_T = 2e-5 #m/s
            print('      no turbulent temperature exchange velocity set, setting value to 2e-5')
        if np.isnan(self.groundedice_melting_rate):
            self.groundedice_melting_rate=np.zeros((md.mesh.numberofvertices,1))
            print('      no basalforcings.groundedice_melting_rate specified: values set as zero')

        return self
    # }}}

    def setdefaultparameters(self): # {{{

        self.maxboxcount       = 5
        self.overturning_coeff = 1e6 #m^3/s
        self.gamma_T           = 2e-5 #m/s
        self.isplume           = False

    # }}}

    def checkconsistency(self,md,solution,analyses): # {{{

        md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0);
        md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements, 1])
        md = checkfield(md,'fieldname','basalforcings.maxboxcount','numel',1,'NaN',1,'Inf',1,'>',0)
        if np.size(self.overturning_coeff)==1:
            md = checkfield(md,'fieldname','basalforcings.overturning_coeff','numel',1,'NaN',1,'Inf',1,'>',0)
        else:
            md = checkfield(md,'fieldname','basalforcings.overturning_coeff','size',[md.mesh.numberofvertices, 1],'NaN',1,'Inf',1,'>',0)
        md = checkfield(md,'fieldname','basalforcings.gamma_T','numel',1,'NaN',1,'Inf',1,'>',0)
        md = checkfield(md,'fieldname','basalforcings.farocean_temperature','NaN',1,'Inf',1,'size',[md.basalforcings.num_basins+1, np.nan])
        md = checkfield(md,'fieldname','basalforcings.farocean_salinity','NaN',1,'Inf',1,'>',0,'size',[md.basalforcings.num_basins+1, np.nan])
        md = checkfield(md,'fieldname','basalforcings.isplume','values',[0, 1])
        md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1)
        md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)

        return md
    # }}}

    def marshall(self,prefix,md,fid): # {{{

        yts=md.constants.yts

        WriteData(fid,prefix,'name','md.basalforcings.model','data',5,'format','Integer') 
        WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer') 
        WriteData(fid,prefix,'object',self,'fieldname','maxboxcount','format','Integer') 
        WriteData(fid,prefix,'object',self,'fieldname','overturning_coeff','format','DoubleMat','mattype',1) 
        WriteData(fid,prefix,'object',self,'fieldname','gamma_T','format','Double') 
        WriteData(fid,prefix,'object',self,'fieldname','farocean_temperature','format','DoubleMat','name','md.basalforcings.farocean_temperature','timeserieslength',md.basalforcings.num_basins+1,'yts',md.constants.yts) 
        WriteData(fid,prefix,'object',self,'fieldname','farocean_salinity','format','DoubleMat','name','md.basalforcings.farocean_salinity','timeserieslength',md.basalforcings.num_basins+1,'yts',md.constants.yts) 
        WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2)    #Change to 0-indexing
        WriteData(fid,prefix,'object',self,'fieldname','isplume','format','Boolean') 
        WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts) 
        WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts) 

    # }}}
