#!/usr/bin/env python3
import numpy as np

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData
from xy2ll import xy2ll

class basalforcingsismip7(object):
    """ISMIP7 BASAL FORCINGS class definition

    Usage:
        basalforcings = basalforcingsismip7()
    """

    def __init__(self,*args):  # {{{
        self.num_basins                = 0
        self.basin_id                  = 0
        self.gamma                     = 0
        self.coriolis_f                = np.nan

        self.salinity                  = np.nan
        self.tf                        = np.nan
        self.tf_depths                 = np.nan

        self.geothermalflux = np.nan
        self.groundedice_melting_rate = np.nan

        if len(args) == 0:
            self.setdefaultparameters()
        elif len(args) == 1:
            self.setdefaultparameters()

            constructor = args[0]
            for field in vars(self).keys():
                if field in constructor.__dict__.keys():
                    setattr(self,field,getattr(constructor,field))
        else:
            raise Exception('constructuor not supported')
    # }}}
    def __repr__(self):  # {{{
        s = '   ISMIP7 basal melt rate parameterization:\n'
        s += '{}\n'.format(fielddisplay(self,'num_basins','[TODO] number of basins the model domain is partitioned into [unitless]'))
        s += '{}\n'.format(fielddisplay(self,'basin_id','[TODO] basin number assigned to each node (unitless)'))
        s += '{}\n'.format(fielddisplay(self,'gamma','melt rate coefficient (m/yr)'))
        s += '{}\n'.format(fielddisplay(self,'tf_depths','elevation of vertical layers in ocean thermal forcing dataset'))
        s += '{}\n'.format(fielddisplay(self,'tf','thermal forcing (ocean temperature minus freezing point) (degrees C)'))
        s += '{}\n'.format(fielddisplay(self,'salinity','salinity (psu)'))
        s += '{}\n'.format(fielddisplay(self,'coriolis_f','Coriolis parameter (s^-1)'))
        s += '{}\n'.format(fielddisplay(self,'geothermalflux','geothermal heat flux (W/m^2)'))
        s += '{}\n'.format(fielddisplay(self,'groundedice_melting_rate','basal melting rate (positive if melting) (m/yr)'))
        return s
    # }}}
    def extrude(self, md):  # {{{

        self.tf = project3d(md,'vector',self.tf,'type','node')
        self.salinity = project3d(md,'vector',self.salinity,'type','node')

        self.geothermalflux = project3d(md, 'vector', self.geothermalflux, 'type', 'node', 'layer', 1) # Bedrock only gets geothermal flux
        self.groundedice_melting_rate = project3d(md, 'vector', self.groundedice_melting_rate, 'type', 'node', 'layer', 1)
        return self
    # }}}
    def initialize(self, md):  # {{{
        # Update fixed-coriolis parameter
        if ~np.any(md.mesh.lat):
            print('      no md.mesh.lat specified.')
            if md.mesh.epsg == 3031: # For Antarctica
                lat, lon = xy2ll(md.mesh.x,md.mesh.y,-1)
            elif md.mesh.epsg == 3413: # For Greenland
                lat, lon = xy2ll(md.mesh.x,md.mesh.y,1)
            else:
                raise Exception('      md.mesh.lat not specified and cannot be calculated from md.mesh.epsg')
        else:
            lat = md.mesh.lat[:]

        omega=7.2921e-5 #angular velocity of the Earth (rad/s)
        self.coriolis_f=2*omega*np.sin(lat/180*np.pi)

        if self.gamma == 0:
            self.gamma = 14477
            print('      no basalforcings.gamma specified: value set to 14477 m/yr')
        if np.all(np.isnan(self.groundedice_melting_rate)):
            self.groundedice_melting_rate = np.zeros((md.mesh.numberofvertices, ))
            print('      no basalforcings.groundedice_melting_rate specified: values set as zero')
        return self
    # }}}
    def setdefaultparameters(self):  # {{{
        self.gamma = 0.0
        return self
    # }}}
    def checkconsistency(self, md, solution, analyses):  # {{{

        md = checkfield(md,'fieldname','basalforcings.num_basins','numel',1,'NaN',1,'Inf',1,'>',0)
        md = checkfield(md,'fieldname','basalforcings.basin_id','Inf',1,'>=',0,'<=',md.basalforcings.num_basins,'size',[md.mesh.numberofelements, 1])

        md = checkfield(md,'fieldname','basalforcings.gamma','numel',1,'NaN',1,'Inf',1,'>',0)

        md = checkfield(md,'fieldname','basalforcings.coriolis_f','size',[md.mesh.numberofvertices, 1],'NaN',1,'Inf',1)
        md = checkfield(md,'fieldname','basalforcings.tf_depths','NaN',1,'Inf',1,'size',[1,np.nan],'<=',0)
        md = checkfield(md,'fieldname','basalforcings.geothermalflux','NaN',1,'Inf',1,'>=',0,'timeseries',1)
        md = checkfield(md,'fieldname','basalforcings.groundedice_melting_rate','NaN',1,'Inf',1,'timeseries',1)

        ndepths = np.shape(self.tf_depths)[1]
        md = checkfield(md,'fieldname','basalforcings.tf','size',[ndepths])
        md = checkfield(md,'fieldname','basalforcings.salinity','size',[ndepths])
        for i in range(len(self.tf_depths)):
            md = checkfield(md,'fieldname','basalforcings.tf[' + str(i) + ']','field',md.basalforcings.tf[i],'size',[md.mesh.numberofvertices+1, np.nan],'NaN',1,'Inf',1,'>=',0,'timeseries',1)
            md = checkfield(md,'fieldname','basalforcings.salinity[' + str(i) + ']','field',md.basalforcings.salinity[i],'size',[md.mesh.numberofvertices+1, np.nan],'NaN',1,'Inf',1,'>=',0,'timeseries',1)

        return md
    # }}}
    def marshall(self, prefix, md, fid):  # {{{

        yts = md.constants.yts

        WriteData(fid,prefix,'name','md.basalforcings.model','data',10,'format','Integer')
        WriteData(fid,prefix,'object',self,'fieldname','num_basins','format','Integer')
        WriteData(fid,prefix,'object',self,'fieldname','basin_id','data',self.basin_id-1,'name','md.basalforcings.basin_id','format','IntMat','mattype',2)   #0-indexed
        WriteData(fid,prefix,'object',self,'fieldname','gamma','format','Double','scale',1/yts)
        WriteData(fid,prefix,'object',self,'fieldname','coriolis_f','format','DoubleMat','name','md.basalforcings.coriolis_f','mattype',1)
        WriteData(fid,prefix,'object',self,'fieldname','tf_depths','format','DoubleMat','name','md.basalforcings.tf_depths')
        WriteData(fid,prefix,'object',self,'fieldname','tf','format','MatArray','name','md.basalforcings.tf','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'fieldname','salinity','format','MatArray','name','md.basalforcings.salinity','timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'fieldname','geothermalflux','format','DoubleMat','name','md.basalforcings.geothermalflux','mattype',1,'timeserieslength',md.mesh.numberofelements+1,'yts',md.constants.yts)
        WriteData(fid,prefix,'object',self,'fieldname','groundedice_melting_rate','format','DoubleMat','mattype',1,'scale',1./yts,'timeserieslength',md.mesh.numberofvertices+1,'yts',md.constants.yts)

    # }}}
