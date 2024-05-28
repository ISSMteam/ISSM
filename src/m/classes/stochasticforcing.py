import numpy as np
from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class stochasticforcing(object):
    """stochasticforcing class definition

    Usage:
        stochasticforcing = stochasticforcing()
    """

    def __init__(self, *args):  # {{{
        self.isstochasticforcing = 0
        self.fields = np.nan
        self.defaultdimension = 0
        self.default_id = np.nan
        self.covariance = np.nan
        self.timecovariance = np.nan
        self.stochastictimestep = 0
        self.randomflag = 1

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise RuntimeError('constructor not supported for stochasticforcing')

    def __repr__(self):  # {{{
        s = '   stochasticforcing parameters:\n'
        s += '{}\n'.format(fielddisplay(self, 'isstochasticforcing', 'is stochasticity activated?'))
        s += '{}\n'.format(fielddisplay(self, 'fields', 'fields with stochasticity applied, ex: [\'SMBautoregression\'], or [\'SMBforcing\',\'DefaultCalving\']'))
        s += '{}\n'.format(fielddisplay(self, 'defaultdimension', 'dimensionality of the noise terms (does not apply to fields with their specific dimension)'))
        s += '{}\n'.format(fielddisplay(self, 'default_id', 'id of each element for partitioning of the noise terms (does not apply to fields with their specific partition)'))
        s += '{}\n'.format(fielddisplay(self, 'covariance', 'covariance matrix for within- and between-fields covariance (units must be squared field units),multiple matrices can be concatenated along 3rd dimension to apply different covariances in time'))
        s += '{}\n'.format(fielddisplay(self, 'timecovariance', 'starting dates at which covariances apply (only applicabe if multiple covariance matrices are prescribed)'))
        s += '{}\n'.format(fielddisplay(self, 'stochastictimestep', 'timestep at which new stochastic noise terms are generated (default: md.timestepping.time_step)'))
        s += '{}\n'.format(fielddisplay(self, 'randomflag', 'whether to apply real randomness (true) or pseudo-randomness with fixed seed (false)'))
        s += 'Available fields:\n'
        s += '   BasalforcingsDeepwaterMeltingRatearma\n'
        s += '   BasalforcingsSpatialDeepwaterMeltingRate\n'
        s += '   DefaultCalving\n'
        s += '   FloatingMeltRate\n'
        s += '   FrictionWaterPressure\n'
        s += '   FrictionCoulombWaterPressure\n'
        s += '   FrictionSchoofWaterPressure\n'
        s += '   FrontalForcingsRignotarma (thermal forcing)\n'
        s += '   FrontalForcingsSubglacialDischargearma\n'
        s += '   SMBarma\n'
        s += '   SMBforcing\n'
        return s
    # }}}

    def setdefaultparameters(self):  # {{{
        # Type of stabilization used
        self.isstochasticforcing = 0  # stochasticforcing is turned off by default
        self.fields = []  # Need to initialize to list to avoid "RuntimeError: object of type 'float' has no len()" on import of class
        self.randomflag = 1  # true randomness is implemented by default
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        # Early return
        if not self.isstochasticforcing:
            return md

        num_fields = len(self.fields)
        if self.stochastictimestep == 0:
            md.stochasticforcing.stochastictimestep = md.timestepping.time_step #by default: stochastictimestep set to ISSM time step
            print('      stochasticforcing.stocahstictimestep not specified: set to md.timestepping.time_step')

        if(len(np.shape(self.covariance))==3):
           numtcovmat = np.shape(self.covariance)[2] #number of covariance matrices in time
           lsCovmats = []
           for ii in range(numtcovmat):
               lsCovmats.append(self.covariance[:,:,ii])
               try:
                   np.linalg.cholesky(self.covariance[:,:,ii])
               except:
                   raise TypeError('an entry in md.stochasticforcing.covariance is not positive definite')
        elif(len(np.shape(self.covariance))==2):
            numtcovmat = 1
            lsCovmats = [self.covariance]
            # Check that covariance matrix is positive definite (this is done internally by linalg)
            try:
                np.linalg.cholesky(self.covariance)
            except:
                raise TypeError('md.stochasticforcing.covariance is not positive definite')

        # Check that all fields agree with the corresponding md class and if any field needs the default params
        checkdefaults = False  # Need to check defaults only if one of the fields does not have its own dimensionality
        structstoch = self.structstochforcing()
        # Check if hydrologyarmapw is used
        if((type(md.hydrology).__name__ == 'hydrologyarmapw') and md.transient.ishydrology==1):
            ispwHydroarma = 1
        else:
            ispwHydroarma = 0
        for field in self.fields:
            # Checking agreement of classes
            if 'SMBarma' in field:
                mdname = structstoch[field]
                if type(md.smb).__name__ != mdname:
                    raise TypeError('md.smb does not agree with stochasticforcing field {}'.format(field))
            if 'SMBforcing' in field:
                mdname = structstoch[field]
                if type(md.smb).__name__ != mdname:
                    raise TypeError('md.smb does not agree with stochasticforcing field {}'.format(field))
            if 'FrontalForcings' in field:
                mdname = structstoch[field]
                if type(md.frontalforcings).__name__ != mdname:
                    raise TypeError('md.frontalforcings does not agree with stochasticforcing field {}'.format(field))
            if 'Calving' in field:
                mdname = structstoch[field]
                if type(md.calving).__name__ != mdname:
                    raise TypeError('md.calving does not agree with stochasticforcing field {}'.format(field))
            if 'BasalforcingsFloatingice' in field:
                mdname = structstoch[field]
                if type(md.basalforcings).__name__ != mdname:
                    raise TypeError('md.basalforcings does not agree with stochasticforcing field {}'.format(field))
            if 'BasalforcingsSpatialDeepwaterMeltingRate' in field:
                mdname = structstoch[field]
                if type(md.basalforcings).__name__ != mdname:
                    raise TypeError('md.basalforcings does not agree with stochasticforcing field {}'.format(field))
            if 'BasalforcingsDeepwaterMeltingRatearma' in field:
                mdname = structstoch[field]
                if type(md.basalforcings).__name__ != mdname:
                    raise TypeError('md.basalforcings does not agree with stochasticforcing field {}'.format(field))
            if 'WaterPressure' in field:
                #mdname = structstoch[field]
                mdnames = ['friction','frictioncoulomb','frictionschoof']
                found   = 0
                for ii in range(len(mdnames)):
                    if type(md.friction).__name__ == mdnames[ii]:
                        found = 1
                if not found:
                    raise TypeError('md.friction does not agree with stochasticforcing field {}'.format(field))
                #if (type(md.friction).__name__ != mdname):
                #    raise TypeError('md.friction does not agree with stochasticforcing field {}'.format(field))
                if type(md.friction).__name__ == 'friction' or type(md.friction).__name__ == 'frictionschoof' or type(md.friction).__name__=='frictioncoulomb':
                    if md.friction.coupling not in[0, 1, 2]:
                        raise TypeError('stochasticforcing field {} is only implemented for cases md.friction.coupling 0 or 1 or 2'.format(field))
                if type(md.friction).__name__ == 'friction':
                    if (np.any(md.friction.q == 0)):
                        raise TypeError('stochasticforcing field {} requires non-zero q exponent'.format(field))

            # Checking for specific dimensions
            if field not in ['SMBarma', 'FrontalForcingsRignotarma','BasalforcingsDeepwaterMeltingRatearma']  and not (field == 'FrictionWaterPressure' and ispwHydroarma == True):
                checkdefaults = True  # field with non-specific dimensionality

        # Retrieve sum of all the field dimensionalities
        dimensions = self.defaultdimension * np.ones((num_fields))
        indSMBarma   = -1  # About to check for index of SMBarma
        indTFarma    = -1  # About to check for index of FrontalForcingsRignotarma
        indSdarma    = -1  # About to check for index of FrontalForcingsSubglacialDischargearma
        indBDWarma   = -1  # About to check for index of BasalforcingsDeepwaterMeltingRatearma
        indPwarma    = -1  # About to check for index of hydrologyarmapw
        if 'SMBarma' in self.fields:
            indSMBarma = self.fields.index('SMBarma')  # Index of SMBarma, now check for consistency with other timesteps
            dimensions[indSMBarma] = md.smb.num_basins
            if(md.smb.arma_timestep<self.stochastictimestep):
                raise TypeError('SMBarma cannot have a timestep shorter than stochastictimestep')
        if 'FrontalForcingsRignotarma' in self.fields:
            indTFarma = self.fields.index('FrontalForcingsRignotarma')  # Index of TFarma, now check for consistency with other timesteps
            dimensions[indTFarma] = md.frontalforcings.num_basins
            if md.frontalforcings.arma_timestep < self.stochastictimestep:
                raise TypeError('FrontalForcingsRignotarma cannot have a timestep shorter than stochastictimestep')
        if 'FrontalForcingsSubglacialDischargearma' in self.fields:
            indSdarma = self.fields.index('FrontalForcingsSubglacialDischargearma')  # Index of Sdarma, now check for consistency with other timesteps
            dimensions[indSdarma] = md.frontalforcings.num_basins
            if md.frontalforcings.sd_arma_timestep < self.stochastictimestep:
                raise TypeError('FrontalForcingsSubglacialDischargearma cannot have a timestep shorter than stochastictimestep')
        if 'BasalforcingsDeepwaterMeltingRatearma' in self.fields:
            indBDWarma = self.fields.index('BasalforcingsDeepwaterMeltingRatearma')  # Index of BDWarma, now check for consistency with other timesteps
            dimensions[indTFarma] = md.basalforcings.num_basins
            if md.basalforcings.arma_timestep < self.stochastictimestep:
                raise TypeError('BasalforcingsDeepwaterMeltingRatearma cannot have a timestep shorter than stochastictimestep')
        if 'FrictionWaterPressure' in self.fields and ispwHydroarma:
            indPwarma = self.fields.index('FrictionWaterPressure')  # Index of Pwarma, now check for consistency with other timesteps
            dimensions[indPwarma] = md.hydrology.num_basins
            if md.hydrology.arma_timestep < self.stochastictimestep:
                raise TypeError('hydrologyarmapw cannot have a timestep shorter than stochastictimestep')
        size_tot = np.sum(dimensions)

        if indBDWarma != -1:
            if indPwarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indBDWarma]).astype(int):np.sum(dimensions[0:indBDWarma + 1]).astype(int), np.sum(dimensions[0:indPwarma]).astype(int):np.sum(dimensions[0:indPwarma + 1]).astype(int)]
                    if md.smb.arma_timestep != md.hydrology.arma_timestep and np.any(covsum != 0):
                        raise IOError('BasalforcingsDeepwaterMeltingRatarma and hydrologyarmapw have different arma_timestep and non-zero covariance')
            elif indSdarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSdarma]).astype(int):np.sum(dimensions[0:indSdarma + 1]).astype(int), np.sum(dimensions[0:indBDWarma]).astype(int):np.sum(dimensions[0:indBDWarma + 1]).astype(int)]
                    if md.frontalforcings.sd_arma_timestep != md.basalforcings.arma_timestep and np.any(covsum != 0):
                        raise IOError('FrontalForcingsSubglacialDischargearma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance')
            elif indSMBarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSMBarma]).astype(int):np.sum(dimensions[0:indSMBarma + 1]).astype(int), np.sum(dimensions[0:indBDWarma]).astype(int):np.sum(dimensions[0:indBDWarma + 1]).astype(int)]
                    if md.smb.arma_timestep != md.basalforcings.arma_timestep and np.any(covsum != 0):
                        raise IOError('SMBarma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance')
            elif indTFarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indTFarma]).astype(int):np.sum(dimensions[0:indTFarma + 1]).astype(int), np.sum(dimensions[0:indBDWarma]).astype(int):np.sum(dimensions[0:indBDWarma + 1]).astype(int)]
                    if md.frontalforcings.arma_timestep != md.basalforcings.arma_timestep and np.any(covsum != 0):
                        raise IOError('FrontalForcingsRignotarma and BasalforcingsDeepwaterMeltingRatearma have different arma_timestep and non-zero covariance')
        elif indPwarma != -1:
            if indSdarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSdarma]).astype(int):np.sum(dimensions[0:indSdarma + 1]).astype(int), np.sum(dimensions[0:indPwarma]).astype(int):np.sum(dimensions[0:indPwarma + 1]).astype(int)]
                    if md.frontalforcings.sd_arma_timestep != md.hydrology.arma_timestep and np.any(covsum != 0):
                        raise IOError('FrontalForingsSubglacialDischargearma and hydrologyarmapw have different arma_timestep and non-zero covariance')
            elif indSMBarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSMBarma]).astype(int):np.sum(dimensions[0:indSMBarma + 1]).astype(int), np.sum(dimensions[0:indPwarma]).astype(int):np.sum(dimensions[0:indPwarma + 1]).astype(int)]
                    if md.smb.arma_timestep != md.hydrology.arma_timestep and np.any(covsum != 0):
                        raise IOError('SMBarma and hydrologyarmapw have different arma_timestep and non-zero covariance')
            elif indTFarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indTFarma]).astype(int):np.sum(dimensions[0:indTFarma + 1]).astype(int), np.sum(dimensions[0:indPwarma]).astype(int):np.sum(dimensions[0:indPwarma + 1]).astype(int)]
                    if md.frontalforcings.arma_timestep != md.hydrology.arma_timestep and np.any(covsum != 0):
                        raise IOError('FrontalForcingsRignotarma and hydrologyarmapw have different arma_timestep and non-zero covariance')
        elif indSdarma != -1:
            if indSMBarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSMBarma]).astype(int):np.sum(dimensions[0:indSMBarma + 1]).astype(int), np.sum(dimensions[0:indSdarma]).astype(int):np.sum(dimensions[0:indSdarma + 1]).astype(int)]
                    if md.smb.arma_timestep != md.frontalforcings.sd_arma_timestep and np.any(covsum != 0):
                        raise IOError('SMBarma and FrontalForcingsSubglacialDischargearma have different arma_timestep and non-zero covariance')
            elif indTFarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSdarma]).astype(int):np.sum(dimensions[0:indSdarma + 1]).astype(int), np.sum(dimensions[0:indTFarma]).astype(int):np.sum(dimensions[0:indTFarma + 1]).astype(int)]
                    if md.frontalforcings.sd_arma_timestep != md.frontalforcings.arma_timestep and np.any(covsum != 0):
                        raise IOError('FrontalForcingsSubglacialDischargearma and FrontalForcingsRignotarma have different arma_timestep and non-zero covariance')
        elif indSMBarma != -1:
            if indTFarma != -1:
                for ii in range(len(lsCovmats)):
                    covm = lsCovmats[ii]
                    covsum = covm[np.sum(dimensions[0:indSMBarma]).astype(int):np.sum(dimensions[0:indSMBarma + 1]).astype(int), np.sum(dimensions[0:indTFarma]).astype(int):np.sum(dimensions[0:indTFarma + 1]).astype(int)]
                    if md.smb.arma_timestep != md.frontalforcings.arma_timestep and np.any(covsum != 0):
                        raise IOError('SMBarma and FrontalForcingsRignotarma have different arma_timestep and non-zero covariance')

        md = checkfield(md, 'fieldname', 'stochasticforcing.isstochasticforcing', 'values', [0, 1])
        md = checkfield(md, 'fieldname', 'stochasticforcing.fields', 'numel', num_fields, 'cell', 1, 'values', self.supportedstochforcings())
        md = checkfield(md, 'fieldname', 'stochasticforcing.covariance', 'NaN', 1, 'Inf', 1, 'size', [size_tot, size_tot, numtcovmat])  # global covariance matrix
        md = checkfield(md, 'fieldname', 'stochasticforcing.stochastictimestep', 'NaN', 1,'Inf', 1, '>=', md.timestepping.time_step)
        md = checkfield(md, 'fieldname', 'stochasticforcing.randomflag', 'numel', [1], 'values', [0, 1])
        if(numtcovmat>1): #check the time steps at which each covariance matrix starts to be applie
            md = checkfield(md, 'fieldname', 'stochasticforcing.timecovariance', 'NaN', 1, 'Inf', 1, '>=',md.timestepping.start_time,'<=',md.timestepping.final_time,'size',[1,numtcovmat])  # global covariance matrix
        if checkdefaults:
            md = checkfield(md, 'fieldname', 'stochasticforcing.defaultdimension', 'numel', 1, 'NaN', 1, 'Inf', 1, '>', 0)
            md = checkfield(md, 'fieldname', 'stochasticforcing.default_id', 'Inf', 1, 'NaN', 1, '>=', 0, '<=', self.defaultdimension, 'size', [md.mesh.numberofelements])
        return md
    # }}}

    def extrude(self, md):  # {{{
        self.default_id = project3d(md,'vector',self.default_id,'type','element')
        return self
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        yts = md.constants.yts

        WriteData(fid, prefix, 'object', self, 'fieldname', 'isstochasticforcing', 'format', 'Boolean')
        if not self.isstochasticforcing:
            return md
        else:
            num_fields = len(self.fields)
            if self.stochastictimestep == 0:
                md.stochasticforcing.stochastictimestep = md.timestepping.time_step #by default: stochastictimestep set to ISSM time step
            # Check if hydroarmapw is used
            if((type(md.hydrology).__name__ == 'hydrologyarmapw') and md.transient.ishydrology==1):
                ispwHydroarma = 1
            else:
                ispwHydroarma = 0

            # Retrieve dimensionality of each field
            dimensions = self.defaultdimension * np.ones((num_fields))
            for ind, field in enumerate(self.fields):
                # Checking for specific dimensions
                if field == 'SMBarma':
                    dimensions[ind] = md.smb.num_basins
                elif field == 'FrontalForcingsRignotarma':
                    dimensions[ind] = md.frontalforcings.num_basins
                elif field == 'FrontalForcingsSubglacialDischargearma':
                    dimensions[ind] = md.frontalforcings.num_basins
                elif field == 'BasalforcingsDeepwaterMeltingRatearma':
                    dimensions[ind] = md.basalforcings.num_basins
                elif field == 'FrictionWaterPressure' and ispwHydroarma:
                    dimensions[ind] = md.hydrology.num_basins
            
            if(len(np.shape(self.covariance))==3):
                nrow,ncol,numtcovmat = np.shape(self.covariance)
                lsCovmats = []
                for ii in range(numtcovmat):
                    lsCovmats.append(self.covariance[:,:,ii])
                if(md.timestepping.interp_forcing==1):
                    print('WARNING: md.timestepping.interp_forcing is 1, but be aware that there is no interpolation between covariance matrices')
                    print('         the changes between covariance matrices occur at the time steps specified in md.stochasticforcing.timecovariance')
            elif(len(np.shape(self.covariance))==2):
                nrow,ncol = np.shape(self.covariance)
                numtcovmat = 1
                lsCovmats = [self.covariance]

            # Scaling covariance matrix (scale column-by-column and row-by-row)
            scaledfields = ['BasalforcingsDeepwaterMeltingRatearma','BasalforcingsSpatialDeepwaterMeltingRate','DefaultCalving', 'FloatingMeltRate', 'SMBarma', 'SMBforcing']  # list of fields that need scaling * 1 / yts
            tempcovariance2d = np.zeros((numtcovmat,nrow*ncol))
            # Loop over covariance matrices #
            for kk in range(numtcovmat):
                kkcov = lsCovmats[kk]
                # Loop over the fields #
                for i in range(num_fields):
                    if self.fields[i] in scaledfields:
                        inds = range(int(np.sum(dimensions[0:i])), int(np.sum(dimensions[0:i + 1])))
                        for row in inds:  # scale rows corresponding to scaled field
                            kkcov[row, :] = 1 / yts * kkcov[row, :]
                        for col in inds:  # scale columns corresponding to scaled field
                            kkcov[:, col] = 1 / yts * kkcov[:, col]
                # Save scaled covariance #
                for rr in range(nrow):
                    ind0 = rr*ncol
                    tempcovariance2d[kk,ind0:ind0+ncol] = np.copy(kkcov[rr,:])
                    
            # Set dummy default_id vector if defaults not used
            if np.any(np.isnan(self.default_id)):
                self.default_id = np.zeros(md.mesh.numberofelements)
            # Set dummy timecovariance vector if a single covariance matrix is used
            if(numtcovmat==1):
                self.timecovariance = np.array([md.timestepping.start_time])
            # Reshape dimensions as column array for marshalling
            dimensions = dimensions.reshape(1, len(dimensions))

            WriteData(fid, prefix, 'data', num_fields, 'name', 'md.stochasticforcing.num_fields', 'format', 'Integer')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'fields', 'format', 'StringArray')
            WriteData(fid, prefix, 'data', dimensions, 'name', 'md.stochasticforcing.dimensions', 'format', 'IntMat', 'mattype', 2)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'default_id', 'data', self.default_id - 1, 'format', 'IntMat', 'mattype', 2)  #12Nov2021 make sure this is zero-indexed!
            WriteData(fid, prefix, 'object', self, 'fieldname', 'defaultdimension', 'format', 'Integer')
            WriteData(fid, prefix, 'data', numtcovmat, 'name', 'md.stochasticforcing.num_timescovariance', 'format', 'Integer')
            WriteData(fid, prefix, 'data', tempcovariance2d, 'name', 'md.stochasticforcing.covariance', 'format', 'DoubleMat')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'timecovariance', 'format', 'DoubleMat', 'scale', yts)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'stochastictimestep', 'format', 'Double', 'scale', yts)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'randomflag', 'format', 'Boolean')
    # }}}

    def supportedstochforcings(self):  # {{{
        """Defines list of fields supported by the class md.stochasticforcing
        """
        list1 = self.structstochforcing()
        list1 = list1.keys()
        return list(list1)
    # }}}

    def structstochforcing(self):  # {{{
        """Defines dictionary with list of fields
           supported and corresponding md names
        """
        structure = {'BasalforcingsDeepwaterMeltingRatearma': 'linearbasalforcingsarma',
                     'BasalforcingsSpatialDeepwaterMeltingRate': 'spatiallinearbasalforcings',
                     'DefaultCalving': 'calving',
                     'FloatingMeltRate': 'basalforcings',
                     'FrictionWaterPressure': 'friction',
                     'FrictionWaterPressure': 'frictioncoulomb',
                     'FrictionWaterPressure': 'frictionschoof',
                     'FrontalForcingsRignotarma': 'frontalforcingsrignotarma',
                     'FrontalForcingsSubglacialDischargearma': 'frontalforcingsrignotarma',
                     'SMBarma': 'SMBarma',
                     'SMBforcing': 'SMBforcing'}
        return structure
    # }}}
