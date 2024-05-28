import numpy as np
from project3d import project3d
from fielddisplay import fielddisplay
from checkfield import checkfield
from WriteData import WriteData


class misfit(object):
    """
    MISFIT class definition

    Usage:
        misfit = misfit()
        misfit = misfit(name = 'SurfaceAltimetry',
                    definitionstring = 'Outputdefinition1',
            model_string = 'Surface',
                    observation_string = 'SurfaceObservations',
                     observation = md.geometry.surface,
                        timeinterpolation = 'nearestneighbor',
                        local = 1,
                        weights = np.ones((md.mesh.numberofvertices, 1)),
                        weights_string = 'WeightsSurfaceObservations')
    """

    def __init__(self, name=None, definitionstring=None, model_string=None, observation=None, observation_string=None, timeinterpolation=None, local=None, weights=None, weights_string=None, cumulated=None):  # {{{
        self.name = name if name is not None else ''
        #string that identifies this output definition uniquely, from 'Outputdefinition[1 - 100]'
        self.definitionstring = definitionstring if definitionstring is not None else ''
        #string for field that is modeled
        self.model_string = model_string if model_string is not None else ''
        #observed field that we compare the model against
        self.observation = observation if observation is not None else float('NaN')
        #string for observed field.
        self.observation_string = observation_string if observation_string is not None else ''
        self.timeinterpolation = timeinterpolation if timeinterpolation is not None else 'nearestneighbor'
        self.local = local if local is not None else 1
        #weight coefficients for every vertex
        self.weights = weights if weights is not None else float('NaN')
        #string to identify this particular set of weights
        self.weights_string = weights_string if weights_string is not None else ''
        #do we cumulate misfit through time?
        self.cumulated = cumulated if cumulated is not None else float('NaN')
    # }}}

    def __repr__(self):  # {{{
        string = '   Misfit:'

        string = "%s\n%s" % (string, fielddisplay(self, 'name', 'identifier for this misfit response'))
        string = "%s\n%s" % (string, fielddisplay(self, 'definitionstring', 'string that identifies this output definition uniquely, from "Outputdefinition[1 - 10]"'))
        string = "%s\n%s" % (string, fielddisplay(self, 'model_string', 'string for field that is modeled'))
        string = "%s\n%s" % (string, fielddisplay(self, 'observation', 'observed field that we compare the model against'))
        string = "%s\n%s" % (string, fielddisplay(self, 'observation_string', 'observation string'))
        string = "%s\n%s" % (string, fielddisplay(self, 'local', 'is the response local to the elements, or global? (default is 1)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'timeinterpolation', 'interpolation routine used to interpolate misfit between two time steps (default is "nearestneighbor"'))
        string = "%s\n%s" % (string, fielddisplay(self, 'weights', 'weights (at vertices) to apply to the misfit'))
        string = "%s\n%s" % (string, fielddisplay(self, 'weights_string', 'string for weights for identification purposes'))
        return string
    # }}}

    def extrude(self, md):  # {{{
        if not np.any(np.isnan(self.weights)):
            self.weights = project3d(md, 'vector', self.weights, 'type', 'node')
        if not np.any(np.isnan(self.observation)):
            self.observation = project3d(md, 'vector', self.observation, 'type', 'node')
        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if type(self.name) != str:
            raise TypeError('misfit error message: "name" field should be a string!')

        OutputdefinitionStringArray = []
        for i in range(100):
            OutputdefinitionStringArray.append('Outputdefinition' + str(i))

        md = checkfield(md, 'fieldname', 'self.definitionstring', 'field', self.definitionstring, 'values', OutputdefinitionStringArray)
        if type(self.timeinterpolation) != str:
            raise TypeError('misfit error message: "timeinterpolation" field should be a string!')

        md = checkfield(md, 'fieldname', 'self.observation', 'field', self.observation, 'timeseries', 1, 'NaN', 1, 'Inf', 1)
        md = checkfield(md, 'fieldname', 'self.timeinterpolation', 'field', self.timeinterpolation, 'values', ['nearestneighbor'])
        md = checkfield(md, 'fieldname', 'self.weights', 'field', self.weights, 'timeseries', 1, 'NaN', 1, 'Inf', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  #  {{{
        WriteData(fid, prefix, 'data', self.name, 'name', 'md.misfit.name', 'format', 'String')
        WriteData(fid, prefix, 'data', self.definitionstring, 'name', 'md.misfit.definitionstring', 'format', 'String')
        WriteData(fid, prefix, 'data', self.model_string, 'name', 'md.misfit.model_string', 'format', 'String')
        WriteData(fid, prefix, 'data', self.observation, 'name', 'md.misfit.observation', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'data', self.observation_string, 'name', 'md.misfit.observation_string', 'format', 'String')
        WriteData(fid, prefix, 'data', self.local, 'name', 'md.misfit.local', 'format', 'Integer')
        WriteData(fid, prefix, 'data', self.timeinterpolation, 'name', 'md.misfit.timeinterpolation', 'format', 'String')
        WriteData(fid, prefix, 'data', self.weights, 'name', 'md.misfit.weights', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
        WriteData(fid, prefix, 'data', self.weights_string, 'name', 'md.misfit.weights_string', 'format', 'String')
    # }}}
