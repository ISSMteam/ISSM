import numpy as np

import math

from checkfield import checkfield
from fielddisplay import fielddisplay
from project3d import project3d
from WriteData import WriteData


class sampling(object):
    """sampling class definition

    Usage:
        sampling = sampling()
    """

    def __init__(self, *args):  # {{{
        self.kappa = np.nan
        self.tau = 0
        self.beta = np.nan
        self.phi = np.nan
        self.alpha = 0
        self.robin = 0
        self.seed = 0
        self.requested_outputs = []

        if len(args) == 0:
            self.setdefaultparameters()
        else:
            raise RuntimeError('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = '   Sampling parameters::\n'
        s += '      Parameters of PDE operator (kappa^2 I-Laplacian)^(alpha/2)(tau):\n'
        s += '{}\n'.format(fielddisplay(self, 'kappa', 'coefficient of the identity operator'))
        s += '{}\n'.format(fielddisplay(self, 'tau', 'scaling coefficient of the solution (default: 1.0)'))
        s += '{}\n'.format(fielddisplay(self, 'alpha', 'exponent in PDE operator, (default: 2.0, BiLaplacian covariance operator)'))

        s += '      Parameters of Robin boundary conditions nabla () \\cdot normvec + beta ():\n'
        s += '{}\n'.format(fielddisplay(self, 'robin', 'Apply Robin boundary conditions (1 if applied and 0 for homogenous Neumann boundary conditions) (default: 0)'))
        s += '{}\n'.format(fielddisplay(self, 'beta', 'Coefficient in Robin boundary conditions (to be defined for robin = 1)'))

        s += '      Parameters for first-order autoregressive process (X_t = phi X_{t-1} + noise) (if transient):\n'
        s += '{}\n'.format(fielddisplay(self, 'phi', 'Temporal correlation factor (|phi|<1 for stationary process, phi = 1 for random walk process) (default 0)'))

        s += '      Other parameters of stochastic sampler:\n'
        s += '{}\n'.format(fielddisplay(self, 'seed', 'Seed for pseudorandom number generator (given seed if >=0 and random seed if <0) (default: -1)'))
        s += '{}\n'.format(fielddisplay(self, 'requested_outputs', 'additional outputs requested (not implemented yet)'))

        return s
    # }}}

    def setdefaultparameters(self):  # {{{

        # Apply Robin boundary conditions
        self.robin = 0

        # Exponent in fraction SPDE (default: 2, biLaplacian covariance operator)
        self.alpha = 2 # Default

        # Seed for pseudorandom number generator (default: -1, for random seed)
        self.seed = -1

        # Default output
        self.requested_outputs = ['default']

        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        return []
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if ('SamplingAnalysis' not in analyses):
            return md

        md = checkfield(md, 'fieldname', 'sampling.kappa', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices], '>', 0)
        md = checkfield(md, 'fieldname', 'sampling.tau', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'sampling.robin', 'numel', 1, 'values', [0, 1])
        if md.sampling.robin:
            md = checkfield(md, 'fieldname', 'sampling.beta', 'NaN', 1, 'Inf', 1, 'size', [md.mesh.numberofvertices], '>', 0)
        end
        md = checkfield(md, 'fieldname', 'sampling.alpha', 'NaN', 1, 'Inf', 1, 'numel', 1, '>', 0)
        md = checkfield(md, 'fieldname', 'sampling.seed', 'NaN', 1, 'Inf', 1, 'numel', 1)
        md = checkfield(md, 'fieldname', 'sampling.requested_outputs', 'stringrow', 1)

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'kappa', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'tau', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'beta', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'phi', 'format', 'DoubleMat', 'mattype', 1)
        WriteData(fid, prefix, 'object', self, 'fieldname', 'alpha', 'format', 'Integer')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'robin', 'format', 'Boolean')
        WriteData(fid, prefix, 'object', self, 'fieldname', 'seed', 'format', 'Integer')

        # Process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        WriteData(fid, prefix, 'data', outputs, 'name', 'md.sampling.requested_outputs', 'format', 'StringArray')
    # }}}

    def setparameters(self, md, lc, sigma):  # {{{
        nu = self.alpha - 1
        KAPPA = pow((8 * nu), 0.5) / lc
        TAU = pow((math.gamma(nu) / math.gamma(self.alpha) * (4 *  np.pi) * pow(KAPPA, 2 * nu) * pow(sigma, 2)), 0.5)
        md.sampling.kappa = KAPPA * np.ones((md.mesh.numberofvertices, 1))
        md.sampling.tau = TAU

        return md
    # }}}
