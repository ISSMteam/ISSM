from fielddisplay import fielddisplay
from project3d import project3d
from checkfield import checkfield
from WriteData import WriteData


class damage(object):
    """
    DAMAGE class definition

       Usage:
          damage = damage()
    """

    def __init__(self, *args):  # {{{
        #damage:
        self.isdamage = 0.
        self.D = float('NaN')
        self.law = float('NaN')
        self.spcdamage = float('NaN')
        self.max_damage = float('NaN')

        #numerical
        self.stabilization = float('NaN')
        self.maxiter = float('NaN')
        self.elementinterp = ''

        #general parameters for evolution law:
        self.stress_threshold = float('NaN')
        self.stress_ubound = float('NaN')
        self.kappa = float('NaN')
        self.c1 = float('NaN')
        self.c2 = float('NaN')
        self.c3 = float('NaN')
        self.c4 = float('NaN')
        self.healing = float('NaN')
        self.equiv_stress = float('NaN')
        self.requested_outputs = []

        if not len(args):
            self.setdefaultparameters()
        else:
            raise RuntimeError("constructor not supported")
    # }}}

    def __repr__(self):  # {{{
        s = '   Damage:\n'
        s += "%s\n" % fielddisplay(self, "isdamage", "is damage mechanics being used? [0 (default) or 1]")
        if self.isdamage:
            s += "%s\n" % fielddisplay(self, "D", "damage tensor (scalar for now)")
            s += "%s\n" % fielddisplay(self, "law", "damage law ['0: analytical', '1: pralong']")
            s += "%s\n" % fielddisplay(self, "spcdamage", "damage constraints (NaN means no constraint)")
            s += "%s\n" % fielddisplay(self, "max_damage", "maximum possible damage (0 <=max_damage < 1)")
            s += "%s\n" % fielddisplay(self, "stabilization", "0: no stabilization, 1: artificial diffusion, 2: SUPG (not working), 4: flux corrected transport")
            s += "%s\n" % fielddisplay(self, "maxiter", "maximum number of non linear iterations")
            s += "%s\n" % fielddisplay(self, "elementinterp", "interpolation scheme for finite elements [''P1'', ''P2'']")
            s += "%s\n" % fielddisplay(self, "stress_threshold", "stress threshold for damage initiation (Pa)")
            s += "%s\n" % fielddisplay(self, "stress_ubound", "stress upper bound for damage healing (Pa)")
            s += "%s\n" % fielddisplay(self, "kappa", "ductility parameter for stress softening and damage [ > 1]")
            s += "%s\n" % fielddisplay(self, "c1", "damage parameter 1 ")
            s += "%s\n" % fielddisplay(self, "c2", "damage parameter 2 ")
            s += "%s\n" % fielddisplay(self, "c3", "damage parameter 3 ")
            s += "%s\n" % fielddisplay(self, "c4", "damage parameter 4 ")
            s += "%s\n" % fielddisplay(self, "healing", "damage healing parameter")
            s += "%s\n" % fielddisplay(self, "equiv_stress", "0: von Mises, 1: max principal")
            s += "%s\n" % fielddisplay(self, 'requested_outputs', 'additional outputs requested')

        return s
    # }}}

    def extrude(self, md):  # {{{
        if self.isdamage:
            self.D = project3d(md, 'vector', self.D, 'type', 'node')
            self.spcdamage = project3d(md, 'vector', self.spcdamage, 'type', 'node')
            return self
    # }}}

    def setdefaultparameters(self):  # {{{
        #damage parameters:
        self.isdamage = 0
        self.D = 0
        self.law = 0
        self.max_damage = 1 - 1e-5  #if damage reaches 1, solve becomes singular, as viscosity becomes nil
        #Type of stabilization used
        self.stabilization = 4
        #Maximum number of iterations
        self.maxiter = 100
        #finite element interpolation
        self.elementinterp = 'P1'
        #damage evolution parameters
        self.stress_threshold = 1.3e5
        self.kappa = 2.8
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        self.c4 = 0
        self.healing = 0
        self.equiv_stress = 0
        #output default:
        self.requested_outputs = ['default']

        return self
    # }}}

    def defaultoutputs(self, md):  # {{{
        if md.mesh.domaintype().lower() == '2dhorizontal':
            list = ['DamageDbar']
        else:
            list = ['DamageD']
        return list
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        md = checkfield(md, 'fieldname', 'damage.isdamage', 'numel', [1], 'values', [0, 1])
        if self.isdamage:
            md = checkfield(md, 'fieldname', 'damage.D', '>=', 0, '<=', self.max_damage, 'size', [md.mesh.numberofvertices])
            md = checkfield(md, 'fieldname', 'damage.max_damage', '<', 1, '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.law', 'numel', [1], 'values', [0, 1, 2, 3])
            md = checkfield(md, 'fieldname', 'damage.spcdamage', 'Inf', 1, 'timeseries', 1)
            md = checkfield(md, 'fieldname', 'damage.stabilization', 'numel', [1], 'values', [0, 1, 2, 4])
            md = checkfield(md, 'fieldname', 'damage.maxiter', ' >= 0', 0)
            md = checkfield(md, 'fieldname', 'damage.elementinterp', 'values', ['P1', 'P2'])
            md = checkfield(md, 'fieldname', 'damage.stress_threshold', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.stress_ubound', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.kappa', '>', 1)
            md = checkfield(md, 'fieldname', 'damage.healing', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.c1', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.c2', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.c3', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.c4', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.healing', '>=', 0)
            md = checkfield(md, 'fieldname', 'damage.equiv_stress', 'numel', [1], 'values', [0, 1])
            md = checkfield(md, 'fieldname', 'damage.requested_outputs', 'stringrow', 1)
        elif self.law != 0:
            if (solution == 'DamageEvolutionSolution'):
                raise RuntimeError('Invalid evolution law (md.damage.law) for a damage solution')

        return md
    # }}}

    def marshall(self, prefix, md, fid):  # {{{
        WriteData(fid, prefix, 'object', self, 'fieldname', 'isdamage', 'format', 'Boolean')
        if self.isdamage:
            WriteData(fid, prefix, 'object', self, 'fieldname', 'D', 'format', 'DoubleMat', 'mattype', 1)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'law', 'format', 'Integer')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'spcdamage', 'format', 'DoubleMat', 'mattype', 1, 'timeserieslength', md.mesh.numberofvertices + 1, 'yts', md.constants.yts)
            WriteData(fid, prefix, 'object', self, 'fieldname', 'max_damage', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'stabilization', 'format', 'Integer')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'maxiter', 'format', 'Integer')
            WriteData(fid, prefix, 'name', 'md.damage.elementinterp', 'data', self.elementinterp, 'format', 'String')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_threshold', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'stress_ubound', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'kappa', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'c1', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'c2', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'c3', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'c4', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'healing', 'format', 'Double')
            WriteData(fid, prefix, 'object', self, 'fieldname', 'equiv_stress', 'format', 'Integer')

    #process requested outputs
        outputs = self.requested_outputs
        indices = [i for i, x in enumerate(outputs) if x == 'default']
        if len(indices) > 0:
            outputscopy = outputs[0:max(0, indices[0] - 1)] + self.defaultoutputs(md) + outputs[indices[0] + 1:]
            outputs = outputscopy
        if self.isdamage:
            WriteData(fid, prefix, 'data', outputs, 'name', 'md.damage.requested_outputs', 'format', 'StringArray')
    # }}}
