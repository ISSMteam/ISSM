from fielddisplay import fielddisplay
from iluasmoptions import iluasmoptions
from IssmConfig import IssmConfig
from issmgslsolver import issmgslsolver
from issmmumpssolver import issmmumpssolver
from mumpsoptions import mumpsoptions


class toolkits(object):
    """toolkits class definition

    Usage:
        self = toolkits()
    """

    def __init__(self, *args):  # {{{
        self.DefaultAnalysis = None
        self.RecoveryAnalysis = None

        nargs = len(args)
        if nargs == 0:
            self.setdefaultparameters()
        elif nargs == 1:
            # TODO: Replace the following with constructor
            self.setdefaultparameters()
        else:
            raise Exception('constructor not supported')
    # }}}

    def __repr__(self):  # {{{
        s = "List of toolkits options per analysis:\n\n"
        for analysis in list(vars(self).keys()):
            s += "{}\n".format(fielddisplay(self, analysis, ''))

        return s
    # }}}

    def addoptions(self, analysis, *args):  # {{{
        """addoptions - add analysis to md.toolkits.analysis

        Optional third parameter adds toolkits options to analysis.
        
        Usage:
            md.toolkits = addoptions(md.toolkits, 'StressbalanceAnalysis', FSoptions())
            md.toolkits = addoptions(md.toolkits, 'StressbalanceAnalysis')
        """

        # Create dynamic property if property does not exist yet
        if not hasattr(self, analysis):
            setattr(self, analysis, None)

        # Add toolkits options to analysis
        if len(args) == 1:
            setattr(self, analysis, args[0])

        return self
    # }}}

    def setdefaultparameters(self):  # {{{
        # Default toolkits
        if IssmConfig('_HAVE_PETSC_')[0]:
            # MUMPS is the default toolkits
            if IssmConfig('_HAVE_MUMPS_')[0]:
                self.DefaultAnalysis = mumpsoptions()
            else:
                self.DefaultAnalysis = iluasmoptions()
        else:
            if IssmConfig('_HAVE_MUMPS_')[0]:
                self.DefaultAnalysis = issmmumpssolver()
            elif IssmConfig('_HAVE_GSL_')[0]:
                self.DefaultAnalysis = issmgslsolver()
            else:
                raise IOError('ToolkitsFile error: need at least MUMPS or GSL to define ISSM solver type, no default solver assigned')

        # Use same solver for Recovery mode
        self.RecoveryAnalysis = self.DefaultAnalysis

        return self
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        supported_analyses = [
            'DefaultAnalysis',
            'RecoveryAnalysis',
            'StressbalanceAnalysis',
            'StressbalanceVerticalAnalysis',
            'GLheightadvectionAnalysis',
            'MasstransportAnalysis',
            'ThermalAnalysis',
            'EnthalpyAnalysis',
            'AdjointBalancethicknessAnalysis',
            'BalancethicknessAnalysis',
            'Balancethickness2Analysis',
            'BalancethicknessSoftAnalysis',
            'BalancevelocityAnalysis',
            'DamageEvolutionAnalysis',
            'LoveAnalysis',
            'EsaAnalysis',
            'SealevelchangeAnalysis',
            'FreeSurfaceBaseAnalysis',
            'FreeSurfaceTopAnalysis',
            'LevelsetAnalysis',
            'DebrisAnalysis',
            'L2ProjectionBaseAnalysis',
            'ExtrudeFromBaseAnalysis',
            'ExtrudeFromTopAnalysis'
        ]
        analyses = list(vars(self).keys())
        for analysis in analyses:
            if analysis not in supported_analyses:
                md.checkmessage('md.toolkits.{} not supported yet'.format(analysis))

            if not getattr(self, analysis):
                md.checkmessage('md.toolkits.{} is empty'.format(analysis))

        return md
    # }}}

    def ToolkitsFile(self, filename):  # {{{
        """ToolkitsFile - build toolkits file

        Build a PETSc compatible options file, from the toolkits model field and return options string.
        This file will also be used when the toolkit used is 'issm' instead of 'petsc'.

        Usage:
            ToolkitsFile(toolkits, filename)
        """

        # Open file for writing
        try:
            fid = open(filename, 'w')
        except IOError as e:
            raise IOError('ToolkitsFile error: could not open {} for writing due to {}'.format(filename), e)

        # Write header
        fid.write('{}{}{}\n'.format('%Toolkits options file: ', filename, ' written from Python toolkits array'))

        # Start writing options
        for analysis in list(vars(self).keys()):
            options = getattr(self, analysis)

            # First write analysis
            fid.write('\n+{}\n'.format(analysis))  # Append a + to recognize it's an analysis enum

            # Now, write options
            for optionname, optionvalue in list(options.items()):

                if not optionvalue:
                    # This option has only one argument
                    fid.write('-{}\n'.format(optionname))
                else:
                    # Option with value. Value can be string or scalar.
                    if isinstance(optionvalue, (bool, int, float)):
                        fid.write('-{} {}\n'.format(optionname, optionvalue))
                    elif isinstance(optionvalue, str):
                        fid.write('-{} {}\n'.format(optionname, optionvalue))
                    else:
                        raise TypeError('ToolkitsFile error: option {} is not well formatted'.format(optionname))

        fid.close()
    # }}}
