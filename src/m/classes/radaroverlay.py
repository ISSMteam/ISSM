from fielddisplay import fielddisplay


class radaroverlay(object):
    """
    RADAROVERLAY class definition

       Usage:
          radaroverlay = radaroverlay()
    """

    def __init__(self):  # {{{
        self.pwr = float('NaN')
        self.x = float('NaN')
        self.y = float('NaN')

    #set defaults
        self.setdefaultparameters()

    # }}}

    def __repr__(self):  # {{{
        string = '   radaroverlay parameters:'
        string = "%s\n%s" % (string, fielddisplay(self, 'pwr', 'radar power image (matrix)'))
        string = "%s\n%s" % (string, fielddisplay(self, 'x', 'corresponding x coordinates [m]'))
        string = "%s\n%s" % (string, fielddisplay(self, 'y', 'corresponding y coordinates [m]'))
        return string
    # }}}

    def setdefaultparameters(self):  # {{{
        return self
    # }}}
