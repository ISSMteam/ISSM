import os.path
from collections import OrderedDict
import pairoptions
from loadvars import loadvars
from loadmodel import loadmodel
from savevars import savevars
from model import model
#hack to keep python 2 compatipility
try:
    #py3 import
    import dbm
except ImportError:
    #py2 import
    from whichdb import whichdb

import MatlabFuncs as m


class organizer(object):
    """
    ORGANIZER class definition

       Supported options:
          repository: directory where all models will be saved
          prefix:     prefix for saved model names
          steps:      requested steps
          trunkprefix:prefix of previous run with a different prefix. Used to branch.

       Usage:
          org = organizer(varargin)

       Examples:
          org = organizer('repository', 'Models/', 'prefix', 'AGU2015', 'steps', 0);  %build an empty organizer object with a given repository
    """

    def __init__(self, *args):  # {{{
        self._currentstep = 0
        self.repository = './'
        self.prefix = 'model.'
        self.trunkprefix = ''
        self.steps = []
        self.requestedsteps = [0]

        #process options
        options = pairoptions.pairoptions(*args)

        #Get prefix
        prefix = options.getfieldvalue('prefix', 'model.')
        if not isinstance(prefix, str):
            raise TypeError("prefix is not a string")
        if not m.strcmp(prefix, prefix.strip()) or len(prefix.split()) > 1:
            raise TypeError("prefix should not have any white space")
        self.prefix = prefix

        #Get repository
        repository = options.getfieldvalue('repository', './')
        if not isinstance(repository, str):
            raise TypeError("repository is not a string")
        if not os.path.isdir(repository):
            raise IOError("Directory '%s' not found" % repository)
        self.repository = repository

        #Get steps
        self.requestedsteps = options.getfieldvalue('steps', [0])

        #Get trunk prefix (only if provided by user)
        if options.exist('trunkprefix'):
            trunkprefix = options.getfieldvalue('trunkprefix', '')
            if not isinstance(trunkprefix, str):
                raise TypeError("trunkprefix is not a string")
            if not m.strcmp(trunkprefix, trunkprefix.strip()) or len(trunkprefix.split()) > 1:
                raise TypeError("trunkprefix should not have any white space")
            self.trunkprefix = trunkprefix
    # }}}

    def __repr__(self):  # {{{
        s = ""
        s += "%s\n" % "   Repository: '%s'" % self.repository
        s += "%s\n" % "   Prefix:     '%s'" % self.prefix
        if not self.steps:
            s += "%s\n" % "   no step"
        else:
            for step in self.steps:
                s += "%s\n" % "   step  #%2i: '%s'", step['id'], step['string']
    # }}}

    def load(self, string):  # {{{
        #Get model path
        if not isinstance(string, str):
            raise TypeError("argument provided is not a string")
        path = os.path.join(self.repository, self.prefix + string)

        #figure out if the model is there
        if os.path.exists(path):
            struc = loadvars(path)
            name = name = [key for key in list(struc.keys())]
            md = struc.name[0]
        else:
            raise IOError("Could not find '%s'" % path)

        return md
    # }}}

    def loadmodel(self, string):  # {{{
        #Get model path
        if not isinstance(string, str):
            raise TypeError("argument provided is not a string")
        path1 = os.path.join(self.repository, self.prefix + string + '.python')
        path2 = os.path.join(self.repository, string)

        #figure out if the model is there, otherwise, we have to use the default path supplied by user.
        if dbm.whichdb(path1):
            md = loadmodel(path1)
            return md
        elif dbm.whichdb(path2):
            md = loadmodel(path2)
            return md

        #If we are here, the model has not been found. Try trunk prefix if provided
        if self.trunkprefix:
            path2 = os.path.join(self.repository, self.trunkprefix + string)
            if not os.path.exists(path2):
                raise IOError("Could find neither '%s' nor '%s'" % (path1, path2))
            else:
                print(("--> Branching '%s' from trunk '%s'" % (self.prefix, self.trunkprefix)))
                md = loadmodel(path2)
                return md
        else:
            raise IOError("Could not find '%s'" % path1)
    # }}}

    def perform(self, string):  # {{{
        bool = False

        #Some checks
        if not isinstance(string, str):
            raise TypeError("Step provided should be a string")
        if not m.strcmp(string, string.strip()) or len(string.split()) > 1:
            raise TypeError("Step provided should not have any white space")
        if self._currentstep > 0 and string in [step['string'] for step in self.steps]:
            raise RuntimeError("Step '%s' already present. Change name" % string)

        #Add step
        self.steps.append(OrderedDict())
        self.steps[-1]['id'] = len(self.steps)
        self.steps[-1]['string'] = string
        self._currentstep += 1

        #if requestedsteps = 0, print all steps in self
        if 0 in self.requestedsteps:
            if self._currentstep == 1:
                print(("   prefix: %s" % self.prefix))
            print(("   step  #%i : %s" % (self.steps[self._currentstep - 1]['id'], self.steps[self._currentstep - 1]['string'])))

        #Ok, now if _currentstep is a member of steps, return true
        if self._currentstep in self.requestedsteps:
            print(("\n   step  #%i : %s\n" % (self.steps[self._currentstep - 1]['id'], self.steps[self._currentstep - 1]['string'])))
            bool = True

            #last check: is this step locked?
            s = self.steps[self._currentstep - 1]['string']
            if len(s) > 7:
                if s[-6:] == 'Locked':
                    raise Exception('organizer: you are trying to run a locked step! Unlock it first!')

        return bool
    # }}}

    def savemodel(self, md, name='default'):  # {{{
        #check
        if self._currentstep == 0:
            raise RuntimeError("Cannot save model because organizer (org) is empty! Make sure you did not skip any perform call")
        if self._currentstep > len(self.steps):
            raise RuntimeError("Cannot save model because organizer (org) is not up to date!")

        if (name == 'default'):
            name = os.path.join(self.repository, self.prefix + self.steps[self._currentstep - 1]['string'] + '.python')
        else:
            name = os.path.join(self.repository, name)
        print(("saving model as: '%s'" % name))

    #check that md is a model
        if not isinstance(md, model):
            print("second argument is not a model")
        if self._currentstep > len(self.steps):
            raise RuntimeError("organizer error message: element with id %d not found" % self._currentstep)

    #save model
        savevars(name, 'md', md)
    # }}}
