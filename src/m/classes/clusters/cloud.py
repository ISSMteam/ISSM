import subprocess

try:
    from cloud_settings import cloud_settings
except ImportError:
    print('You need cloud_settings.py to proceed, check presence and sys.path')

class cloud(object):
    """CLOUD cluster class definition

    Usage:
        cluster = cloud('name', 'astrid', 'np', 3)
        cluster = cloud('name', oshostname(), 'np', 3, 'login', 'username')
    """

    def __init__(self, *args):  # {{{
        self.name = ''
        self.login = ''
        self.np = 1
        self.codepath = ''
        self.executionpath = ''
        self.interactive = 0

        # Initialize cluster using user settings if provided
        try:
            self = cloud_settings(self)
        except NameError:
            print('cloud_settings.py not found, using default settings')

        # OK get other fields
        options = pairoptions(*args)
        self = options.AssignObjectFields(self)
    # }}}

    def __repr__(self):  # {{{
        # Display the object
        s = 'class \'{}\' object \'{}\' = \n'.format(type(self), 'self')
        s += '    name: {}\n'.format(self.name)
        s += '    login: {}\n'.format(self.login)
        s += '    np: {}\n'.format(self.np)
        s += '    codepath: {}\n'.format(self.codepath)
        s += '    executionpath: {}\n'.format(self.executionpath)
        s += '    interactive: {}\n'.format(self.interactive)
        return s
    # }}}

    def checkconsistency(self, md, solution, analyses):  # {{{
        if self.np < 1:
            md = md.checkmessage('number of processors should be at least 1')

        if np.isnan(self.np):
            md = md.checkmessage('number of processors should not be NaN!')

        return self
    # }}}

    def BuildQueueScript(self, dirname, modelname, solution, io_gather, isvalgrind, isgprof, isdakota, isoceancoupling):  # {{{
        # Write queuing script
        fid = open(modelname + '.queue', 'w')

        fid.write('#/bin/bash\n')
        fid.write('source {}{}\n'.format(self.codepath, '/../etc/environment.sh'))
        fid.write('cd {}/{}\n'.format(self.executionpath, dirname))
        fid.write('mpiexec -np {} -f /home/mpich2.hosts {}/issm.exe {} {}/{} {} 2> {}.errlog > /dev/stdout | tee {}.outlog\n'.format(self.np, self.codepath, solution, self.executionpath, dirname, modelname, modelname, modelname))
    # }}}

    def UploadQueueJob(self, modelname, dirname, filelist):  # {{{
        # Compress the files into one zip
        compressstring = 'tar -zcf {}.tar.gz'.format(dirname)
        for file in filelist:
            compressstring += ' {}'.format(file)
        subprocess.call(compressstring, shell=True)

        if isempty(self.login):
            raise Exception('cloud BuildQueueScript: login should be supplied!')

        #upload input files
        issmstscpout(self.name, self.executionpath, self.login, '{}.tar.gz'.format(dirname))
    # }}}

    def LaunchQueueJob(self, modelname, dirname, filelist, restart, batch):  # {{{
        if self.interactive:
            print('sending files to remote cluster. once done, please log into cluster and launch job')
            if not isempty(restart):
                launchcommand = 'cd {} && cd {}'.format(self.executionpath, dirname)
            else:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname)
        else:
            #Execute Queue job
            if not isempty(restart):
                launchcommand = 'cd {} && cd {} && qsub {}.queue'.format(self.executionpath, dirname, modelname)
            else:
                launchcommand = 'cd {} && rm -rf ./{} && mkdir {} && cd {} && mv ../{}.tar.gz ./ && tar -zxf {}.tar.gz && qsub {}.queue'.format(self.executionpath, dirname, dirname, dirname, dirname, dirname, modelname)
        issmstssh(self.name, self.login, launchcommand)
    # }}}

    def Download(self, dirname, filelist):  # {{{
        # Copy files from cluster to current directory
        directory = '{}/{}/'.format(self.executionpath, dirname)
        issmstscpin(self.name, self.login, directory, filelist)
    # }}}
