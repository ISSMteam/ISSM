function issmstssh(host,login,command)
%ISSMSTSSH - wrapper for OS independent ssh command, run through starcluster
%
%   Usage:
%      issmstssh(host,login,command)

%just use starcluster command to pipe an ssh command through
system([starcluster() ' sshmaster ' host ' --user ' login ' ''' command '''']);

