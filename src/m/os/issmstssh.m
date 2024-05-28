function issmstssh(host,login,command)
%ISSMSSH - wrapper for OS independent ssh command.
%
%   usage: 
%      issmstssh(host,command)

%just use starcluster command to pipe an ssh command through
system([starcluster() ' sshmaster ' host ' --user ' login ' ''' command '''']);

