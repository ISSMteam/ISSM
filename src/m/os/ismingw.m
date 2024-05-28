function result=ismingw()
%ISMINGW - Returns 1 if machine is running MinGW.
%
%   Usage:
%      if ismingw,
%         [...]
%      end
%

[status,result]=system('uname -rs | grep "MINGW" | wc -l');
result=str2num(result);
