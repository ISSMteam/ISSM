%
%  function to find and read the first or last positive integer
%  in a character string.
%
%  function [aint]=str2int(astr,cfl)
%
function [aint]=str2int(astr,cfl)

aint=[];

if     ~exist('cfl','var') || strncmpi(cfl,'f',1)
    i=1;

    while (i <= length(astr))
        if (astr(i) >= '0' && astr(i) <= '9')
            aint=sscanf(astr(i:length(astr)),'%d',[1,1]);
            return
        else
            i=i+1;
        end
	end

elseif strncmpi(cfl,'l',1)
    i=length(astr);
    ifound=false;

    while (i >= 1)
        if     (astr(i) >= '0' && astr(i) <= '9')
            ifound=true;
            i=i-1;
        elseif ~ifound
            i=i-1;
        else
            aint=sscanf(astr(i+1:length(astr)),'%d',[1,1]);
            return
        end
	end

    if ifound
        aint=sscanf(astr,'%d',[1,1]);
        return
    end
end

end
