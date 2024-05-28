function splittedstring = strsplit_strict(inpstr,delimiter)
%STRSPLIT_STRICT - split a string of delimiter separated values
%
%   Usage:
%      output = strsplit_strict(inpstr,delimiter)

%Check input arguments
if(nargin ~= 2)
    error('There is no argument defined');
end

%deblank string
deblank(inpstr);

%Get number of substrings
idx  = findstr(inpstr,delimiter);
if size(idx) == 0
    splittedstring = {inpstr};
else
    sz = size(idx,2);
    splittedstring = {};
    %Loop through string and itinerate from delimiter to delimiter
    for i = 1:sz
        strtpos = 1;
        endpos = idx(i)-1;
        if i ~= 1
            strtpos = idx(i-1)+1;
        end
        if i == sz
            endpos = size(inpstr,2); 
            splittedstring(i+1) = {inpstr(idx(i)+1 : endpos)};
            endpos = idx(i)-1;
        end
        splittedstring(i) = {inpstr(strtpos : endpos)};   
    end
end
