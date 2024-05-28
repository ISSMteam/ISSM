function [xi yi but] = exp_ginput(numclicks,options);
%EXP_GINPUT - equivalent to MATLAB's ginput function but with more options
%
%   Usage:
%      [xi yi] = exp_ginput(numclicks,options);

%ginputtype = getfieldvalue(options,'ginputtype','default');
ginputtype = getfieldvalue(options,'ginputtype','myginput');

switch ginputtype
	case 'default'
		[xi yi but] = ginput(numclicks);
	case 'myginput'
		[xi yi but] = myginput(numclicks,'arrow');
	otherwise
		error('not supported yet');
end
