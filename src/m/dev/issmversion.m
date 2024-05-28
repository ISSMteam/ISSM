function vers = issmversion(),
%ISSMVERSION - display ISSM version
%
%   Usage:
%      issmversion()
%      version = issmversion()


if exist('IssmConfig_matlab')~=3,
	error('ISSM not correctly installed. "IssmConfig_matlab" not found');
end

if nargout==1
	vers = IssmConfig('PACKAGE_VERSION');
	return;
end

disp([' ']);
disp([IssmConfig('PACKAGE_NAME') ' Version ' IssmConfig('PACKAGE_VERSION')]);
disp(['(website: ' IssmConfig('PACKAGE_URL') ' contact: ' IssmConfig('PACKAGE_BUGREPORT') ')']);
disp([' ']);
disp(['Build date: ' IssmConfig('PACKAGE_BUILD_DATE')]);
disp(['Compiled on ' IssmConfig('HOST_VENDOR') ' ' IssmConfig('HOST_OS') ' ' IssmConfig('HOST_ARCH') ' by ' IssmConfig('USER_NAME')]);
disp([' ']);
disp(['Copyright (c) 2009-2024 California Institute of Technology']);
disp([' ']);
disp(['    to get started type: issmdoc']);
disp([' ']);
