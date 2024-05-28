function STARCLUSTER=starcluster()
%STARCLUSTER - Get path to STARCLUSTER command
%
%   Usage:
%      STARCLUSTER=STARCLUSTER()

STARCLUSTER=[issmdir() '/externalpackages/python/install/bin/starcluster -c ' jplsvn '/proj-group/CloudComputing/starcluster.config'];
