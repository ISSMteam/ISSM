function expxy2ll(filename,sgn,central_meridian,standard_parallel)  
%EXPLL2XY: switch exp argus file from lat,long to x,y
%   Usage:
%      expxy2ll(filename,sgn,central_meridian,standard_parallel)
%      - sgn = Sign of latitude +1 : north latitude (default is mer=45 lat=70)
%                               -1 : south latitude (default is mer=0  lat=71)

%Get central_meridian and standard_parallel depending on hemisphere
if nargin==4,
	delta = central_meridian;
	slat  = standard_parallel;
elseif nargin==2
	if sgn == 1,
		delta = 45; slat = 70;
		disp('Info: creating coordinates in polar stereographic (Std Latitude: 70ºN Meridian: 45º)');
	elseif sgn==-1,
		delta = 0;  slat = 71;
		disp('Info: creating coordinates in polar stereographic (Std Latitude: 71ºS Meridian: 0º)');
	else
		error('Sign should be either +1 or -1');
	end
else
	help expxy2ll
	error('bad usage');
end

%read filename: 
domain=expread(filename);

%change to x,y: 
for i=1:length(domain),
	[domain(i).y domain(i).x]= xy2ll(domain(i).x,domain(i).y,sgn,delta,slat); %watch out to swap lat and long
end

%write back to filename: 
expwrite(domain,filename);
