function excludeoutliers(newcontourname,contourname,domainname)
%EXCLUDEOUTLIERS - exclude points of contour that are not within the domain
%contour. return new contours in a different file.
%
%        Usage: excludeoutliers('NewContour.exp','Contour.exp','DomainOutline.exp');
%
%
%   See also EXPMASTER, EXPDOC

contour=expread(contourname);

for i=1:length(contour),
	flags=ContourToNodes(contour(i).x,contour(i).y,domainname,0);
	contour(i).x=contour(i).x(find(flags));
	contour(i).y=contour(i).y(find(flags));
end

expwrite(contour,newcontourname);
