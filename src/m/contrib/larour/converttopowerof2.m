function converttopowerof2(tiffname,pngname)
%CONVERTTOPOWEROF2: read it a tiff, resize it so its xy dimensions are multiples of 2, and convert to png 
%
%  Usage:    converttopowerof2('temp.tif','temp.png')
%

	if ismac,
		dyld_library_path_old=getenv('DYLD_LIBRARY_PATH');
		setenv('DYLD_LIBRARY_PATH','/opt/local/lib:/usr/lib');
	end

	%figure out the size of the tiff
	[status,width]=system(['tiffinfo ' tiffname ' 2>/dev/null | command grep "Image Width" | awk ''{printf("%s\n",$3);}''']); 
	[status,length]=system(['tiffinfo ' tiffname ' 2>/dev/null | command grep "Image Width" | awk ''{printf("%s\n",$6);}''']);
	width=str2num(width); length=str2num(length);

	
	%Now, figure out the highest multiple of 2 for both width and length:
	width=2^nextpow2(width); length=2^nextpow2(length);

	%make sure the width and length are < 2000: 
	if width>2^11, width=2^11; end
	if length>2^11, length=2^11; end

	%convert image to that size: 
	setenv('DYLD_LIBRARY_PATH','/opt/local/lib:/usr/lib');
	
	[status,result]=system(sprintf('convert %s -resize %ix%i! %s',tiffname,width,length,pngname));
	system(sprintf('rm -rf %s',tiffname));

	%reset DYLD_LIBRARY_PATH to what it was:
	if ismac, setenv('DYLD_LIBRARY_PATH',dyld_library_path_old); end
end
