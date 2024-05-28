function gmtOce(num) 

%---------------------------------------------------------------------
% gmtOce :: a function to write a GMT code for selecting pixels (e.g., ocean function) 
%---------------------------------------------------------------------
% This code is written as a part of the ISSM-PSL project
% (c) S. Adhikari 
%     Jet Propulsion Laboratory, Caltech 
%     November 3, 2014
%---------------------------------------------------------------------

text0 = ''; 
text1 = '#!/bin/sh'; 
text2 = '# GMT code for selecting ocean and continental data (from Selen)'; 
text3 = '# This code is written as a part of the ISSM-PSL project'; 
text4 = '# (c) Surendra Adhikari'; 
text5 = '#     Jet Propulsion Laboratory, Caltech'; 
text6 = '#     November 3, 2014'; 
text7 = sprintf('pixAll="./gmtdir/pixels/pix%dall.txt"',num); 
text8 = sprintf('pixOce="./gmtdir/pixels/pix%doce.txt"',num); 
text9 = sprintf('pixCon="./gmtdir/pixels/pix%dcon.txt"',num); 
text10 = '# Select ocean (and continent) pixels'; 
text11 = 'gmt gmtselect $pixAll -h0 -Df -R0/360/-90/90  -A0 -JQ180/200 -Nk/s/k/s/k > $pixOce'; 
text12 = '#gmt gmtselect $pixAll -h0 -Df -R              -A0 -JQ        -Ns/k/s/k/s > $pixCon'; 
% 
fid = fopen(sprintf('./gmtdir/oceanFunc.sh'),'w');  % "./" = ISSM_PSL directory! 
fprintf(fid,'%s\n',text1); 
fprintf(fid,'%s\n',text0); 
fprintf(fid,'%s\n',text2); 
fprintf(fid,'%s\n',text3); 
fprintf(fid,'%s\n',text4); 
fprintf(fid,'%s\n',text5); 
fprintf(fid,'%s\n',text6); 
fprintf(fid,'%s\n',text0); 
fprintf(fid,'%s\n',text7); 
fprintf(fid,'%s\n',text8); 
fprintf(fid,'%s\n',text9); 
fprintf(fid,'%s\n',text0); 
fprintf(fid,'%s\n',text10); 
fprintf(fid,'%s\n',text11); 
fprintf(fid,'%s\n',text12); 
fclose(fid); 

