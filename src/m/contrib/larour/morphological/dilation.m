function mask=dilation(mask,neighboorhood)
%deletion algorithm using 4 pixel neighboors.

%convolve: 

matrix=ones(3,3); 

%4 neighboorhood: 
%corners
if neighboorhood==4,
	matrix(1,1)=0;
	matrix(1,3)=0;
	matrix(3,1)=0;
	matrix(3,3)=0;
end
%center
matrix(2,2)=0;

%convolve mask: 
convol=filter2(matrix,mask,'same');

pos=find(mask==0);
pos2=find(convol(pos)~=0);

mask(pos(pos2))=1;
