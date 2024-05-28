function mask=erosion(mask,neighboorhood)
%erosion algorithm using neighboorhood pixel neighboors.

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

pos=find(mask==1);
pos2=find(convol(pos)<neighboorhood);

mask(pos(pos2))=0;
