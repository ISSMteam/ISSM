function [mask]=nunataks(mask)
%NUNATAKS - bias mask towards increased 0 coverage
% 
%  mask is an image of arbitrary size, format binary, with values 1 for foreground, and 0 for background
% 
%  Usage:   mask=nunataks(mask)
%           [mask]=aggregation(mask);
%
%  See also CLOSING, OPENING, DILATION, EROSION, AGGREGATION

rocks=~mask;

%matrices for convolution: 
matrix=[0 1 0; 1 0 1; 0 1 0];

%do not exist, i.e. locations that stand pretty much alone. 
mask=filter2(matrix1,mask,'same');
pos=find(~crocks & rocks);

mask(pos)=0;
