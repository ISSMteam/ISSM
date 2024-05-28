function surface = f_surface_triangle(vec1, vec2, vec3)

% Copyright 2010-2011 by Lee Samuel Finn.
% This program file is part of MEALPix. It is licensed under the Apache
% License, Version 2.0 (the  "License"); you may not use MEALPix except in
% compliance with the License. You may obtain a copy of the License at
% <http://www.apache.org/licenses/LICENSE-2.0>
%
% Unless required by applicable law or agreed to in writing, MEALPix
% software distributed under the License is distributed on an "AS IS"
% BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
% implied. See the License for the specific language governing permissions
% and limitations under the License.

% $Id: f_surface_triangle.m 5795 2011-02-20 04:01:32Z lsf@GRAVITY.PSU.EDU $

% use pix_tools
%=======================================================================
% returns the surface in steradians
%  of the spherical triangle with vertices vec1, vec2, vec3
%
% algorithm : finds triangle sides and uses l'Huilier formula to compute
% 'spherical excess' = surface area of triangle on a sphere of radius one
% see, eg Bronshtein, Semendyayev Eq 2.86
%=======================================================================
%   real(kind=dp), dimension(1:3) :: v1, v2, v3
% persistent side x0 x1 x2 x3 ;
%
% if isempty(side), side=zeros(1,3); end;
% %  real(kind=dp) :: hp
% if isempty(x0), x0=0; end;
% if isempty(x1), x1=0; end;
% if isempty(x2), x2=0; end;
% if isempty(x3), x3=0; end;
%=======================================================================
% half perimeter
%   hp = 0.5 * (side1 + side2 + side3)
%   ! l'Huilier formula
%   x0 = tan( hp          * 0.5)
%   x1 = tan((hp - side1) * 0.5)
%   x2 = tan((hp - side2) * 0.5)
%   x3 = tan((hp - side3) * 0.5)
% find triangle sides
side(1)=f_angdist(vec2, vec3);
side(2)=f_angdist(vec3, vec1);
side(3)=f_angdist(vec1, vec2);
% divide by 4
side([1:3]) = side([1:3]) .* 0.25;
% l'Huilier formula
x0 = tan( side(1) + side(2) + side(3) );
x1 = tan(-side(1) + side(2) + side(3) );
x2 = tan( side(1) - side(2) + side(3) );
x3 = tan( side(1) + side(2) - side(3) );
surface = 4.0 .* atan( sqrt(x0 .* x1 .* x2 .* x3) );

return