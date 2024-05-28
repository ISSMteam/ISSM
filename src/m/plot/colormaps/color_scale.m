function map = color_scale(n, theta, r, dir)
%COLOR_SCALE Colormap with luminance and hue ramps and constant chroma.
%   MAP = COLOR_SCALE(N, THETA, R, DIR) computes a colormap that works well 
%   on color displays and also works well when printed on a grayscale 
%   printer.
%   
%   The color map is computed using a simple path through L*a*b* space that
%   approximates a uniform ramp in the L* direction, and a semicircular
%   path in the a*-b* plane.
%
%   N is the number of colormap colors.  THETA is the angle (in degrees) in
%   the a*-b* plane of the first color.  THETA is measured clockwise from
%   the a* axis.  R is the radius of the semicircular path in the a*-b*
%   plane.  DIR is either 'cw', for a clockwise traversal, or 'ccw', for a
%   counterclockwise traversal.
%
%   All of the arguments are optional.  The default values are N = 256,
%   THETA = 0, R = 50, and DIR = 'cw'.
%
%   COLOR_SCALE requires the Image Processing Toolbox version 4 or later.
%
%   Example
%   -------
%   Display a Radon transform image with colormaps produced by color_scale.
%
%       I = zeros(100,100);
%       I(25:75, 25:75) = 1;
%       theta = 0:180;
%       [R,xp] = radon(I,theta);
%       imshow(R,[],'InitialMag','fit')
%       colormap(color_scale)
%
%       % Try it with different parameters.
%       colormap(color_scale(256,0,88,'ccw'))
%
%   See also COLOR_SCALE_TOOL.

%   Steve Eddins
%   $Revision$  $Date$

if nargin < 4
    dir = 'cw';
end

if nargin < 3
    r = 50;
end

if nargin < 2
    theta = 0;
end

if nargin < 1
    n = 256;
end

if strcmp(dir, 'cw')
    angle_offset = -pi/2;
else
    angle_offset = pi/2;
end

theta = pi * theta / 180;
theta_vec = linspace(theta, theta + angle_offset, n).';

a = r*cos(theta_vec);
b = r*sin(theta_vec);
L = linspace(0, 100, n).';

Lab = [L, a, b];

cform = makecform('lab2srgb');

map = applycform(Lab, cform);
