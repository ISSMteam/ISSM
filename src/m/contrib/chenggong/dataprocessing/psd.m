function [psdx, freq] = psd(x)
%psd - Power Spectral Density (using fft)
%
%   x: time series
%   
%   psdx: Power Spectral Density(right half plane, value doubled)
%   freq: frequency(non-negative half)
%

% length of data
Nx = length(x);
% fft
xdft = fft(x);
xdft = xdft(1:floor(Nx/2+1));
% power specturm
psdx = (1/(2*pi*Nx)) * abs(xdft).^2;
%  In order to conserve the total power, multiply all frequencies that 
%  occur in both sets - the positive and negative frequencies — by a factor of 2
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/Nx:pi;

end
