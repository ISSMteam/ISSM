function w2 = wfDistance(x, y)
%wfDistance - compute Wasserstein-Fourier Distance, according to Cazelles, E.,
%         Robert, A. & Tobar, F. The Wasserstein-Fourier Distance for 
%         Stationary Time Series. Ieee T Signal Proces 69, 709–721 (2019).
%
%   x: first time series
%   y: second time series
%
%   w2:  Wasserstein-Fourier Distance, or W2 distance of sx and sy
%
N = length(x);
Sinv=linspace(0,1,N);

[csx, ~, freqx] = npsd(x);
[csy, ~, freqy] = npsd(y);

% inverse cummulative function
[~, qx] = inverseHistograms(csx, freqx, Sinv);
[~, qy] = inverseHistograms(csy, freqy, Sinv);

%
w2 = sqrt(trapz(Sinv ,(qx-qy).^2));

end
