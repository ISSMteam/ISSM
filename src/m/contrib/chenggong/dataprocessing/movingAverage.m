function outdata =  movingAverage(rawdata, varargin)
%
% rawdata is an MxN matrix, where the first M-1 rows are the raw data, the row M is the time

options    = pairoptions(varargin{:});
timeWindow	= getfieldvalue(options, 'time window', 30); % unit in days
resamp = getfieldvalue(options, 'resample', 1); % resample the data set with the size of window

% get data info {{{
[M, N] = size(rawdata);
if M < 2
	error('The data has to be at least two rows!');
end
if N < 2
	error('The data has to be at least two columns!');
end

% sort the time
time = rawdata(M, :);
[time, tI] = sort(time);
data = rawdata(1:M-1, :);
data = data(:, tI);

dt = time(2) - time(1);
if dt < 0
	error('The time need to be in ascending order.');
end
windowSize = ceil(timeWindow/365/dt);
%}}}
%% filter {{{
if resamp
	ind = [1:timeWindow:N];
else
	ind = [1:N];
end
b = (1/windowSize)*ones(1,windowSize);
a = 1;
disp(['Applying filter of width ', num2str(timeWindow)]);
smoothdata = filter(b,a,data,[],2);
time = time(ind);
outdata = [smoothdata(:,ind);time];
%}}}
