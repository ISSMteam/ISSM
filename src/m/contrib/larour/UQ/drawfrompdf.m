function s=drawfrompdf(x,y,nsamples,mode)
% Samples a 1-D pdf histogram
% x: variable to be sampled
% y: pdf histogram of x, best if y(1)=0 in continuous sample mode
% nsamples: how many samples must be generated
% mode: 'continuous' or 'discrete'. 'discrete' will generate only values provided in x, 'continuous' will interpolate such that values provided be < x(end), or in [x(1) x(end)] if y(1)=0 
 y=[cumsum(y/sum(y))]; %compute the CDF
 s=zeros(nsamples,1);

if y(1)~=0 & strcmpi(mode,'continuous')
    disp('WARNING: y(1)~=0, values < x(1) might be generated extrapolating y(2)-y(1)')
end

 for i=1:nsamples
     r=rand(1); %uniform [0 1] sampling
    j=2;    
    while r>y(j)
        j=j+1;
    end
    if strcmpi(mode,'continuous')
    d=(r-y(j-1))/(y(j)-y(j-1));
    s(i)=x(j-1)*(1.d0-d)+x(j)*d; %linear interpolation of the cdf to map r into x
    elseif strcmpi(mode, 'discrete')
        s(i)=x(j);
    else 
        error(' ''mode'' must be set to ''discrete'' or ''continuous'' ')
    end
 end
 
