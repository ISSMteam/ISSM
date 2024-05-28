%
%  plot a normal probability plot of the response functions.
%
%  []=plot_normplot(rfunc)
%
function []=plot_qmunormplot(rfunc,width,ii)

if ~nargin
    help plot_normplot
    return
end

%%  assemble the data into a matrix

desc=cell (1,length(rfunc));
for i=1:length(rfunc)
    ldata(i)=length(rfunc(i).sample);
end
data=zeros(max(ldata),length(rfunc));

for i=1:length(rfunc)
    desc(i)=cellstr(rfunc(i).descriptor);
    data(1:ldata(i),i)=rfunc(i).sample;
end

%standard plot:
subplot(width,width,ii);

%%  draw the plot

%  draw normal probability plot

normplot(data)
ax1=gca;

%  add the annotation

title('Normal Probability Plot of Design Variables and/or Response Functions')
xlabel('Value')
ylabel('Probability')

hleg1=legend(ax1,desc,'Location','EastOutside',...
             'Orientation','vertical','Interpreter','none');

end
