function plot_qmuhistnorm(md,options,nlines,ncols,index)

%prepare plot
subplot(nlines,ncols,index); 
hold on

%recover histnorm data
if ~exist(options,'qmudata')
	error('plot_qmuhistnorm error message: option qmudata is required');
else
	qmudata=getfieldvalue(options,'qmudata');
end

%process options for the qmu plot: 

%    hmin          (numeric, minimum for histogram)
%    hmax          (numeric, maximum for histogram)
%    hnint         (numeric, number of intervals for histogram)
%    ymin1         (numeric, minimum of histogram y-axis)
%    ymax1         (numeric, maximum of histogram y-axis)
%    ymin2         (numeric, minimum of cdf y-axis)
%    ymax2         (numeric, maximum of cdf y-axis)
%    cdfplt        (char, 'off' to turn off cdf line plots)
%    cdfleg        (char, 'off' to turn off cdf legends)
%

qmuoptions='';

if exist(options,'hmin'), hmin=getfieldvalue(options,'hmin'); qmuoptions=[qmuoptions ',''hmin'',' num2str(hmin)]; end
if exist(options,'hmax'), hmax=getfieldvalue(options,'hmax'); qmuoptions=[qmuoptions ',''hmax'',' num2str(hmax)]; end
if exist(options,'hnint'), hnint=getfieldvalue(options,'hnint'); qmuoptions=[qmuoptions ',''hnint'',' num2str(hnint)]; end
if exist(options,'ymin1'), ymin1=getfieldvalue(options,'ymin1'); qmuoptions=[qmuoptions ',''ymin1'',' num2str(ymin1)]; end
if exist(options,'ymax1'), ymax1=getfieldvalue(options,'ymax1'); qmuoptions=[qmuoptions ',''ymax1'',' num2str(ymax1)]; end
if exist(options,'ymin2'), ymin2=getfieldvalue(options,'ymin2'); qmuoptions=[qmuoptions ',''ymin2'',' num2str(ymin2)]; end
if exist(options,'ymax2'), ymax2=getfieldvalue(options,'ymax2'); qmuoptions=[qmuoptions ',''ymax2'',' num2str(ymax2)]; end
if exist(options,'cdfplt'), cdfplt=getfieldvalue(options,'cdfplt'); qmuoptions=[qmuoptions ',''cdfplt'',''' cdfplt '''']; end
if exist(options,'cdfleg'), cdfleg=getfieldvalue(options,'cdfleg'); qmuoptions=[qmuoptions ',''cdfleg'',''' cdfleg '''']; end
if exist(options,'nrmplt'), nrmplt=getfieldvalue(options,'nrmplt'); qmuoptions=[qmuoptions ',''nrmplt'',''' nrmplt '''']; end
if exist(options,'EdgeColor'), EdgeColor=getfieldvalue(options,'EdgeColor'); qmuoptions=[qmuoptions ',''EdgeColor'',''' EdgeColor '''']; end
if exist(options,'FaceColor'), FaceColor=getfieldvalue(options,'FaceColor'); qmuoptions=[qmuoptions ',''FaceColor'',''' FaceColor '''']; end

%use qmu plot
eval(['plot_hist_norm(qmudata' qmuoptions ');']);

%apply options
options=changefieldvalue(options,'colorbar','off');
applyoptions(md,[],options);
