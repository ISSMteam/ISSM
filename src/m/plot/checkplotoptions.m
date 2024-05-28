function options=checkplotoptions(md,options)
%PARSE_OPTIONS - build a structure that holds all plot options
%
%   Usage:
%      options=checkplotoptions(md,options);
%
%   See also: PLOTMODEL

%units
if exist(options,'unit'),
	if strcmpi(getfieldvalue(options,'unit'),'km')
		options=changefieldvalue(options,'unit',10^-3);
	end
	if strcmpi(getfieldvalue(options,'unit'),'100km')
		options=changefieldvalue(options,'unit',10^-5);
	end

end

%density
if exist(options,'density'),
	density=getfieldvalue(options,'density');
	options=changefieldvalue(options,'density',abs(ceil(density)));
end

%Show section
if exist(options,'showsection'),
	if strcmpi(getfieldvalue(options,'showsection'),'on')
		options=changefieldvalue(options,'showsection',4);
	end
end

%smooth values
if exist(options,'smooth'),
	if strcmpi(getfieldvalue(options,'smooth'),'on')
		options=changefieldvalue(options,'smooth',0);
	end
end

%contouronly values
if exist(options,'contouronly'),
	if strcmpi(getfieldvalue(options,'contouronly'),'on')
		options=changefieldvalue(options,'contouronly',1);
	end
end

%Colorbar;
if exist(options,'colorbar'),
	if strcmpi(getfieldvalue(options,'colorbar'),'on')
		options=changefieldvalue(options,'colorbar',1);
	elseif strcmpi(getfieldvalue(options,'colorbar'),'off')
			options=changefieldvalue(options,'colorbar',0);
	end
end

%text
if exist(options,'text'),
	%1: textvalue
	textvalues=getfieldvalue(options,'text');
	%ischar if only one expstyle -> create a cell
	if ischar(textvalues),
		textvalues={textvalues};
		numtext=1;
	elseif iscell(textvalues),
		numtext=length(textvalues);
	else
		error('plot error message: ''text'' option should be either a string or a cell');
	end

	%2: textweight
	if exist(options,'textweight'),
		textweightvalues=getfieldvalue(options,'textweight');
		%ischar if only one textweight -> create a cell
		if ischar(textweightvalues),
			textweightvalues={textweightvalues};
		elseif ~iscell(textweightvalues);
			error('plot error message: ''textweight'' option should be either a string or a cell');
		end
	else
		textweightvalues={'n'};
	end
	textweightvalues=repmat(textweightvalues,1,numtext); textweightvalues(numtext+1:end)=[];

	%3: textsize
	if exist(options,'textsize'),
		textsizevalues=getfieldvalue(options,'textsize');
		%ischar if only one textsize -> create a cell
		if isnumeric(textsizevalues),
			textsizevalues={textsizevalues};
		elseif ~iscell(textsizevalues);
			error('plot error message: ''textsize'' option should be either a number or a cell');
		end
	else
		textsizevalues={14};
	end
	textsizevalues=repmat(textsizevalues,1,numtext); textsizevalues(numtext+1:end)=[];
	%4: textcolor
	if exist(options,'textcolor'),
		textcolorvalues=getfieldvalue(options,'textcolor');
		%ischar if only one textcolor -> create a cell
		if ischar(textcolorvalues),
			textcolorvalues={textcolorvalues};
		elseif ~iscell(textcolorvalues);
			error('plot error message: ''textcolor'' option should be either a string or a cell');
		end
	else
		textcolorvalues={'k'};
	end
	textcolorvalues=repmat(textcolorvalues,1,numtext); textcolorvalues(numtext+1:end)=[];
	%5: textposition
	if exist(options,'textposition'),
		textpositionvalues=getfieldvalue(options,'textposition');
		%ischar if only one textposition -> create a cell
		if isnumeric(textpositionvalues),
			textpositionvalues={textpositionvalues};
		elseif ~iscell(textpositionvalues);
			error('plot error message: ''textposition'' option should be either a string or a cell');
		end
	else
		error('plot error message: ''textposition'' option is missing');
	end
	%6: textrotation
	if exist(options,'textrotation'),
		textrotationvalues=getfieldvalue(options,'textrotation');
		%ischar if only one textsize -> create a cell
		if isnumeric(textrotationvalues),
			textrotationvalues={textrotationvalues};
		elseif ~iscell(textrotationvalues);
			error('plot error message: ''textrotation'' option should be either a number or a cell');
		end
	else
		textrotationvalues={0};
	end
	textrotationvalues=repmat(textrotationvalues,1,numtext); textrotationvalues(numtext+1:end)=[];
	options=changefieldvalue(options,'text',textvalues);
	options=changefieldvalue(options,'textsize',textsizevalues);
	options=changefieldvalue(options,'textweight',textweightvalues);
	options=changefieldvalue(options,'textcolor',textcolorvalues);
	options=changefieldvalue(options,'textposition',textpositionvalues);
	options=changefieldvalue(options,'textrotation',textrotationvalues);
end

%expdisp
expdispvaluesarray=cell(0,0);
expstylevaluesarray=cell(0,0);
expstylevalues=cell(0,0);
if exist(options,'expstyle'),
	expstylevalues=getfieldvalue(options,'expstyle');
	%ischar if only one expstyle -> create a cell
	if ischar(expstylevalues),
		expstylevalues={expstylevalues};
	end
end
if exist(options,'expdisp'),
	expdispvalues=getfieldvalue(options,'expdisp');
	%ischar if only one expstyle -> create a cell
	if ischar(expdispvalues),
		expdispvalues={expdispvalues};
	end
	for i=1:length(expdispvalues)
		expdispvaluesarray{end+1}=expdispvalues{i};
		if (length(expstylevalues)>=i),
			expstylevaluesarray{end+1}=expstylevalues{i};
		else
			expstylevaluesarray{end+1}='g-';
		end
	end
end
options=changefieldvalue(options,'expstyle',expstylevaluesarray);
options=changefieldvalue(options,'expdisp',expdispvaluesarray);

%latlonnumbering
if exist(options,'latlonclick'),
	if strcmpi(getfieldvalue(options,'latlonclick'),'on')
		options=changefieldvalue(options,'latlonclick',1);
	end
end

%north arrow
if exist(options,'northarrow'),
	if strcmpi(getfieldvalue(options,'northarrow'),'on')
		%default values
		Lx=max(md.mesh.y)-min(md.mesh.y);
		Ly=max(md.mesh.y)-min(md.mesh.y);
		%default values
		options=changefieldvalue(options,'northarrow',[min(md.mesh.x)+1/6*Lx   min(md.mesh.y)+5/6*Ly   1/15*Ly   0.25   1/250*Ly]);
	end
end

%scale ruler
if exist(options,'scaleruler'),
	if strcmpi(getfieldvalue(options,'scaleruler'),'on')
		%default values
		Lx=max(md.mesh.x)-min(md.mesh.x);
		Ly=max(md.mesh.y)-min(md.mesh.y);
		%default values
		options=changefieldvalue(options,'scaleruler',[min(md.mesh.x)+6/8*Lx   min(md.mesh.y)+1/10*Ly   10^(ceil(log10(Lx)))/5 floor(Lx/100) 5]);
	end
end

%Log scale (LOTS of changes to be performed
if exist(options,'log'),
	if exist(options,'caxis')
		options=addfield(options,'caxis_pre',getfieldvalue(options,'caxis'));
		options=changefieldvalue(options,'caxis',log(getfieldvalue(options,'caxis'))/log(getfieldvalue(options,'log')));
	end
	options=changefieldvalue(options,'cutoff',log(getfieldvalue(options,'cutoff',1.5))/log(getfieldvalue(options,'log')));
end
