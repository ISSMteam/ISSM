function plot_transient_movie(md,options,width,i)
%PLOT_TRANSIENT_MOVIE - plot a transient result as a movie
%   Usage:
%      plot_transient_movie(md,options,width,i);
%
%   See also: PLOTMODEL, PLOT_UNIT, PLOT_MANAGER

	%prepare subplot
	subplot(width,width,i);

	%xlim
	if exist(options,'transient_movie_field'),
		field=getfieldvalue(options,'transient_movie_field');
	elseif ischar(getfieldvalue(options,'data')) && ~strcmp(getfieldvalue(options,'data'),'transient_movie')
		field=getfieldvalue(options,'data');
	else
		disp('List of available fields:');
		F=fieldnames(md.results.TransientSolution(1));
		num = [];
		for i=1:numel(F),
			if ~strcmp(F{i},'time') & ...
				~strcmp(F{i},'step') & ...
				~strcmp(F{i},'errlog') & ...
				~strcmp(F{i},'outlog') & ...
				~strcmp(F{i},'MaxIterationConvergenceFlag') & ...
				~strcmp(F{i},'SolutionType'),
				disp(['   ' num2str(i) ': ' F{i} ]);
				num = [num i];
			end
		end
		choice=input(['please enter the field number? (between ' num2str(min(num)) ' and ' num2str(max(num)) ')  ']);
		field =  F{choice};
	end

	results=md.results.TransientSolution;
	%loop over the time steps
	if exist(options,'transient_movie_limit'),
		limit=getfieldvalue(options,'transient_movie_limit');
		steps=[limit(1):limit(end)];
	elseif exist(options,'transient_movie_steps'),
		warning('option ''transient_movie_steps'' is now ''steps'', please update your script');
		steps = getfieldvalue(options,'transient_movie_steps');
	elseif exist(options,'steps'),
		steps = getfieldvalue(options,'steps');
	elseif exist(options,'step'),
		steps = getfieldvalue(options,'step');
	else
		steps=1:length(results);
	end

	%Do we have an output?
	isavi = 0;
	isgif = 0;
	ismp4 = 0;
	if exist(options,'transient_movie_output'),
		filename=getfieldvalue(options,'transient_movie_output');
		[pathstr,name,ext] = fileparts(filename);
		if strcmp(ext,'.gif')
			isgif = 1;
		elseif strcmp(ext,'.mp4')
			ismp4 = 1;
		elseif strcmp(ext,'.avi')
			isavi = 1;
		end
	end
	if isavi || ismp4
		vid=VideoWriter([filename(1:end-4),'.avi'],'Motion JPEG AVI');
		vid.FrameRate = 2; 
		open(vid); 
	end

	%calculate caxis
	if ~exist(options,'caxis'),
		range = [Inf -Inf];
		for i=steps
			if isfield(results(i), 'MeshElements')
				options=changefieldvalue(options,'amr', i);
			end
			[data datatype]=processdata(md,results(i).(field), options);
			range(1) = min(range(1),min(data));
			range(2) = max(range(2),max(data));
		end
		options=addfielddefault(options,'caxis',range);
	end

	%Process mesh once for all
	[x y z elements is2d isplanet]=processmesh(md,results(i).(field),options);

	%display movie
	nstep=1;
	deltat = getfieldvalue(options,'pause',.5);
	for i=steps

		if ~isempty(results(i).(field)),
			%Process mesh if necessary
			if isfield(results(i), 'MeshElements')
				options=changefieldvalue(options,'amr', i);
				[x y z elements is2d isplanet]=processmesh(md,results(i).(field),options);
			end

			%process data
			[data datatype]=processdata(md,results(i).(field),options);

			clf;
			titlestring=[field ' (time ' num2str(results(i).time,'%7.2f') ' yr)'];
			plot_unit(x,y,z,elements,data,is2d,isplanet,datatype,options)
			apply_options_movie(md,data,options,titlestring);

			%Add grounding line
			if exist(options,'groundingline') 
				if dimension(md.mesh)==2
					contours=isoline(md, results(i).MaskOceanLevelset,'output','matrix');
				else
					ocean = project2d(md, results(i).MaskOceanLevelset, 1);
					contours=isoline(md, ocean,'output','matrix');
				end
				hold on
				plot(contours(:,1),contours(:,2),getfieldvalue(options,'groundingline'));
			end

			%Add ice front
			if exist(options,'icefront')
				if dimension(md.mesh)==2
					if exist(options, 'amr')
						contours=isoline(md, results(i).MaskIceLevelset,'output','matrix', 'amr', results(i));
					else
						contours=isoline(md, results(i).MaskIceLevelset,'output','matrix');
					end
				else
					ice = project2d(md, results(i).MaskIceLevelset, 1);
					contours=isoline(md, ice,'output','matrix');
				end
				hold on
				plot(contours(:,1),contours(:,2),getfieldvalue(options,'icefront'));
			end

			if isgif
				frame=getframe(gcf);
				im = frame2im(frame);
				[imind,cmap] = rgb2ind(im,256);
				if i==1
					imwrite(imind, cmap, filename, 'DelayTime',getfieldvalue(options,'transient_movie_time',.5), 'LoopCount',inf)
				else
					imwrite(imind, cmap, filename, 'WriteMode','append');
				end
			elseif isavi || ismp4
				F=getframe(gcf);
				writeVideo(vid, F);
			end

			if exist(options,'transient_movie_output')
				%set(gcf,'Renderer','zbuffer','color','white'); %fixes a bug on Mac OS X (not needed in future Matlab version)
				if nstep==1,
					%initialize images and frame
					frame=getframe(gcf);
					[images,map]=rgb2ind(frame.cdata,256,'nodither');
					images(1,1,1,length(steps))=0;
				else
					frame=getframe(gcf);
					images(:,:,1,nstep) = rgb2ind(frame.cdata,map,'nodither');
				end
			else
				pause(deltat)
			end
			nstep=nstep+1;
		end
	end

	%output movie if requested.
	if isavi || ismp4
		close(vid);
	end
	if ismp4
		filename = filename(1:end-4);
		while(~exist([filename '.avi']))
			disp(['Waiting for ' filename '.avi']);
			pause(1)
		end
		command=sprintf('ffmpeg -y -i %s.avi -c:v libx264 -crf 19 -preset slow -c:a libfaac -b:a 192k -ac 2 -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" %s.mp4',filename,filename);
		system(command);
	end
	if isgif
		imwrite(images,map,filename,'DelayTime',getfieldvalue(options,'transient_movie_time',.5),'LoopCount',inf)
	end

end %function

function apply_options_movie(md,data,options,titlestring)
	%apply options
	options=changefieldvalue(options,'title',titlestring);
	options=addfielddefault(options,'colorbar',1);
	applyoptions(md,data,options);
end
