function plot_qmu_mass_flux_segments(md,options,nlines,ncols,i)
%PLOT_QMU_MASS_FLUX_SEGMENTS - plot segments from the qmu analysis of mass fluxes
%
%   Usage:
%      plot_qmu_mass_flux_segments(md,options,nlines,ncols,i);
%

subplot(nlines,ncols,i); 

%process mesh and data
[x y z elements is2d isplanet]=processmesh(md,[],options);

allsegments=md.qmu.mass_flux_segments;

if dimension(md.mesh)==2,

	%recover segments
	hold on
	for i=1:length(allsegments),
		segments=allsegments{i};

		%plot semgnets
		for j=1:length(segments),
			plot([segments(j,1) segments(j,3)],[segments(j,2) segments(j,4)]);
		end
		text(segments(j,1),segments(j,2),['Profile #' num2str(i)]);

		%plot normals

		for j=1:length(segments),
			xstart=mean([segments(j,1) segments(j,3)]);
			ystart=mean([segments(j,2) segments(j,4)]);
			length1=sqrt((segments(j,1)-segments(j,3)).^2 + (segments(j,2)-segments(j,4)).^2);
			normal(:,1)=cos(atan2(segments(j,1)-segments(j,3) , segments(j,4)-segments(j,2)));
			normal(:,2)=sin(atan2(segments(j,1)-segments(j,3) , segments(j,4)-segments(j,2)));
			xend=xstart+length1.*normal(:,1);
			yend=ystart+length1.*normal(:,2);
			plot([xstart xend],[ystart yend],'r-');
		end

	end
else
	error('plot_qmu_mass_flux_segments: 3d plot of segments not supported yet!');
end

%apply options
options=addfielddefault(options,'title','Mass Flux segments and normals');
options=addfielddefault(options,'colorbar',0);
applyoptions(md,[],options);
