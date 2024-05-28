function lsf=exp2levelsetfunction(md, exp_icedomain)
	%EXP2LEVELSETFUNCTION: returns signed distance function from EXP-file
	%
	%   This routine computes a signed distance function from an EXP-file given in the input.
	%   It can be used with the level-set method.
	%
	%   USAGE:
	%      levelsetfunction=exp2levelsetfunction(md, exp_icedomain)
	%

	mesh=md.mesh;
	profiles=expread(exp_icedomain);

	min_dist=NaN(size(mesh.x));
	for p=1:size(profiles,2)
		profile=profiles(p);

		%construct ice domain segments
		inds_v=1:profile.nods-1;
		inds_w=2:profile.nods;
		segments.v.x=profile.x(inds_v);	segments.v.y=profile.y(inds_v);
		segments.w.x=profile.x(inds_w);	segments.w.y=profile.y(inds_w);
		segments.numsegments=length(segments.v.x);

		% compute minimum distance to segments
		for s=1:segments.numsegments
			segment.v.x=segments.v.x(s);	segment.v.y=segments.v.y(s);
			segment.w.x=segments.w.x(s);	segment.w.y=segments.w.y(s);
			min_dist=min(min_dist, compute_distance_to_segment(segment, mesh.x, mesh.y));
		end
	end

	% set sign of lsf
	sign_lsf=ones(mesh.numberofvertices,1);
	isice=ContourToMesh(mesh.elements,mesh.x,mesh.y,exp_icedomain,'node',2);
	sign_lsf(find(isice))=-1;

	lsf=sign_lsf.*min_dist;

	function dist=compute_distance_to_segment(segment,x,y)
		%compute horizontal euclidean distance to segment
		v=[segment.v.x segment.v.y];
		w=[segment.w.x segment.w.y];
		verts=[x y];
		dist_vw2=norm(w-v)^2;
		if(dist_vw2==0.),	t=zeros(size(x)); %cover case where segment has length 0
		else t=[x-v(1) y-v(2)]*(w-v)'/dist_vw2; end %projection of verts on line defined by v and w
		dist_vec=(ones(length(x),1)*v+max(0,min(1,t))*(w-v))-verts; %vector of shortest distance between verts and segment v-w
		dist=sqrt(sum(abs(dist_vec).^2,2));
	end
end
