function md=meshyamsrecreateriftsegments(md)

	%recreate rift segments: just used for yams. temporaroy routine.
	pos_record=[];
	if md.rifts.numrifts,
		for i=1:md.rifts.numrifts,
			rift=md.rifts.riftstruct(i);

			%closed rifts first:
			if length(rift.tips)==2,

				%find tip1 and tip2 for this rift, in the new mesh created by yams.
				pos=find_point(md.mesh.x(md.mesh.segments(:,1)),md.mesh.y(md.mesh.segments(:,1)),rift.tip1coordinates(1),rift.tip1coordinates(2));
				tip1=md.mesh.segments(pos,1);
				pos=find_point(md.mesh.x(md.mesh.segments(:,1)),md.mesh.y(md.mesh.segments(:,1)),rift.tip2coordinates(1),rift.tip2coordinates(2));
				tip2=md.mesh.segments(pos,1);

				%start from tip1, and build segments of this rift. 
				pos=find_point(md.mesh.x(md.mesh.segments(:,1)),md.mesh.y(md.mesh.segments(:,1)),rift.tip1coordinates(1),rift.tip1coordinates(2));
				pos_record=[pos_record; pos];
				riftsegs=md.mesh.segments(pos,:);
				while 1,
					A=riftsegs(end,1); B=riftsegs(end,2); el=riftsegs(end,3);
					%find other segment that holds B.
					pos=find(md.mesh.segments(:,1)==B);
					pos_record=[pos_record; pos];
					riftsegs=[riftsegs; md.mesh.segments(pos,:)];
					if riftsegs(end,2)==tip1, 
						break;
					end
				end
				md.rifts.riftstruct(i).segments=riftsegs;
				md.rifts.riftstruct(i).tips=[tip1 tip2];

			else
				%ok, this is a rift that opens up to the domain outline.  One tip is going to be 
				%double, the other one, single. We are going to start from the single tip, towards the two 
				%other doubles

				%find tip1 and tip2 for this rift, in the new mesh created by yams.
				pos1=find_point(md.mesh.x(md.mesh.segments(:,1)),md.mesh.y(md.mesh.segments(:,1)),rift.tip1coordinates(1),rift.tip1coordinates(2));
				tip1=md.mesh.segments(pos1,1);
				pos2=find_point(md.mesh.x(md.mesh.segments(:,1)),md.mesh.y(md.mesh.segments(:,1)),rift.tip2coordinates(1),rift.tip2coordinates(2));
				tip2=md.mesh.segments(pos2,1);
				if length(tip1)==2,
					%swap.
					temp=tip1; tip1=tip2; tip2=temp;
					temp=pos1; pos1=pos2; pos2=temp;
					pos=pos1;
				else
					pos=pos1;
				end

				pos_record=[pos_record; pos];
				riftsegs=md.mesh.segments(pos,:);
				while 1,
					A=riftsegs(end,1); B=riftsegs(end,2); el=riftsegs(end,3);
					%find other segment that holds B.
					pos=find(md.mesh.segments(:,1)==B);
					pos_record=[pos_record; pos];
					riftsegs=[riftsegs; md.mesh.segments(pos,:)];
					if ((riftsegs(end,2)==tip2(1)) | (riftsegs(end,2)==tip2(2))), 
						%figure out which tip we reached
						if riftsegs(end,2)==tip2(1), index=2; else index=1; end
						break;
					end
				end

				%ok, now, we start from the other tip2, towards tip1
				pos=pos2(index);
				pos_record=[pos_record; pos];
				riftsegs=[riftsegs; md.mesh.segments(pos,:)];
				while 1,
					A=riftsegs(end,1); B=riftsegs(end,2); el=riftsegs(end,3);
					%find other segment that holds B.
					pos=find(md.mesh.segments(:,1)==B);
					pos_record=[pos_record; pos];
					riftsegs=[riftsegs; md.mesh.segments(pos,:)];
					if riftsegs(end,2)==tip1, 
						break;
					end
				end
				md.rifts.riftstruct(i).segments=riftsegs;
				md.rifts.riftstruct(i).tips=[tip1 tip2(1) tip2(2)];

			end
		end
	end
	%take out rift segments from segments
	md.mesh.segments(pos_record,:)=[];
