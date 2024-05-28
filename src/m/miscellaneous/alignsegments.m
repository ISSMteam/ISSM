function newsegments=alignsegments(segments)
%ALIGNSEGMENTS: 
% 
%

	nt=length(segments);
	newsegments=zeros(nt,3);
	newsegments(1,:)=segments(1,:);

	for  i=2:nt, 
		last=newsegments(i-1,2); %last vertex of the previous segment: 
		for j=1:nt,
			if last==segments(j,1),
				%we found the next segment: 
				newsegments(i,:)=segments(j,:);
				break;
			end
		end
	end
