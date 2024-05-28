function nodefield=ecco32issm(field,transition,xecco3,yecco3)

	xecco3linear=xecco3(:); yecco3linear=yecco3(:); %linearize
	nodefieldlinear=zeros(length(xecco3linear),1);
	nodefieldlinear(transition(:,1))=field(transition(:,2));
	nodefield=xecco3;
	nodefield(:)=nodefieldlinear;
	%nodefield=nodefield'; %not sure we need that
