function md=readelmermesh(domainfile)

	%Read node file
	filename = [domainfile '.nodes'];
	disp(['Reading ' filename]);
	nodes = load(filename);

	%Get coordinates
	x = nodes(:,3);
	y = nodes(:,4);

	%Read element file
	filename = [domainfile '.elements'];
	disp(['Reading ' filename]);
	elements = load(filename);

	%Get indices
	index = elements(:,4:6);

	%Convert mesh
	md=meshconvert(model,index,x,y);
