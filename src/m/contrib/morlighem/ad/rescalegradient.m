function grad_out = rescalegradient(md,grad_in);
%RESCALEGRADIENT - rescale gradient using mass matrix
%
%   Usage:
%      grad_out = rescalegradient(md,grad_in);

   %Define index
	index = md.mesh.elements;

	%Get surface areas of all elements
	A = GetAreas(index,md.mesh.x,md.mesh.y);

	%Preallocate to speed up computation
	disp('Constructing mass matrix...');
	tic
	row   = zeros(10*md.mesh.numberofvertices,1);
	col   = zeros(10*md.mesh.numberofvertices,1);
	value = zeros(10*md.mesh.numberofvertices,1);

	%Construct mass matrix using MATLAB's sparse function
	count = 0;
	for n=1:md.mesh.numberofelements
		for l=1:3
			for k=1:3
				count=count+1;
				row(count) = index(n,k);
				col(count) = index(n,l);
				if l == k
					value(count) = A(n)/6.;  % \int_E phi_i * phi_i dE = A/6
				else
					value(count) = A(n)/12.; % \int_E phi_i * phi_i dE = A/12
				end 
			end
		end
	end

	%Delete unused elements
	row = row(1:count);
	col = col(1:count);
	value = value(1:count);

	%Make mass matrix
	M=sparse(row,col,value);
	toc

	tic
	disp('Solving...');
	grad_out = M\grad_in;
	toc

	disp('Adjusting output');
	pos = find(grad_in==0);
	grad_out(pos)==0;
