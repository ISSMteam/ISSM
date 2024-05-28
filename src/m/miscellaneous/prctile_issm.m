%
%  wrapper for prctile to avoid using the matlab statistics toolbox.
%
function [y]=prctile_issm(x,p,dim)

	try
		y=prctile(argin{:});

	catch me
		if length(size(x)) > 2
			error('Number of dimensions %d not implemented.',length(size(x)));
		end
		if ~exist('dim','var')
			dim=0;
			for i=1:length(size(x))
				if ~dim && size(x,i)>1
					dim=i;
				end
			end
			if ~dim
				dim=1;
			end
		end

		psize=size(p);
		if size(p,2)>1
			p=transpose(p);
		end

		xsize=size(x);
		if dim==2
			x=transpose(x);
		end

%  check for any NaN in any columns

		if ~any(any((isnan(x))))
			x=sort(x,1);
			n=size(x,1);

%  branch based on number of elements

			if     n>1

%  set up percent values and interpolate

				xi=transpose(100.*([1:n]-0.5)/n);
				y=interp1q(xi,x,p);

%  fill in high and low values
				y(p<xi(1),:)=repmat(x(1,:),nnz(p<xi(1)),1);
				y(p>xi(n),:)=repmat(x(n,:),nnz(p>xi(n)),1);

%  if one value, just copy it

			elseif n==1
				y=repmat(x(1,:),length(p),1);

%  if no values, use NaN

			else
				y=repmat(NaN,size(p,1),size(x,2));
			end

		else

%  must loop over columns, since number of elements could be different

			y=zeros(size(p,1),size(x,2));
			for j=1:size(x,2)

%  remove any NaN and recursively call column

				y(:,j)=prctile_issm(x(~isnan(x(:,j)),j),p);
			end
		end

		if (min(xsize)==1 && xsize(dim)>1 && psize(2)>1) || ...
		   (min(xsize)> 1 && dim==2)
			y=transpose(y);
		end
	end
end
