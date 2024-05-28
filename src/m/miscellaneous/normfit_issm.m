%
%  wrapper for normfit to avoid using the matlab statistics toolbox.
%
function [muhat,sigmahat,muci,sigmaci]=normfit_issm(x,alpha)

	if ~exist('alpha','var')
		alpha=0.05;
	end

%  check for any NaN in any columns

	if ~any(any((isnan(x))))

%  explicitly calculate the moments

		muhat   =mean(x);
		sigmahat=std(x);

		if (nargout>2)
			prob=1.-alpha/2.;

			if (size(x,1) == 1)
				% operate like matlab normfit, mean, std, etc.
				n=length(x);
			else
				n=size(x,1);
			end

			muci    =zeros(2,length(muhat   ));
			sigmaci =zeros(2,length(sigmahat));

			try
				muci(1,:)   =muhat-tinv(prob,n-1)*sigmahat/sqrt(n);
				muci(2,:)   =muhat+tinv(prob,n-1)*sigmahat/sqrt(n);
				sigmaci(1,:)=sigmahat*sqrt((n-1)/chi2inv(prob   ,n-1));
				sigmaci(2,:)=sigmahat*sqrt((n-1)/chi2inv(1.-prob,n-1));
			catch me
				muci(1,:)   =muhat;
				muci(2,:)   =muhat;
				sigmaci(1,:)=sigmahat;
				sigmaci(2,:)=sigmahat;
			end
		end

	else

%  must loop over columns, since number of elements could be different

		muhat   =zeros(1,size(x,2));
		sigmahat=zeros(1,size(x,2));
		muci    =zeros(2,size(x,2));
		sigmaci =zeros(2,size(x,2));

%  remove any NaN and recursively call column

		for j=1:size(x,2)
			[muhat(j),sigmahat(j),muci(:,j),sigmaci(:,j)]=normfit_issm(x(~isnan(x(:,j)),j),alpha);
		end
	end
end
