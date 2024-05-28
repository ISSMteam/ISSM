function pdf=lognormal_pdf(x,mu,sigma)

pdf=1./(x*sigma*sqrt(2*pi)) .* exp( -(log(x)-mu).^2/(2*sigma^2)  );

