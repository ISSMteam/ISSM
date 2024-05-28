function pdf=gaussian_pdf(x,average,sigma)

pdf=1/(sqrt(2*pi*sigma^2))*exp(-(x-average).^2/(2*sigma^2));
