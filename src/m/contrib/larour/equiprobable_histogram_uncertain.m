function [abscissas,counts,pair_per_variable]=equiprobable_histogram_uncertain(numberabscissas)

	pair_per_variable=numberabscissas+1;
	abscissas=(1:(numberabscissas+1))';
	p=ones(numberabscissas,1)/numberabscissas;
	counts=[p;0]; counts=counts/sum(counts);
