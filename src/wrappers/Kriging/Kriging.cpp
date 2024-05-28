/*\file Kriging.c
 *\brief: best linear predictor
 */
#include "./Kriging.h"

void KrigingUsage(void){/*{{{*/
	int num=1;
#ifdef _MULTITHREADING_
	num=_NUMTHREADS_;
#endif
	_printf0_("\n");
	_printf0_("   usage: predictions=" << __FUNCT__ << "(x,y,observations,x_interp,y_interp,'options');\n");
	_printf0_("   available options:\n");
	_printf0_("      -'model': Available variogram models 'gaussian' (default),'spherical','power','exponential'\n");
	_printf0_("         -'nugget': nugget effect (default 0.2)\n");
	_printf0_("         -'range':  for gaussian, spherical and exponential models (default sqrt(3))\n");
	_printf0_("         -'sill':   for gaussian, spherical and exponential models (default 1)\n");
	_printf0_("         -'slope':  for power model (default 1)\n");
	_printf0_("         -'power':  for power model (default 1)\n");
	_printf0_("      -'searchradius': search radius for each prediction (default is observations span)\n");
	_printf0_("      -'boxlength':    minimum length of quadtree boxes (useful to decrease the number of observations)\n");
	_printf0_("      -'maxdata':      minimum number of observations for a prediction (default is 50)\n");
	_printf0_("      -'mindata':      maximum number of observations for a prediction (default is 1)\n");
	_printf0_("      -'maxtrimming':  maximum trimming value (default is -1.e+21)\n");
	_printf0_("      -'mintrimming':  minimum trimming value (default is +1.e+21)\n");
	_printf0_("      -'minspacing':   minimum distance between observation (default is 0.01)\n");
	_printf0_("      -'numthreads':   number of threads, default is "<<num << "\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(Kriging_python){

	/*Outputs*/
	double  *x            = NULL;
	double  *y            = NULL;
	double  *observations = NULL;
	double  *x_interp     = NULL;
	double  *y_interp     = NULL;
	double  *predictions  = NULL;
	double  *error        = NULL;
	Options *options      = NULL;
	int      N_interp,M_interp,M,N,n_obs;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	if (nrhs<NRHS || nlhs>NLHS){
		KrigingUsage(); _error_("Kriging usage error");
	}

	/*Fetch inputs: */
	FetchData(&x,&n_obs,X);
	FetchData(&y,&N,Y);                       if(n_obs!=N) _error_("x and y should have the same size");
	FetchData(&observations,&N,OBSERVATIONS); if(n_obs!=N) _error_("x and observations should have the same size");
	FetchData(&x_interp,&M_interp,&N_interp,XINTERP);
	FetchData(&y_interp,&M,&N,YINTERP);       if(N!=N_interp || M!=M_interp) _error_("x_interp and y_interp should have the same size");
	FetchData(&options,NRHS,nrhs,prhs);

	/*Call x layer*/
	Krigingx(&predictions,&error,x,y,observations,n_obs,x_interp,y_interp,M_interp*N_interp,options);

	/*Generate output Matlab Structures*/
	if(nlhs>=1) WriteData(PREDICTIONS,predictions,M_interp,N_interp);
	if(nlhs==2) WriteData(ERROR,error,M_interp,N_interp);

	/*Free resources: */
	xDelete<double>(x);
	xDelete<double>(y);
	xDelete<double>(observations);
	xDelete<double>(x_interp);
	xDelete<double>(y_interp);
	xDelete<double>(predictions);
	xDelete<double>(error);
	delete options;

	/*end module: */
	MODULEEND();
}
