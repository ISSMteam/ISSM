/*!\file M1qn3.c
 * \brief: data interpolation from a list of (x,y,values) into mesh vertices
*/

#include "./M1qn3.h"

#ifdef _HAVE_M1QN3_
/*m1qn3 prototypes {{{*/
extern "C" void *ctonbe_; // DIS mode : Conversion
extern "C" void *ctcabe_; // DIS mode : Conversion
extern "C" void *euclid_; // Scalar product
typedef void (*SimulFunc) (long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs);
extern "C" void m1qn3_ (void f(long* indic,long* n, double* x, double* pf,double* g,long [],float [],void* dzs),
			void **, void **, void **,
			long *, double [], double *, double [], double*, double *,
			double *, char [], long *, long *, long *, long *, long *, long *, long [], double [], long *,
			long *, long *, long [], float [],void* );
/*Cost function prototype*/
void fakesimul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs);
typedef struct {
	int priorn; 
	int counter; 
	double* Gs; 
	double* Xs;
	double* Js;
} Data; 
/*}}}*/
void M1qn3Usage(void){/*{{{*/
	_printf0_("   usage:\n");
	_printf0_("         X=M1qn3(Xs,Gs);\n");
	_printf0_("   where:\n");
	_printf0_("      Xs are the X values (m x n, where m is the number of independents, n the number of evaluations previously carried out on X)\n");
	_printf0_("      Gs are the G (gradient) values (m x n, where m is the number of independents, n the number of evaluations previously carried out on X,G)\n");
	_printf0_("      X - the new direction.\n");
	_printf0_("\n");
}/*}}}*/
#endif
WRAPPER(M1qn3_python){

	/*Boot module: */
	MODULEBOOT();

#ifdef _HAVE_M1QN3_
	/*input: */
	double* Xs=NULL;
	double* Gs=NULL;
	double* Js=NULL;
	int     intn;
	int     priorn;
	Data    data_struct;

	/* output: */
	double* X_out = NULL;

	/*intermediary: */
	double f;
	double dxmin=0.01;
	double gttol=0.0001;
	long   omode;

	/*checks on arguments on the matlab side: */
	if(nlhs!=NLHS || nrhs!=3){
		M1qn3Usage();
		_error_("M1qn3 usage error");
	}

	/*Input datasets: */
	FetchData(&Xs,&intn,&priorn,XHANDLE);
	FetchData(&Gs,&intn,&priorn,GHANDLE);
	FetchData(&Js,&priorn,JHANDLE);

	/*Initialize M1QN3 parameters*/
	SimulFunc costfuncion  = &fakesimul;  /*Cost function address*/
	void**    prosca       = &euclid_;  /*Dot product function (euclid is the default)*/
	char      normtype[]   = "dfn";     /*Norm type: dfn = scalar product defined by prosca*/
	long      izs[5];                   /*Arrays used by m1qn3 subroutines*/
	long      iz[5];                    /*Integer m1qn3 working array of size 5*/
	float     rzs[1];                   /*Arrays used by m1qn3 subroutines*/
	long      impres       = 0;         /*verbosity level*/
	long      imode[3]     = {0};       /*scaling and starting mode, 0 by default*/
	long      indic        = 4;         /*compute f and g*/
	long      reverse      = 0;         /*reverse or direct mode*/
	long      io           = 6;         /*Channel number for the output*/

	/*Optimization criterions*/
	long niter = long(intn); /*Maximum number of iterations*/
	long nsim  = long(intn); /*Maximum number of function calls*/

	/*Get problem dimension and initialize X, G and f: */
	long n = long(intn);
	IssmPDouble* G = xNew<IssmPDouble>(n); for (int i=0;i<n;i++)G[i]=Gs[i*priorn+0];
	IssmPDouble* X = xNew<IssmPDouble>(n); for (int i=0;i<n;i++)X[i]=Xs[i*priorn+0];
	f = Js[0];
	
	/*Allocate m1qn3 working arrays (see doc)*/
	long      m   = 100;
	long      ndz = 4*n+m*(2*n+1);
	double*   dz  = xNew<double>(ndz);
	double f1=f;

	/*Initialize: */
	data_struct.priorn=priorn;
	data_struct.counter=0;
	data_struct.Gs=Gs;
	data_struct.Js=Js;
	data_struct.Xs=Xs;

	m1qn3_(costfuncion,prosca,&ctonbe_,&ctcabe_,
				&n,X,&f,G,&dxmin,&f1,
				&gttol,normtype,&impres,&io,imode,&omode,&niter,&nsim,iz,dz,&ndz,
				&reverse,&indic,izs,rzs,(void*)&data_struct);

	switch(int(omode)){
		case 0: /* _printf0_("   Stop requested\n");*/ break;
		case 1:  _printf0_("   Convergence reached (gradient satisfies stopping criterion)\n"); break;
		case 2:  _printf0_("   Bad initialization\n"); break;
		case 3:  _printf0_("   Line search failure\n"); break;
		case 4:  _printf0_("   Maximum number of iterations exceeded\n");break;
		case 5:  _printf0_("   Maximum number of function calls exceeded\n"); break;
		case 6:  _printf0_("   stopped on dxmin during line search\n"); break;
		case 7:  _printf0_("   <g,d> > 0  or  <y,s> <0\n"); break;
		default: _printf0_("   Unknown end condition\n");
	}

	/*build output: */
	X_out=xNew<IssmPDouble>(n);
	for(int i=0;i<n;i++)X_out[i]=X[i];

	/*Write data: */
	WriteData(XOUT,X_out,n);

	/*end module: */
	xDelete<double>(Xs);
	xDelete<double>(Gs);
	xDelete<double>(Js);
	xDelete<double>(X_out);
	xDelete<double>(G);
	xDelete<double>(X);
	#else
	_error_("m1qn3 is not installed");
	#endif
	MODULEEND();

}


#ifdef _HAVE_M1QN3_
void fakesimul(long* indic,long* n,double* X,double* pf,double* G,long izs[1],float rzs[1],void* dzs){

	Data* ds=(Data*)dzs;
	double* Xs=ds->Xs;
	double* Gs=ds->Gs;
	double* Js=ds->Js;
	
	/*Are we done? : */
	if(ds->counter+1==ds->priorn){
		*indic=0;
		return;
	}
	else{
		ds->counter++;
		*pf=Js[ds->counter];
		for(int i=0;i<*n;i++) X[i]=Xs[ds->priorn*i+ds->counter];
		for(int i=0;i<*n;i++) G[i]=Gs[ds->priorn*i+ds->counter];
	}
}
#endif
