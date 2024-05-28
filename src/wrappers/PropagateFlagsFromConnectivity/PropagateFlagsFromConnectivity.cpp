/*\file PropagateFlagsFromConnectivity.c
 *\brief: propagate flags onto mesh, element by element, using connectivity.
 */

#include "./PropagateFlagsFromConnectivity.h"

void PropagateFlagsFromConnectivityUsage(void) {/*{{{*/
	_printf_("\n");
	_printf_("   usage: [pool] = " << __FUNCT__ << "(connectivity,pool,index,flags);\n");;
	_printf_("\n");
}/*}}}*/
WRAPPER(PropagateFlagsFromConnectivity_python){

	/*input/output datasets: */
	double* connectivity=NULL;
	int     nel;
	double* pool=NULL;
	double* flags=NULL;
	int     index;
	int     dummy;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments on the matlab side: */
	CheckNumMatlabArguments(nlhs,NLHS,nrhs,NRHS,__FUNCT__,&PropagateFlagsFromConnectivityUsage);

	/*Input datasets: */
	FetchData(&connectivity,&nel,&dummy,CONNECTIVITY);
	FetchData(&pool,&dummy,POOL);
	FetchData(&index,INDEX);
	FetchData(&flags,&dummy,FLAGS);

	/*!Generate internal degree of freedom numbers: */
	PropagateFlagsFromConnectivityx(pool,connectivity,index,flags);

	/*write output datasets: */
	WriteData(POOLOUT,pool,nel);

	/*Free resources: */
	xDelete<double>(connectivity);
	xDelete<double>(pool);
	xDelete<double>(flags);

	/*end module: */
	MODULEEND();
}
