/*\file NodeConnectivity.c
 *\brief: build node connectivity from elements. 
 */

#include "./NodeConnectivity.h"

void NodeConnectivityUsage(void){/*{{{*/
	_printf0_("\n");
	_printf0_("   usage: connectivity = " << __FUNCT__ << "(elements, numnodes);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(NodeConnectivity_python){

	/*inputs: */
	int* elements=NULL;
	int  nels;
	int  nods;

	/*outputs: */
	int* connectivity=NULL;
	int  width;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	CHECKARGUMENTS(NLHS,NRHS,&NodeConnectivityUsage);

	/*Input datasets: */
	FetchData(&elements,&nels,NULL,ELEMENTS);
	FetchData(&nods,NUMNODES);

	/*!Generate internal degree of freedom numbers: */
	NodeConnectivityx(&connectivity,&width,elements,nels,nods);

	/*write output datasets: */
	WriteData(CONNECTIVITY,connectivity,nods,width);

	/*end module: */
	xDelete<int>(elements);
	xDelete<int>(connectivity);
	MODULEEND();
}
