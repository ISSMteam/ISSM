/*\file ElementConnectivity.c
 *\brief: build element connectivity using node connectivity and elements. 
 */

#include "./ElementConnectivity.h"

void ElementConnectivityUsage(void) {/*{{{*/
	_printf0_("\n");
	_printf0_("   usage: elementconnectivity = " << __FUNCT__ << "(elements, nodeconnectivity);\n");
	_printf0_("\n");
}/*}}}*/
WRAPPER(ElementConnectivity_python){

	/*inputs: */
	int* elements=NULL;
	int* nodeconnectivity=NULL;
	int  nels,nods;
	int  width;

	/*outputs: */
	int* elementconnectivity=NULL;

	/*Boot module: */
	MODULEBOOT();

	/*checks on arguments: */
	CHECKARGUMENTS(NLHS,NRHS,&ElementConnectivityUsage);

	/*Input datasets: */
	FetchData(&elements,&nels,NULL,ELEMENTS);
	FetchData(&nodeconnectivity,&nods,&width,NODECONNECTIVITY);

	/*!Generate internal degree of freedom numbers: */
	ElementConnectivityx(&elementconnectivity,elements,nels,nodeconnectivity,nods,width);

	/*write output datasets: */
	WriteData(ELEMENTCONNECTIVITY,elementconnectivity,nels,3);

	/*Clean up*/
	xDelete<int>(elements);
	xDelete<int>(nodeconnectivity);
	xDelete<int>(elementconnectivity);

	/*end module: */
	MODULEEND();
}
