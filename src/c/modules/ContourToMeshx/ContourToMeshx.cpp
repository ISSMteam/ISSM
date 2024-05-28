/*! \file  ContourToMeshx.c
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./ContourToMeshx.h"

int ContourToMeshx(double** pin_nod,double** pin_elem, double* index, double* x, double* y,Contours* contours,char* interptype,int nel,int nods, int edgevalue) {

	/*Contour:*/
	double value;

	/*output: */
	double*  in_nod;
	double*  in_elem;
	in_nod   = xNewZeroInit<double>(nods);
	in_elem  = xNewZeroInit<double>(nel);

	/*initialize thread parameters: */
	ContourToMeshxThreadStruct gate;
	gate.contours  = contours;
	gate.nods      = nods;
	gate.edgevalue = edgevalue;
	gate.in_nod    = in_nod;
	gate.x         = x;
	gate.y         = y;

	/*launch the thread manager with ContourToMeshxt as a core: */
	LaunchThread(ContourToMeshxt,(void*)&gate,_NUMTHREADS_);

	/*Take care of the case where an element interpolation has been requested: */
	if ((strcmp(interptype,"element")==0) || (strcmp(interptype,"element and node")==0)){
		for(int n=0;n<nel;n++){
			if ( (in_nod[ (int)*(index+3*n+0) -1] == 1) && (in_nod[ (int)*(index+3*n+1) -1] == 1) && (in_nod[ (int)*(index+3*n+2) -1] == 1) ){
				value=1.; in_elem[n]=value;
			}
		}
	}

	/*Assign output pointers: */
	*pin_nod=in_nod;
	*pin_elem=in_elem;

	return 1;
}
