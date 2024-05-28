/*! \file  ContourToNodesx.c
 */

#include "./ContourToNodesx.h"

int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contour<IssmPDouble>** contours,int numcontours,int edgevalue){

	int i;

	/*Contour:*/
	Contour<IssmPDouble>* contouri=NULL;
	int      numnodes;
	double*  xc=NULL;
	double*  yc=NULL;

	/*output: */
	IssmPDouble* flags=NULL;
	flags=xNew<IssmPDouble>(nods);

	/*Loop through all contours: */
	for (i=0;i<numcontours;i++){
		contouri=*(contours+i);
		numnodes=contouri->nods;
		xc=contouri->x;
		yc=contouri->y;
		IsInPoly(flags,xc,yc,numnodes,x,y,0,nods,edgevalue);
	}

	/*Assign output pointers: */
	*pflags=flags;
	return 1;
}

int ContourToNodesx(IssmPDouble** pflags,double* x, double* y, int nods, Contours* contours, int edgevalue){

	/*output: */
	IssmPDouble* flags=NULL;
	flags=xNewZeroInit<IssmPDouble>(nods);

	/*Loop through all contours: */
	if(contours){
		for(Object* & object:contours->objects){
			Contour<IssmPDouble>* contour=(Contour<IssmPDouble>*)object;
			IsInPoly(flags,contour->x,contour->y,contour->nods,x,y,0,nods,edgevalue);
		}
	}

	/*Assign output pointers: */
	*pflags=flags;
	return 1;
}
