/*!\file:  Exp.cpp
 * \brief Exp.cpp: write the vertex coordinates defined in a domain 
 * outline from Argus (.exp file). The first contour in the file is for 
 * the outside domain outline. 
 */
#include <stdio.h>
#include <math.h>

#include "../MemOps/MemOps.h"
#include "../Exceptions/exceptions.h"
#include "./exp.h"

int ExpWrite(int nprof,int* profnvertices,double** pprofx,double** pprofy,char* domainname){/*{{{*/

	/*I/O: */
	FILE* fid=NULL;

	/*open domain outline file for writing: */
	if((fid=fopen(domainname,"w"))==NULL) _error_("could not open domain file " << domainname); 

	/*Start writing profiles: */
	for(int counter=0;counter<nprof;counter++){

		/*Write header: */
		fprintf(fid,"## Name:%s\n",domainname);
		fprintf(fid,"## Icon:0\n");
		fprintf(fid,"# Points Count	Value\n");
		fprintf(fid,"%u %s\n",profnvertices[counter]  ,"1.");
		fprintf(fid,"# X pos	Y pos\n");

		/*Write vertices: */
		for(int i=0;i<profnvertices[counter];i++){
			fprintf(fid,"%lf\t%lf\n",pprofx[counter][i],pprofy[counter][i]);
		}

		/*Write blank line: */
		if(counter<nprof-1) fprintf(fid,"\n");
	}

	/*close Exp file: */
	fclose(fid);

	return 1;
}/*}}}*/
int IsInPolySerial(double* in,double* xc,double* yc,int numvertices,double* x,double* y,int nods, int edgevalue){ /*{{{*/

	double x0,y0;

	/*Go through all vertices of the mesh:*/
	for(int i=0;i<nods;i++){
		if (in[i]){
			/*this vertex already is inside one of the contours, continue*/
			continue;
		}
		/*pick up vertex: */
		x0=x[i];
		y0=y[i];
		if (pnpoly(numvertices,xc,yc,x0,y0,edgevalue)){
			in[i]=1.;
		}
		else{
			in[i]=0.;
		}
	}

	return 1;
}
/*}}}*/
/*pnpoly{{{*/
int pnpoly(int npol, double *xp, double *yp, double x, double y, int edgevalue) {
	int i, j, c = 0;
	double n1, n2, normp, scalar;

	/*first test, are they colinear? if yes, is the point in the middle of the segment*/
	if (edgevalue != 2 ){
		for (i = 0, j = npol-1; i < npol; j = i++) {
			n1=pow(yp[i]-yp[j],2.0)+pow(xp[i]-xp[j],2.0);
			n2=pow(y-yp[j],2.0)+pow(x-xp[j],2.0);
			normp=pow(n1*n2,0.5);
			scalar=(yp[i]-yp[j])*(y-yp[j])+(xp[i]-xp[j])*(x-xp[j]);

			if (scalar == normp){
				if (n2<=n1){
					c = edgevalue;
					return c;
				}
			}
		}
	}
	/*second test : point is neither on a vertex, nor on a side, where is it ?*/
	for (i = 0, j = npol-1; i < npol; j = i++) {
		if ((((yp[i]<=y) && (y<yp[j])) ||
					((yp[j]<=y) && (y<yp[i]))) &&
				(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i])){
			c = !c;
		}
	}
	return c;
}/*}}}*/
