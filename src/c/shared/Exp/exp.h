/* \file exp.h
 * \brief: header file for contour (argus type, files in .exp extension) operations
 */

#ifndef _EXP_H_
#define _EXP_H_

#include <cstring>
#include "../Numerics/recast.h"
#include "../Exceptions/exceptions.h"
#include "../MemOps/MemOps.h"

int IsInPolySerial(double* in,double* xc,double* yc,int numvertices,double* x,double* y,int nods, int edgevalue);
int ExpWrite(int nprof,int* profnvertices,double** pprofx,double** pprofy,char* domainname);
int pnpoly(int npol, double *xp, double *yp, double x, double y, int edgevalue);

template <class doubletype> int IsInPoly(doubletype* in,double* xc,double* yc,int numvertices,double* x,double* y,int i0,int i1, int edgevalue){ /*{{{*/

	int i;
	double x0,y0;
	doubletype value;
	double xmin=xc[0];
	double xmax=xc[0];
	double ymin=yc[0];
	double ymax=yc[0];

	/*Get extrema*/
	for (i=1;i<numvertices;i++){
		if(xc[i]<xmin) xmin=xc[i];
		if(xc[i]>xmax) xmax=xc[i];
		if(yc[i]<ymin) ymin=yc[i];
		if(yc[i]>ymax) ymax=yc[i];
	}

	/*Go through all vertices of the mesh:*/
	for (i=i0;i<i1;i++){

		//Get current value of value[i] -> do not change it if != 0
		value=in[i];
		if (reCast<bool,doubletype>(value)){
			/*this vertex already is inside one of the contours, continue*/
			continue;
		}

		/*pick up vertex (x[i],y[i]) and figure out if located inside contour (xc,yc)*/
		x0=x[i]; y0=y[i];
		if(x0<xmin || x0>xmax || y0<ymin || y0>ymax){
			value=0;
		}
		else{
			value=pnpoly(numvertices,xc,yc,x0,y0,edgevalue);
		}
		in[i]=value;
	}
	 return 1;
}/*}}}*/
template <class doubletype> int ExpRead(int* pnprof,int** pprofnvertices,doubletype*** ppprofx,doubletype*** ppprofy,bool** pclosed,char* domainname){ /*{{{*/

	/*indexing: */
	int i,counter;

	/*I/O: */
	FILE   *fid = NULL;
	char    chardummy[256];
	double  ddummy;

	/*output: */
	int          nprof;                //number of profiles in the domainname file
	int         *profnvertices = NULL; //array holding the number of vertices for the nprof profiles
	doubletype **pprofx        = NULL; //array of profiles x coordinates
	doubletype **pprofy        = NULL; //array of profiles y coordinates
	bool        *closed        = NULL; //array holding closed flags for the nprof profiles

	/*For each profile: */
	int         n;
	doubletype *x  = NULL;
	doubletype *y  = NULL;
	bool        cl;

	/*open domain outline file for reading: */
	if ((fid=fopen(domainname,"r"))==NULL){
		_error_("could not find file \"" << domainname<<"\". Make sure that the file and path provided exist.");
	}

	/*Do a first pass through the domainname file, to figure out how many profiles we need to read: */
	nprof=1;
	for(;;){
		//## Name:filename
		if(fscanf(fid,"%255s %255s\n",chardummy,chardummy)!=2) _error_("Could not read " << domainname);
		//## Icon:0
		if(fscanf(fid,"%255s %255s\n",chardummy,chardummy)!=2) _error_("Could not read " << domainname<<"(Expecting ## Icon:0 and read "<<chardummy<<")");
		//# Points Count Value
		if(fscanf(fid,"%255s %255s %255s %255s\n",chardummy,chardummy,chardummy,chardummy)!=4) _error_("Could not read " << domainname);
		if(fscanf(fid,"%20i %255s\n",&n,chardummy)!=2) _error_("Could not read number of points in "<<domainname);
		//# X pos Y pos
		if(fscanf(fid,"%255s %255s %255s %255s %255s\n",chardummy,chardummy,chardummy,chardummy,chardummy)!=5) _error_("Could not read " << domainname);
		for (i=0;i<n;i++){
			if(fscanf(fid,"%30lf %30lf\n",&ddummy,&ddummy)!=2){
				_error_("Could not read coordinate of vertex "<< i <<" of "<<domainname);
			}
		}
		/*check whether we are at the end of the file, otherwise, keep reading next profile:*/
		if(feof(fid)) break;
		nprof++;
	}

	/*Allocate and initialize all the profiles: */
	profnvertices = xNew<int>(nprof);
	pprofx        = xNew<doubletype*>(nprof);
	pprofy        = xNew<doubletype*>(nprof);
	for (i=0;i<nprof;i++){
		pprofx[i] = NULL;
		pprofy[i] = NULL;
	}
	closed=xNew<bool>(nprof);

	/*Reset file pointer to beginning of file: */
	fseek(fid,0,SEEK_SET);

	/*Start reading profiles: */
	for(counter=0;counter<nprof;counter++){

		/*Skip header: */
		//## Name:filename
		if(fscanf(fid,"%255s %255s\n",chardummy,chardummy)!=2) _error_("Could not read " << domainname);
		//## Icon:0
		if(fscanf(fid,"%255s %255s\n",chardummy,chardummy)!=2) _error_("Could not read " << domainname);
		//# Points Count Value
		if(fscanf(fid,"%255s %255s %255s %255s\n",chardummy,chardummy,chardummy,chardummy)!=4) _error_("Could not read " << domainname);

		/*Get number of profile vertices: */
		if(fscanf(fid,"%20i %255s\n",&n,chardummy)!=2) _error_("Could not read number of points in "<<domainname);

		/*Skip next line: */
		//# X pos Y pos
		if(fscanf(fid,"%255s %255s %255s %255s %255s\n",chardummy,chardummy,chardummy,chardummy,chardummy)!=5) _error_("Could not read " << domainname);

		/*Allocate vertices: */
		x=xNew<doubletype>(n);
		y=xNew<doubletype>(n);

		/*Read vertices: */
		for (i=0;i<n;i++){
			if(fscanf(fid,"%30lf %30lf\n",&x[i],&y[i])!=2){
				_error_("Could not read coordinate of vertex "<<i<<" of "<<domainname);
			}
		}

		/*Now check that we are dealing with open contours: */
		cl=false;
		if((x[0]==x[n-1]) && (y[0]==y[n-1])){
			cl=true;
		}

		/*Assign pointers: */
		profnvertices[counter]=n;
		pprofx[counter]=x;
		pprofy[counter]=y;
		closed[counter]=cl;
	}

	/*close domain outline file: */
	fclose(fid);

	/*Assign output pointers: */
	*pnprof=nprof;
	*pprofnvertices=profnvertices;
	*ppprofx=pprofx;
	*ppprofy=pprofy;
	if(pclosed) *pclosed=closed;
	else         xDelete<bool>(closed);
	return 1;

} /*}}}*/

#endif
