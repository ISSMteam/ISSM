/*
	MeshProfileIntersectionx.h
*/

#ifndef _MESHPROFILEINTERSECTIONX_H
#define _MESHPROFILEINTERSECTIONX_H

#include "../../shared/shared.h"
#include "../../classes/classes.h"

/* local prototypes: */
void MeshProfileIntersectionx(double** psegments, int* pnumseg, int* index, double* x, double* y, int nel, int nods,  Contour<IssmPDouble>** contours,int numcontours);
void MeshSegmentsIntersection(double** psegments, int* pnumsegs,int* index, double* x, double* y, int nel, int nods, double* xc, double* yc, int numnodes);
void ElementSegmentsIntersection(DataSet* segments_dataset,int el, double* xnodes,double* ynodes,double* xc,double* yc,int numnodes);
void ElementSegment(DataSet* segments_dataset,int el,int contouri, double* xnodes,double* ynodes,double* xsegment,double* ysegment);
int  SegmentIntersect(double* palpha, double* pbeta, double* x1, double* y1, double* x2, double* y2);
bool NodeInElement(double* xnodes, double* ynodes, double x, double y);
bool IsIdenticalNode(double x1, double y1, double x2, double y2, double tolerance);

#endif /* _MESHPROFILEINTERSECTIONX_H */
