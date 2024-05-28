/*!\file:  triangle.h
 * \brief
 */ 

#ifndef _SHARED_TRIANGLE_H
#define _SHARED_TRIANGLE_H

#include <stdio.h>
#include <math.h>

//#define REAL double //took  it out because it may conflict with stdlib.h defines. put back if necessary
int AssociateSegmentToElement(int** psegments,int nseg,int* index,int nel);
int OrderSegments(int** psegments,int nseg, int* index,int nel);
int GridInsideHole(double* px0,double* py0,int n,double* x,double* y);
int FindElement(int A,int B,int* index,int nel);
int SplitMeshForRifts(int* pnel,int** pindex,int* pnods,double** px,double** py,int* pnsegs,int** psegments,int** psegmentmarkerlist);
int IsGridOnRift(int* riftsegments, int nriftsegs, int node);
int GridElementsList(int** pGridElements, int* pNumGridElements,int node,double * index,int nel);
int IsNeighbor(int el1,int el2,int* index);
int IsOnRift(int el,int nriftsegs,int* riftsegments);
void RiftSegmentsFromSegments(int* pnriftsegs, int** priftsegments, int nel,int* index, int nsegs,int* segments);
int DetermineGridElementListOnOneSideOfRift(int* pNumGridElementListOnOneSideOfRift, int** pGridElementListOnOneSideOfRift,int segmentnumber, int nriftsegs,int* riftsegments, int node,int* index,int nel);
int UpdateSegments(int** psegments,int** psegmentmarkerlist, int* pnsegs,int* index, double* x,double* y,int* riftsegments,int nriftsegs,int nods,int nel);
int FindElement(double A,double B,int* index,int nel);
int IsRiftPresent(int* priftflag,int* pnumrifts,int* segmentmarkerlist,int nsegs);
int SplitRiftSegments(int** psegments,int** psegmentmarkerlist, int* pnumsegs, int* pnumrifts,int** priftsnumsegs,int*** priftssegments,int numrifts,int nods,int nels);
int OrderRifts(int** priftstips,int** riftssegments,int** riftspairs,int numrifts,int* riftsnumsegments,double* x,double* y,int nods,int nels);
int PenaltyPairs(double*** priftspenaltypairs,int** priftsnumpenaltypairs,int numrifts,int**  riftssegments,
		int* riftsnumsegments,int** riftspairs,int* riftstips,double* x,double* y);
int RemoveCornersFromRifts(int** pindex,int* pnel,double** px,double** py,int* pnods,int* segments,int* segmentmarkers,int num_seg);
int PairRiftElements(int** priftsnumpairs,int*** priftspairs,int numrifts,int* riftsnumsegments,int** riftssegments,double* x,double* y);

#endif  /* _SHARED_TRIANGLE_H */
