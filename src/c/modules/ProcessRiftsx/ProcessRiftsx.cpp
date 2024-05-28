/*!\file:  ProcessRifts.cpp
 * \brief split a mesh where a rift (or fault) is present
 */ 

#include "./ProcessRiftsx.h"
#include "../../classes/RiftStruct.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ProcessRiftsx(int** pindex, int* pnel,double** px,double** py,int* pnods,int** psegments,int** psegmentmarkers,int *pnum_seg,RiftStruct **priftstruct){

	/*Output*/
	int      numrifts,numrifts0;
	int     *riftsnumsegments     = NULL;
	int    **riftssegments        = NULL;
	int     *riftsnumpairs        = NULL;
	int    **riftspairs           = NULL;
	int     *riftstips            = NULL;
	double **riftspenaltypairs    = NULL;
	int     *riftsnumpenaltypairs = NULL;

	/*Recover initial mesh*/
	int     nel            = *pnel;
	int    *index          = *pindex;
	double *x              = *px;
	double *y              = *py;
	int     nods           = *pnods;
	int    *segments       = *psegments;
	int    *segmentmarkers = *psegmentmarkers;
	int     num_seg        = *pnum_seg;

	/*Intermediary*/
	int     riftflag;

	/*First, do some fixing on the existing mesh: we do not want any element belonging entirely to the segment list (ie: 
	 *all the nodes of this element belong to the segments (tends to happen when there are corners: */
	RemoveCornersFromRifts(&index,&nel,&x,&y,&nods,segments,segmentmarkers,num_seg);

	/*Figure out if we have rifts, and how many: */
	IsRiftPresent(&riftflag,&numrifts0,segmentmarkers,num_seg);

	if(!riftflag) _error_("No rift present in mesh");

	/*Split mesh*/
	SplitMeshForRifts(&nel,&index,&nods,&x,&y,&num_seg,&segments,&segmentmarkers);

	/*Order segments so that their normals point outside the domain: */
	OrderSegments(&segments,num_seg, index,nel);

	/*We do not want to output segments mixed with rift segments: wring out the rifts from the segments, using the 
	 *segmentmarkerlist:*/
	SplitRiftSegments(&segments,&segmentmarkers,&num_seg,&numrifts,&riftsnumsegments,&riftssegments,numrifts0,nods,nel);

	/*Using rift segments, associate rift faces in pairs, each pair face representing opposite flanks of the rifts facing one another directly: */
	PairRiftElements(&riftsnumpairs,&riftspairs,numrifts,riftsnumsegments,riftssegments,x,y);

	/*Order rifts so that they start from one tip, go to the other tip, and back: */
	OrderRifts(&riftstips,riftssegments,riftspairs,numrifts,riftsnumsegments,x,y,nods,nel);

	/*Create penalty pairs, used by Imp: */
	PenaltyPairs(&riftspenaltypairs,&riftsnumpenaltypairs,numrifts,riftssegments,riftsnumsegments,riftspairs,riftstips,x,y);

	/*Create Riftstruct*/
	RiftStruct* riftstruct = new RiftStruct(numrifts,riftsnumsegments,riftssegments,riftsnumpairs,riftspairs,riftsnumpenaltypairs,riftspenaltypairs,riftstips);

	/*Assign output pointers for mesh*/
	*pnel            = nel;
	*pindex          = index;
	*px              = x;
	*py              = y;
	*pnods           = nods;
	*psegments       = segments;
	*psegmentmarkers = segmentmarkers;
	*pnum_seg        = num_seg;

	/*Assign output pointers for rifts*/
	*priftstruct = riftstruct;
}
