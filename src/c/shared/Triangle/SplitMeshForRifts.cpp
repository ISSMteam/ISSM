/*
 * SplitMeshForRifts.c:
 */
#include "./triangle.h"
#include "../MemOps/MemOps.h"

int SplitMeshForRifts(int* pnel,int** pindex,int* pnods,double** px,double** py,int* pnsegs,int** psegments,int** psegmentmarkerlist){

	/*Some notes on dimensions: 
	index  of size nelx3
	x and y of size nodsx1
	segments of size nsegsx3*/

	int i,j,k,l;
	int node;
	int el;
	int  nriftsegs;
	int* riftsegments=NULL; 
	int* flags=NULL;
	int  NumGridElementListOnOneSideOfRift;
	int* GridElementListOnOneSideOfRift=NULL;

	/*Recover input: */
	int     nel               = *pnel;
	int    *index             = *pindex;
	int     nods              = *pnods;
	double *x                 = *px;
	double *y                 = *py;
	int     nsegs             = *pnsegs;
	int    *segments          = *psegments;
	int    *segmentmarkerlist = *psegmentmarkerlist;

	/*Establish list of segments that belong to a rift: */
	/*riftsegments of size nriftsegsx4 (4 for first element on segment,second element,first node and second snode)*/
	RiftSegmentsFromSegments(&nriftsegs,&riftsegments,nel,index,nsegs,segments);

	/*Go through all nodes of the rift segments, and start splitting the mesh: */
	flags=xNewZeroInit<int>(nods); //to make sure we don't split the same nodes twice!
	for (i=0;i<nriftsegs;i++){
		for (j=0;j<2;j++){

			node=riftsegments[4*i+j+2];
			if(flags[node-1]){
				/*This node was already split, skip:*/
				continue;
			}
			else{
				flags[node-1]=1;
			}

			if(IsGridOnRift(riftsegments,nriftsegs,node)){

				DetermineGridElementListOnOneSideOfRift(&NumGridElementListOnOneSideOfRift,&GridElementListOnOneSideOfRift,i,nriftsegs,riftsegments,node,index,nel);

				/*Summary: we have for node, a list of elements
				 * (GridElementListOnOneSideOfRift, of size
				 * NumGridElementListOnOneSideOfRift) that all contain node 
				 *and that are on the same side of the rift. For all these
				 elements, we clone node into another node, and we swap all
				 instances of node in the triangulation *for those elements, to the
				 new node.*/

				//create new node
				x=xReNew<double>(x,nods,nods+1);
				y=xReNew<double>(y,nods,nods+1);
				x[nods]=x[node-1]; //matlab indexing
				y[nods]=y[node-1]; //matlab indexing

				//augment number of nodes 
				nods++;

				//change elements owning this node
				for (k=0;k<NumGridElementListOnOneSideOfRift;k++){
					el=GridElementListOnOneSideOfRift[k];
					for (l=0;l<3;l++){
						if (index[3*el+l]==node) index[3*el+l]=nods; //again, matlab indexing.
					}
				}
			}// if(IsGridOnRift(riftsegments,nriftsegs,node))
		} //for(j=0;j<2;j++)
	} //for (i=0;i<nriftsegs;i++)

	/*update segments: they got modified completely by adding new nodes.*/
	UpdateSegments(&segments,&segmentmarkerlist, &nsegs,index,x,y,riftsegments,nriftsegs,nods,nel);

	/*Assign output pointers: */
	*pnel=nel;
	*pindex=index;
	*pnods=nods;
	*px=x;
	*py=y;
	*pnsegs=nsegs;
	*psegments=segments;
	*psegmentmarkerlist=segmentmarkerlist;
	return 1;
}
