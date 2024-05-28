/*
 * TriangleUtils: mesh manipulation routines: 
 */

#include <stdio.h>

#include "./triangle.h"
#include "../Exceptions/exceptions.h"
#include "../MemOps/MemOps.h"

#define RIFTPENALTYPAIRSWIDTH 8
int IsGridOnRift(int* riftsegments, int nriftsegs, int node){/*{{{*/

	/*Does this node belong to 4 elements, or just 2? If it belongs to 4 elements, it is inside a rift, 
	 *if it belongs to 2 elements, it is on the tip of a rift, or it has already been split across the rift (see below).*/

	int i;
	int j;
	int count;

	count=0;
	for (i=0;i<nriftsegs;i++){
		for (j=0;j<2;j++){
			if ((*(riftsegments+4*i+2+j))==node) count++;
		}
	}
	if (count==2){
		return 1;
	}
	else{
		return 0;
	}
}/*}}}*/
int GridElementsList(int** pGridElements, int* pNumGridElements,int node,int* index,int nel){/*{{{*/

	/*From a node, recover all the elements that are connected to it: */
	int i,j;
	int noerr=1;

	int max_number_elements=12;
	int current_size;
	int NumGridElements;
	int* GridElements=NULL;
	int* GridElementsRealloc=NULL;

	/*From a mesh with 30 degrees minimum angle, we get 12 possible elements that own 
	 * the node. We start by allocating GridElements with that size, and realloc 
	 * more if needed.*/

	current_size=max_number_elements;
	NumGridElements=0;
	GridElements=xNew<int>(max_number_elements);

	for (i=0;i<nel;i++){
		for (j=0;j<3;j++){
			if (index[3*i+j]==node){
				if (NumGridElements<=(current_size-1)){
					GridElements[NumGridElements]=i;
					NumGridElements++;
					break;
				}
				else{
					/*Reallocate another max_number_elements slots in the GridElements: */
					GridElementsRealloc=xReNew<int>(GridElements,current_size,(current_size+max_number_elements));
					if (!GridElementsRealloc){
						noerr=0;
						goto cleanup_and_return;
					}
					current_size+=max_number_elements;
					GridElements=GridElementsRealloc;
					GridElements[NumGridElements]=i;
					NumGridElements++;
					break;
				}
			}
		}
	}
	cleanup_and_return:
	if(!noerr){
		xDelete<int>(GridElements);
	}
	/*Allocate return pointers: */
	*pGridElements=GridElements;
	*pNumGridElements=NumGridElements;
	return noerr;
}/*}}}*/
int IsNeighbor(int el1,int el2,int* index){/*{{{*/
	/*From a triangulation held in index, figure out if elements 1 and 2 have two nodes in common: */
	int i,j;
	int count=0;
	for (i=0;i<3;i++){
		for (j=0;j<3;j++){
			if (index[3*el1+i]==index[3*el2+j])count++;
		}
	}
	if (count==2){
		return 1;
	}
	else{
		return 0;
	}
}/*}}}*/
int IsOnRift(int el,int nriftsegs,int* riftsegments){/*{{{*/
	/*From a list of elements segments, figure out if el belongs to it: */
	int i;
	for (i=0;i<nriftsegs;i++){
		if ((*(riftsegments+4*i+0)==el) || (*(riftsegments+4*i+1)==el)){
			return 1;
		}
	}
	return 0;
}/*}}}*/
void RiftSegmentsFromSegments(int* pnriftsegs, int** priftsegments, int nel,int* index,int nsegs,int* segments){/*{{{*/

	int i,counter;
	int el,el2;

	int  nriftsegs;
	int* riftsegments=NULL;
	int* riftsegments_uncompressed=NULL; 
	int  element_nodes[3];

	/*Allocate segmentflags: */
	riftsegments_uncompressed=xNewZeroInit<int>(nsegs*5);

	/*Find the segments that belong to a rift: they are the ones that see two elements. The other ones belong to a boundary 
	 *or a hole: */
	nriftsegs=0;
	for (i=0;i<nsegs;i++){
		el=(int)*(segments+3*i+2)-1; //element found in AssociateSegmentToElements
		/*Temporarily set nodes belonging to the segments to -1 in the triangulation index, and 
		 *then  proceed to find another element that owns the segment. If we don't find it, we know 
		 *we are dealing with a boundary or hole, otherwise, we are dealing with a rift: */
		element_nodes[0]=*(index+3*el+0);
		element_nodes[1]=*(index+3*el+1);
		element_nodes[2]=*(index+3*el+2);

		index[3*el+0]=-1;
		index[3*el+1]=-1;
		index[3*el+2]=-1;

		el2=FindElement(*(segments+3*i+0),*(segments+3*i+1),index,nel); 

		/*Restore index: */
		index[3*el+0]=element_nodes[0];
		index[3*el+1]=element_nodes[1];
		index[3*el+2]=element_nodes[2];

		if (el2!=-1){
			/*el and el2 are on a segment rift, facing one another, plug them into riftsegments_uncompressed: */
		    riftsegments_uncompressed[5*i+0]=1;
		    riftsegments_uncompressed[5*i+1]=el;
		    riftsegments_uncompressed[5*i+2]=el2;
		    riftsegments_uncompressed[5*i+3]=segments[3*i+0];
			 riftsegments_uncompressed[5*i+4]=segments[3*i+1];
			 nriftsegs++;
		}
	}

	/*Compress riftsegments_uncompressed:*/
	riftsegments=xNew<int>(nriftsegs*4);
	counter=0;
	for (i=0;i<nsegs;i++){
		if (riftsegments_uncompressed[5*i+0]){
			riftsegments[counter*4+0]=riftsegments_uncompressed[5*i+1];
			riftsegments[counter*4+1]=riftsegments_uncompressed[5*i+2];
			riftsegments[counter*4+2]=riftsegments_uncompressed[5*i+3];
			riftsegments[counter*4+3]=riftsegments_uncompressed[5*i+4];
			counter++;
		}
	}
	xDelete<int>(riftsegments_uncompressed);

	/*Assign output pointers: */
	*priftsegments=riftsegments;
	*pnriftsegs=nriftsegs;
}/*}}}*/
int DetermineGridElementListOnOneSideOfRift(int* pNumGridElementListOnOneSideOfRift, int** pGridElementListOnOneSideOfRift, int segmentnumber, int nriftsegs, int* riftsegments, int node,int* index,int nel){/*{{{*/

	int noerr=1;
	int k,l,counter;
	int newel;

	int* GridElements=NULL;
	int  NumGridElements;

	/*Output: */
	int NumGridElementListOnOneSideOfRift;
	int* GridElementListOnOneSideOfRift=NULL;

	/*Build a list of all the elements connected to this node: */
	GridElementsList(&GridElements,&NumGridElements,node,index,nel);

	/*Figure out the list of elements  that are on the same side of the rift. To do so, we start from one 
	 * side of the rift and keep rotating in the same direction:*/
	GridElementListOnOneSideOfRift=xNew<int>(NumGridElements);
	//bootstrap the GridElementListOnOneSideOfRift by filling elements from riftsegments: */
	GridElementListOnOneSideOfRift[0]=*(riftsegments+4*segmentnumber+0); /*this one does not belong to the same side, but is just there 
															   for a rotation direction, we 'll take it out when we are 
															   done rotating*/
	GridElementListOnOneSideOfRift[1]=*(riftsegments+4*segmentnumber+1);
	counter=1;
	for (;;){
		/*Find neighbour of element GridElementListOnOneSideOfRift[counter], not 
		 * equal to GridElementListOnOneSideOfRift[counter-1]*/
		for (k=0;k<NumGridElements;k++){
			if(IsNeighbor(GridElements[k],GridElementListOnOneSideOfRift[counter],index)){
				/*Verify this element is not already in our list of element on the same side of the rift: */
				newel=1;
				for (l=0;l<=counter;l++){
					if (GridElements[k]==GridElementListOnOneSideOfRift[l]){
						newel=0;
						break;
					}
				}
				if (newel){
					counter++;
					GridElementListOnOneSideOfRift[counter]=GridElements[k];
					if (IsOnRift(GridElements[k],nriftsegs,riftsegments)){
						break;
					}
					k=-1;
				}
			}
		}
		/*Reduce counter by 1 and get rift of first element in GridElementListOnOneSideOfRift:*/
		NumGridElementListOnOneSideOfRift=counter;
		for (l=0;l<NumGridElementListOnOneSideOfRift;l++){
			GridElementListOnOneSideOfRift[l]=GridElementListOnOneSideOfRift[l+1];
		}
		break;
	}// for (;;)

	/*Free resources: */
	xDelete<int>(GridElements);
	/*Assign output pointers: */
	*pNumGridElementListOnOneSideOfRift=NumGridElementListOnOneSideOfRift;
	*pGridElementListOnOneSideOfRift=GridElementListOnOneSideOfRift;
	return noerr;
}/*}}}*/
int UpdateSegments(int** psegments,int** psegmentmarkerlist, int* pnsegs,int* index, double* x,double* y,int* riftsegments,int nriftsegs,int nods,int nel){/*{{{*/

	int noerr=1;
	int i,j,k;
	int el1,el2;

	int *segments          = NULL;
	int *segmentmarkerlist = NULL;
	int  nsegs;

	/*Recover input: */
	segments          = *psegments;
	segmentmarkerlist = *psegmentmarkerlist;
	nsegs             = *pnsegs;

	/*Reallocate segments: */
	segments         =xReNew<int>(segments,         nsegs*3,(nsegs+nriftsegs)*3);
	segmentmarkerlist=xReNew<int>(segmentmarkerlist,nsegs,(nsegs+nriftsegs));

	/*First, update the existing segments to the new nodes :*/
	for (i=0;i<nriftsegs;i++){
		el1=riftsegments[4*i+0];
		el2=riftsegments[4*i+1];
		for (j=0;j<nsegs;j++){
			if (segments[3*j+2]==(el1+1)){
				/*segment j is the same as rift segment i.Let's update segments[j][:] using  element el1 and the corresponding rift segment.
				 *Because riftsegments does not represent a list of rift segments anymore (it got heavily modified in SplitElementsForRifts, 
				 *we can only rely on the position (x,y) of the rift nodes to create a segment:*/
				for (k=0;k<3;k++){
					if ((x[*(index+el1*3+k)-1]==x[*(segments+3*j+0)-1]) && (y[*(index+el1*3+k)-1]==y[*(segments+3*j+0)-1])){
						*(segments+3*j+0)=*(index+el1*3+k); _assert_(segments[3*j+0]<nods+1);
						break;
					}
				}
				for (k=0;k<3;k++){
					if ((x[*(index+el1*3+k)-1]==x[*(segments+3*j+1)-1])  && (y[*(index+el1*3+k)-1]==y[*(segments+3*j+1)-1])){
						*(segments+3*j+1)=*(index+el1*3+k); _assert_(segments[3*j+1]<nods+1);
						break;
					}
				}
				/*Deal with el2: */
				*(segments+3*(nsegs+i)+2)=el2+1;
				*(segmentmarkerlist+(nsegs+i))=*(segmentmarkerlist+j);
				for (k=0;k<3;k++){
					if ((x[*(index+el2*3+k)-1]==x[*(segments+3*j+0)-1]) && (y[*(index+el2*3+k)-1]==y[*(segments+3*j+0)-1])){
						*(segments+3*(nsegs+i)+0)=*(index+el2*3+k); _assert_(segments[3*(nsegs+i)+0]<nods+1);
						break;
					}
				}
				for (k=0;k<3;k++){
					if ((x[*(index+el2*3+k)-1]==x[*(segments+3*j+1)-1]) && (y[*(index+el2*3+k)-1]==y[*(segments+3*j+1)-1])){
						*(segments+3*(nsegs+i)+1)=*(index+el2*3+k); _assert_(segments[3*(nsegs+i)+1]<nods+1);
						break;
					}
				}
			}
			if (*(segments+3*j+2)==(el2+1)){
				/*segment j is the same as rift segment i.*/
				/*Let's update segments[j][:] using  element el2 and the corresponding rift segment: */
				for (k=0;k<3;k++){
					if ((x[*(index+el2*3+k)-1]==x[*(segments+3*j+0)-1]) && (y[*(index+el2*3+k)-1]==y[*(segments+3*j+0)-1])){
						*(segments+3*j+0)=*(index+el2*3+k); _assert_(segments[3*j+0]<nods+1);
						break;
					}
				}
				for (k=0;k<3;k++){
					if ((x[*(index+el2*3+k)-1]==x[*(segments+3*j+1)-1]) && (y[*(index+el2*3+k)-1]==y[*(segments+3*j+1)-1])){
						*(segments+3*j+1)=*(index+el2*3+k);_assert_(segments[3*j+1]<nods+1);
						break;
					}
				}
				/*Deal with el1: */
				*(segments+3*(nsegs+i)+2)=el1+1;
				*(segmentmarkerlist+(nsegs+i))=*(segmentmarkerlist+j);
				for (k=0;k<3;k++){
					if ((x[*(index+el1*3+k)-1]==x[*(segments+3*j+0)-1]) && (y[*(index+el1*3+k)-1]==y[*(segments+3*j+0)-1])){
						*(segments+3*(nsegs+i)+0)=*(index+el1*3+k);_assert_(segments[3*(nsegs+i)+0]<nods+1);
						break;
					}
				}
				for (k=0;k<3;k++){
					if ((x[*(index+el1*3+k)-1]==x[*(segments+3*j+1)-1]) && (y[*(index+el1*3+k)-1]==y[*(segments+3*j+1)-1])){
						*(segments+3*(nsegs+i)+1)=*(index+el1*3+k);_assert_(segments[3*(nsegs+i)+1]<nods+1);
						break;
					}
				}
			}
		}
	}
	nsegs+=nriftsegs;

	/*Assign output pointers: */
	*psegments=segments;
	*psegmentmarkerlist=segmentmarkerlist;
	*pnsegs=nsegs;

	return noerr;
}/*}}}*/
int FindElement(int A,int B,int* index,int nel){/*{{{*/

	int el=-1;
	for (int n=0;n<nel;n++){
		if(((index[3*n+0]==A) || (index[3*n+1]==A) || (index[3*n+2]==A)) && ((index[3*n+0]==B) || (index[3*n+1]==B) || (index[3*n+2]==B))){
			el=n;
			break;
		}
	}
	return el;
}/*}}}*/
int SplitRiftSegments(int** psegments,int** psegmentmarkerlist, int* pnumsegs, int* pnumrifts,int** priftsnumsegs,int*** priftssegments,int numrifts,int nods,int nel){/*{{{*/

	/*Using segment markers, wring out the rift segments from the segments. Rift markers are 
	 *of the form 2+i where i=0 to number of rifts */

	int noerr=1;
	int i,j,counter;

	/*input: */
	int *segments          = NULL;
	int *segmentmarkerlist = NULL;
	int numsegs;

	/*output: */
	int   new_numsegs;
	int  *riftsnumsegs       = NULL;
	int **riftssegments      = NULL;
	int  *new_segments       = NULL;
	int  *new_segmentmarkers = NULL;

	/*intermediary: */
	int* riftsegment=NULL;

	/*Recover input arguments: */
	segments          = *psegments;
	numsegs           = *pnumsegs;
	segmentmarkerlist = *psegmentmarkerlist;

	/*First, figure out  how many segments will be left in 'segments': */
	counter=0;
	for (i=0;i<numsegs;i++){
		if (segmentmarkerlist[i]==1)counter++; //1 is default marker for non-rifts;
	}
	/*Allocate new segments: */
	new_numsegs=counter;
	new_segments=xNew<int>(new_numsegs*3);
	new_segmentmarkers=xNew<int>(new_numsegs);

	/*Copy new segments info : */
	counter=0;
	for (i=0;i<numsegs;i++){
		if (segmentmarkerlist[i]==1){
			new_segments[3*counter+0]=segments[3*i+0];
			new_segments[3*counter+1]=segments[3*i+1];
			new_segments[3*counter+2]=segments[3*i+2];
			new_segmentmarkers[counter]=segmentmarkerlist[i];
			counter++;
		}
	}

	/*Now deal with rift segments: */
	riftsnumsegs=xNew<int>(numrifts);
	riftssegments=xNew<int*>(numrifts);
	for (i=0;i<numrifts;i++){
		/*Figure out how many segments for rift i: */
		counter=0;
		for (j=0;j<numsegs;j++){
			if (segmentmarkerlist[j]==2+i)counter++;
		}
		riftsnumsegs[i]=counter;
		riftsegment=xNew<int>(counter*3);
		/*Copy new segments info :*/
		counter=0;
		for (j=0;j<numsegs;j++){
			if (segmentmarkerlist[j]==(2+i)){
				riftsegment[3*counter+0]=segments[3*j+0];_assert_(riftsegment[3*counter+0]<nods+1);
				riftsegment[3*counter+1]=segments[3*j+1];_assert_(riftsegment[3*counter+1]<nods+1);
				riftsegment[3*counter+2]=segments[3*j+2];_assert_(riftsegment[3*counter+2]<nel+1);
				counter++;
			}
		}
		*(riftssegments+i)=riftsegment;
	}

	/*Free resources: */
	xDelete<int>(segments);

	/*Assign output pointers: */
	*psegments=new_segments;
	*psegmentmarkerlist=new_segmentmarkers;
	*pnumsegs=new_numsegs;
	*pnumrifts=numrifts;
	*priftssegments=riftssegments;
	*priftsnumsegs=riftsnumsegs;
	return noerr;
}/*}}}*/
int PairRiftElements(int** priftsnumpairs,int*** priftspairs,int numrifts,int* riftsnumsegments,int** riftssegments,double* x,double* y){/*{{{*/

	int noerr=1;
	int i,j,k;

	/*output: */
	int  *riftsnumpairs = NULL;
	int **riftspairs    = NULL;

	/*intermediary :*/
	int  numsegs;
	int* segments=NULL;
	int* pairs=NULL;
	int  node1,node2,node3,node4;

	riftsnumpairs=xNew<int>(numrifts);
	riftspairs=xNew<int*>(numrifts);
	for (i=0;i<numrifts;i++){
		segments=riftssegments[i];
		numsegs =riftsnumsegments[i];
		riftsnumpairs[i]=numsegs;
		pairs=xNew<int>(2*numsegs);
		for (j=0;j<numsegs;j++){
			pairs[2*j+0]=segments[3*j+2]; //retrieve element to which this segment belongs.
			node1=segments[3*j+0]-1; node2=segments[3*j+1]-1;
			/*Find element facing on other side of rift: */
			for (k=0;k<numsegs;k++){
				if (k==j)continue;
				node3=segments[3*k+0]-1; node4=segments[3*k+1]-1;
				/*We are trying to find 2 elements, where position of node3 == position of node1, and position of node4 == position of node2*/
				if (   ((x[node3]==x[node1]) && (y[node3]==y[node1]) && (x[node4]==x[node2]) && (y[node4]==y[node2]))
				    || ((x[node3]==x[node2]) && (y[node3]==y[node2]) && (x[node4]==x[node1]) && (y[node4]==y[node1]))  ){
					/*We found the corresponding element: */
					pairs[2*j+1]=segments[3*k+2];
					break;
				}
			}
		}
		riftspairs[i]=pairs;
	}

	/*Assign output pointers: */
	*priftsnumpairs=riftsnumpairs;
	*priftspairs=riftspairs;
	return noerr;
}/*}}}*/
int IsRiftPresent(int* priftflag,int* pnumrifts,int* segmentmarkerlist,int nsegs){/*{{{*/

	int i;
	int noerr=1;

	/*output: */
	int riftflag=0;
	int numrifts=0;

	int maxmark=1; //default marker for regular segments

	/*Any marker >=2 indicates a certain rift: */
	numrifts=0;
	for (i=0;i<nsegs;i++){
		if (segmentmarkerlist[i]>maxmark){
			numrifts++;
			maxmark=segmentmarkerlist[i];
		}
	}
	if(numrifts)riftflag=1;

	/*Assign output pointers:*/
	*priftflag=riftflag;
	*pnumrifts=numrifts;
	return noerr;
}/*}}}*/
int OrderRifts(int** priftstips,int** riftssegments,int** riftspairs,int numrifts,int* riftsnumsegments,double* x,double* y,int nods,int nels){/*{{{*/

	int noerr=1;
	int i,j,k,counter;

	/*intermediary: */
	int *riftsegments = NULL;
	int *riftpairs    = NULL;
	int numsegs;

	/*ordering and copy: */
	int *order             = NULL;
	int *riftsegments_copy = NULL;
	int *riftpairs_copy    = NULL;

	/*node and element manipulation: */
	int node1,node2,node3,node4,temp_node,tip1,tip2,node;
	int el2;
	int already_ordered=0;

	/*output: */
	int* riftstips=NULL;

	/*Allocate byproduct of this routine, riftstips: */
	riftstips=xNew<int>(numrifts*2);

	/*Go through all rifts: */
	for (i=0;i<numrifts;i++){
		riftsegments = riftssegments[i];
		riftpairs    = riftspairs[i];
		numsegs      = riftsnumsegments[i];

		/*Allocate copy of riftsegments and riftpairs, 
		 *as well as ordering vector: */
		riftsegments_copy=xNew<int>(numsegs*3);
		riftpairs_copy=xNew<int>(numsegs*2);
		order=xNew<int>(numsegs);

		/*First find the tips, using the pairs. If a pair of elements has one node in common, this node is a rift tip: */
		tip1=-1;
		tip2=-1;

		for (j=0;j<numsegs;j++){
			el2=*(riftpairs+2*j+1);
			node1=*(riftsegments+3*j+0);
			node2=*(riftsegments+3*j+1);
			/*Summary, el1 and el2 are facing one another across the rift. node1 and node2 belong to el1 and 
			 *are located on the rift. Find node3 and node4, nodes belonging to el2 and located on the rift: */
			for (k=0;k<numsegs;k++){
				if (*(riftsegments+3*k+2)==el2){
					node3=*(riftsegments+3*k+0);
					node4=*(riftsegments+3*k+1);
					break;
				}
			}
			/* Make sure node3 faces node1 and node4 faces node2: */
			_assert_(node1<nods+1 && node4<nods+1);
			_assert_(node1>0 && node4>0);
			if ((x[node1-1]==x[node4-1]) && (y[node1-1]==y[node4-1])){
				/*Swap node3 and node4:*/
				temp_node=node3;
				node3=node4;
				node4=temp_node;
			}

			/*Figure out if a tip is on this element: */
			if (node3==node1){
				/*node1 is a tip*/
				if (tip1==-1) {
					tip1=node1;
					continue;
				}
				if ((tip2==-1) && (node1!=tip1)){
					tip2=node1;
					break;
				}
			}

			if (node4==node2){
				/*node2 is a tip*/
				if (tip1==-1){
					tip1=node2;
					continue;
				}
				if ((tip2==-1) && (node2!=tip1)){
					tip2=node2;
					break;
				}
			}
		}

		/*Record tips in riftstips: */
		*(riftstips+2*i+0)=tip1;
		*(riftstips+2*i+1)=tip2;

		/*We have the two tips for this rift.  Go from tip1 to tip2, and figure out the order in which segments are sequential. 
		 *Because two elements are connected to tip1, we chose one first, which defines the direction we are rotating along the rift. */
		node=tip1;
		for (counter=0;counter<numsegs;counter++){
			for (j=0;j<numsegs;j++){
				node1=*(riftsegments+3*j+0);
				node2=*(riftsegments+3*j+1);

				if ((node1==node) || (node2==node)){
					/*Ok, this segment is connected to node, plug its index into order, unless we already plugged it before: */
					already_ordered=0;
					for (k=0;k<counter;k++){
						if(order[k]==j){
							already_ordered=1;
							break;
						}
					}
					if (!already_ordered){
						order[counter]=j;
						if(node1==node){
							node=node2;
						}
						else if(node2==node){
							node=node1;
						}
						break;
					}
				}
			}
		}

		/*Using the order vector, and the riftsegments_copy and riftspairs_copy, reorder the segments and the pairs: */
		for (j=0;j<numsegs;j++){
			_assert_(order[j]<numsegs);
			*(riftsegments_copy+3*j+0)=*(riftsegments+3*order[j]+0);
			*(riftsegments_copy+3*j+1)=*(riftsegments+3*order[j]+1);
			*(riftsegments_copy+3*j+2)=*(riftsegments+3*order[j]+2);
			*(riftpairs_copy+2*j+0)=*(riftpairs+2*order[j]+0);
			*(riftpairs_copy+2*j+1)=*(riftpairs+2*order[j]+1);
		}

		for (j=0;j<numsegs;j++){
			*(riftsegments+3*j+0)=*(riftsegments_copy+3*j+0);
			*(riftsegments+3*j+1)=*(riftsegments_copy+3*j+1);
			*(riftsegments+3*j+2)=*(riftsegments_copy+3*j+2);
			*(riftpairs+2*j+0)=*(riftpairs_copy+2*j+0);
			*(riftpairs+2*j+1)=*(riftpairs_copy+2*j+1);
		}

		xDelete<int>(order);
		xDelete<int>(riftsegments_copy);
		xDelete<int>(riftpairs_copy);

	}

	/*Assign output pointer:*/
	*priftstips=riftstips;
	return noerr;
}/*}}}*/
int PenaltyPairs(double*** priftspenaltypairs,int** priftsnumpenaltypairs,int numrifts,int** riftssegments,/*{{{*/
		int* riftsnumsegs,int** riftspairs,int* riftstips,double* x,double* y){

	int noerr=1;
	int i,j,k,k0;

	double el1,el2,node1,node2,node3,node4;
	double temp_node;

	/*output: */
	double **riftspenaltypairs    = NULL;
	double  *riftpenaltypairs     = NULL;
	int     *riftsnumpenaltypairs = NULL;

	/*intermediary: */
	int numsegs;
	int* riftsegments=NULL;
	int* riftpairs=NULL;
	int counter;
	double normal[2];
	double length;
	int    k1,k2;

	/*Allocate: */
	riftspenaltypairs=xNew<double*>(numrifts);
	riftsnumpenaltypairs=xNew<int>(numrifts);

	for(i=0;i<numrifts;i++){
		numsegs=riftsnumsegs[i];
		riftsegments=riftssegments[i];
		riftpairs=riftspairs[i];

		/*allocate riftpenaltypairs, and riftnumpenaltypairs: */
		if((numsegs/2-1)!=0)riftpenaltypairs=xNewZeroInit<double>((numsegs/2-1)*RIFTPENALTYPAIRSWIDTH);

		/*Go through only one flank of the rifts, not counting the tips: */
		counter=0;
		for(j=0;j<(numsegs/2);j++){
			el1=*(riftpairs+2*j+0);
			el2=*(riftpairs+2*j+1);
			node1=*(riftsegments+3*j+0);
			node2=*(riftsegments+3*j+1);
			/*Find segment index to recover node3 and node4, facing node1 and node2: */
			k0=-1;
			for(k=0;k<numsegs;k++){
				if(*(riftsegments+3*k+2)==el2){
					k0=k;
					break;
				}
			}
			node3=*(riftsegments+3*k0+0);
			node4=*(riftsegments+3*k0+1);

			/* Make sure node3 faces node1 and node4 faces node2: */
			if ((x[(int)node1-1]==x[(int)node4-1]) && (y[(int)node1-1]==y[(int)node4-1])){
				/*Swap node3 and node4:*/
				temp_node=node3;
				node3=node4;
				node4=temp_node;
			}	
			/*Ok, we have node1 facing node3, and node2 facing node4. Compute the normal to 
			 *this segment, and its length: */
			normal[0]=cos(atan2(x[(int)node1-1]-x[(int)node2-1],y[(int)node2-1]-y[(int)node1-1]));
			normal[1]=sin(atan2(x[(int)node1-1]-x[(int)node2-1],y[(int)node2-1]-y[(int)node1-1]));
			length=sqrt(pow(x[(int)node2-1]-x[(int)node1-1],(double)2)+pow(y[(int)node2-1]-y[(int)node1-1],(double)2));

			/*Be careful here, we want penalty loads on each node, not on each segment. This means we cannot plug node1,
			 * node2, node3 and node4 directly into riftpenaltypairs. We need to include node1, node2, node3 and node4, 
			 * only once. We'll add the normals and the lengths : */

			if(node1!=node3){ //exclude tips from loads
				k1=-1;
				for(k=0;k<counter;k++){
					if( (*(riftpenaltypairs+k*7+0))==node1){
						k1=k; 
						break;
					}
				}
				if(k1==-1){
					*(riftpenaltypairs+counter*7+0)=node1;
					*(riftpenaltypairs+counter*7+1)=node3;
					*(riftpenaltypairs+counter*7+2)=el1;
					*(riftpenaltypairs+counter*7+3)=el2;
					*(riftpenaltypairs+counter*7+4)=normal[0];
					*(riftpenaltypairs+counter*7+5)=normal[1];
					*(riftpenaltypairs+counter*7+6)=length/2;
					counter++;
				}
				else{
					*(riftpenaltypairs+k1*7+4)+=normal[0];
					*(riftpenaltypairs+k1*7+5)+=normal[1];
					*(riftpenaltypairs+k1*7+6)+=length/2;
				}
			}
			if(node2!=node4){
				k2=-1;
				for(k=0;k<counter;k++){
					if( (*(riftpenaltypairs+k*7+0))==node2){
						k2=k;
						break;
					}
				}
				if(k2==-1){
					*(riftpenaltypairs+counter*7+0)=node2;
					*(riftpenaltypairs+counter*7+1)=node4;
					*(riftpenaltypairs+counter*7+2)=el1;
					*(riftpenaltypairs+counter*7+3)=el2;
					*(riftpenaltypairs+counter*7+4)=normal[0];
					*(riftpenaltypairs+counter*7+5)=normal[1];
					*(riftpenaltypairs+counter*7+6)=length/2;
					counter++;
				}
				else{
					*(riftpenaltypairs+k2*7+4)+=normal[0];
					*(riftpenaltypairs+k2*7+5)+=normal[1];
					*(riftpenaltypairs+k2*7+6)+=length/2;
				}
			}
		}
		/*Renormalize normals: */
		for(j=0;j<counter;j++){
			double magnitude=sqrt(pow( double(riftpenaltypairs[j*7+4]),2) + pow( double(riftpenaltypairs[j*7+5]),2) );
			*(riftpenaltypairs+j*7+4)=*(riftpenaltypairs+j*7+4)/magnitude;
			*(riftpenaltypairs+j*7+5)=*(riftpenaltypairs+j*7+5)/magnitude;
		}

		riftspenaltypairs[i]=riftpenaltypairs;
		riftsnumpenaltypairs[i]=(numsegs/2-1);
	}

	/*Assign output pointers: */
	*priftspenaltypairs=riftspenaltypairs;
	*priftsnumpenaltypairs=riftsnumpenaltypairs;
	return noerr;
}/*}}}*/
int RemoveCornersFromRifts(int** pindex,int* pnel,double** px,double** py,int* pnods,int* segments,int* segmentmarkers,int num_seg){/*{{{*/

	int noerr=1;
	int i,j,k;
	int node1,node2,node3;
	int el;
	double  pair[2];
	int     pair_count=0;
	int     triple=0;

	/*Recover input: */
	int    *index = *pindex;
	int     nel   = *pnel;
	double *x     = *px;
	double *y     = *py;
	int     nods  = *pnods;

	for (i=0;i<num_seg;i++){
		node1=*(segments+3*i+0);
		node2=*(segments+3*i+1);
		/*Find all elements connected to [node1 node2]: */
		pair_count=0;
		for (j=0;j<nel;j++){
			if (*(index+3*j+0)==node1){
				if ((*(index+3*j+1)==node2) || (*(index+3*j+2)==node2)){
					pair[pair_count]=j;
					pair_count++;
				}
			}
			if (*(index+3*j+1)==node1){
				if ((*(index+3*j+0)==node2) || (*(index+3*j+2)==node2)){
					pair[pair_count]=j;
					pair_count++;
				}
			}
			if (*(index+3*j+2)==node1){
				if ((*(index+3*j+0)==node2) || (*(index+3*j+1)==node2)){
					pair[pair_count]=j;
					pair_count++;
				}
			}
		}
		/*Ok, we have pair_count elements connected to this segment. For each of these elements, 
		 *figure out if the third node also belongs to a segment: */
		if ((pair_count==0) || (pair_count==1)){ //we only select the rift segments, which belong to  2 elements
			continue;
		}
		else{
			for (j=0;j<pair_count;j++){
				el=(int)pair[j];
				triple=0;
				/*First find node3: */
				if (*(index+3*el+0)==node1){
					if (*(index+3*el+1)==node2)node3=*(index+3*el+2);
					else node3=*(index+3*el+1);
				}
				if (*(index+3*el+1)==node1){
					if (*(index+3*el+0)==node2)node3=*(index+3*el+2);
					else node3=*(index+3*el+0);
				}
				if (*(index+3*el+2)==node1){
					if (*(index+3*el+0)==node2)node3=*(index+3*el+1);
					else node3=*(index+3*el+0);
				}
				/*Ok, we have node3. Does node3 belong to a segment? : */
				for (k=0;k<num_seg;k++){
					if ((node3==*(segments+3*k+0)) || (node3==*(segments+3*k+1))){
						triple=1;
						break;
					}
				}
				if(triple==1){
					/*el is a corner element: we need to split it in 3 triangles: */
					x=xReNew<double>(x,nods,nods+1);
					y=xReNew<double>(y,nods,nods+1);
					x[nods]=(x[(int)node1-1]+x[(int)node2-1]+x[(int)node3-1])/3;
					y[nods]=(y[(int)node1-1]+y[(int)node2-1]+y[(int)node3-1])/3;
					index=xReNew<int>(index,nel*3,(nel+2*3));
					/*First, reassign element el: */
					*(index+3*el+0)=node1;
					*(index+3*el+1)=node2;
					*(index+3*el+2)=nods+1;
					/*Other two elements: */
					*(index+3*nel+0)=node2;
					*(index+3*nel+1)=node3;
					*(index+3*nel+2)=nods+1;

					*(index+3*(nel+1)+0)=node3;
					*(index+3*(nel+1)+1)=node1;
					*(index+3*(nel+1)+2)=nods+1;
					/*we need  to change the segment elements corresponding to el: */
					for (k=0;k<num_seg;k++){
						if (*(segments+3*k+2)==(el+1)){
							if ( ((*(segments+3*k+0)==node1) && (*(segments+3*k+1)==node2)) || ((*(segments+3*k+0)==node2) && (*(segments+3*k+1)==node1))) *(segments+3*k+2)=el+1;
							if ( ((*(segments+3*k+0)==node2) && (*(segments+3*k+1)==node3)) || ((*(segments+3*k+0)==node3) && (*(segments+3*k+1)==node2))) *(segments+3*k+2)=nel+1;
							if ( ((*(segments+3*k+0)==node3) && (*(segments+3*k+1)==node1)) || ((*(segments+3*k+0)==node1) && (*(segments+3*k+1)==node3))) *(segments+3*k+2)=nel+2;
						}
					}

					nods=nods+1;
					nel=nel+2;
					i=0;
					break;
				}
			} //for (j=0;j<pair_count;j++)
		}
	}// for (i=0;i<num_seg;i++)

	/*Assign output pointers: */
	*pindex=index;
	*pnel=nel;
	*px=x;
	*py=y;
	*pnods=nods;
	return noerr;
}/*}}}*/
