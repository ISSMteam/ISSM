/*!\file:  AssociateSegmentToElement.cpp
 * \brief for each segment, look for the corresponding element.
 */ 

#include "./triangle.h"

int AssociateSegmentToElement(int** psegments,int nseg,int* index,int nel){

	/*node indices: */
	int A,B;

	/*Recover segments: */
	int* segments=*psegments;

	for(int i=0;i<nseg;i++){
		A=segments[3*i+0];
		B=segments[3*i+1];
		segments[3*i+2]=FindElement(A,B,index,nel)+1; //matlab indexing.
	}

	/*Assign output pointers: */
	*psegments=segments;
	return 1;
}
