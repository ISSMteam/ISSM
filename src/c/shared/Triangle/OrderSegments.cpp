/*
 * OrderSegments.c: 
 * reorder segments so that their normals point outside the domain outline.
 */
#include "./triangle.h"

int OrderSegments(int** psegments,int nseg,int* index,int nel){

	/*vertex indices: */
	int A,B;

	/*element index*/
	int el;

	/*Recover segments: */
	int* segments=*psegments;

	for(int i=0;i<nseg;i++){
		A=segments[3*i+0];
		B=segments[3*i+1];
		el=segments[3*i+2]-1; //after AssociateSegmentToElement, el was a matlab index, we need the c index now.

		if (index[3*el+0]==A){
			if (index[3*el+2]==B){
				segments[3*i+0]=B;
				segments[3*i+1]=A;
			}
		}
		else if (index[3*el+1]==A){
			if (index[3*el+0]==B){
				segments[3*i+0]=B;
				segments[3*i+1]=A;
			}
		}
		else{
			if (index[3*el+1]==B){
				segments[3*i+0]=B;
				segments[3*i+1]=A;
			}
		}
	}

	/*Assign output pointers: */
	*psegments=segments;
	return 1;
}
