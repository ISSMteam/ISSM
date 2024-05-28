/*!\file ElementConnectivityx
 * \brief: compute element connectivity table, using node connectivity table and elements.
 *
 * For each element, we want to know which neighbouring elements it connects to (fully, via an entire segment, not by a node).
 * We use the nodeconnectivity to speed up the computation. The nodeconnectivity gives us for each node of the element, 
 * all the neighbouring elements of this node, which are good candidates to be neighbours of the element itself.
 * For now, only triangular elements, ie 3 neighbours max per element.
 */

#include "./ElementConnectivityx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

int hascommondedge(int* element1,int* element2);

void ElementConnectivityx(int** pelementconnectivity,int* elements, int nels,int* nodeconnectivity, int nods, int width){

	int i,j,k,n;

	/*intermediary: */
	int    maxels;
	int  connectedelement;
	int    connectedelementindex;
	int    node;
	int    index;
	int    num_elements;

	/*maxels: */
	maxels=width-1;

	/*Allocate connectivity: */
	int* elementconnectivity=xNewZeroInit<int>(nels*3);

	/*Go through all elements, and for each element, go through its nodes, to get the neighbouring elements. 
	 * Once we get the neighbouring elements, figure out if they share a segment with the current element. If so, 
	 * plug them in the connectivity, unless they are already there.: */
	for(n=0;n<nels;n++){

		//element=n+1; //matlab indexing

		for(i=0;i<3;i++){

			node=elements[n*3+i]; //already matlab indexed, elements comes directly from the workspace.
			index=node-1;

			num_elements=nodeconnectivity[width*index+maxels]; //retrieve number of elements already  plugged into the connectivity of this node.

			for(j=0;j<num_elements;j++){

				/*for each element connected to node, figure out if it has a commond edge with element: */
				connectedelement=nodeconnectivity[width*index+j];
				connectedelementindex=connectedelement-1; //go from matlab indexing to c indexing.

				if(hascommondedge(&elements[n*3+0],&elements[connectedelementindex*3+0])){
					/*Ok, this connected element has a commond edge  with element, plug it into elementconnectivity, unless 
					 *it is already there: */

					for(k=0;k<3;k++){
						if(elementconnectivity[3*n+k]==0){
							elementconnectivity[3*n+k]=connectedelement;
							break;
						}
						else{
							if(connectedelement==elementconnectivity[3*n+k]) break;
						}
					}
				}
			}
		}
	}

	/*Assign output pointers: */
	*pelementconnectivity=elementconnectivity;
}

int hascommondedge(int* el1,int* el2){

	int count=0;
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(el1[i]==el2[j]) count++;
		}
	}
	if(count==2)
	 return 1;
	else
	 return 0;
}
