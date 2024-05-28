/*!\file NodeConnectivityx
 * \brief: compute node connectivity table, using elements connectivity table.
 *
 * For each node, we want to know how many elements are connected to this element, and which they are. 
 * Given that the 2d meshes we create in ISSM are triangular for now, and they are delaunay conforming, 
 * each triangle has a minimum angle of 30 degrees, which implies a connectivity <=6. We therefore return 
 * a nods x 7 connectivity table, with the 7'th column giving us the number of elements connected to each 
 * row node, and the first 6 columns giving us the elements numbers. 
 * Amend that: sounds like some triangles get up to 9 connectivity. Take 10 to be on the safe side.
 * In order to be compatible with matlab output, the connectivity table is given in matlab indexing (starts at 1).
 */

#include "./NodeConnectivityx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void	NodeConnectivityx(int** pconnectivity,int* pwidth,int* elements, int nels, int nods){

	int i,j,n;
	const int maxels=100;
	const int width=maxels+1;

	/*intermediary: */
	int     node;
	int     index;
	int     num_elements;
	int     already_plugged=0;
	int element;

	/*Allocate connectivity: */
	int* connectivity=xNewZeroInit<int>(nods*width);

	/*Go through all elements, and for each elements, plug into the connectivity, all the nodes. 
	 * If nodes are already plugged into the connectivity, skip them.: */
	for(n=0;n<nels;n++){

		element=n+1; //matlab indexing

		for(i=0;i<3;i++){

			node=elements[n*3+i]; //already matlab indexed, elements comes directly from the workspace.
			index=node-1;

			num_elements=connectivity[width*index+maxels]; //retrieve number of elements already  plugged into the connectivity of this node.

			already_plugged=0;
			for(j=0;j<num_elements;j++){
				if (element==*(connectivity+width*index+j)){
					already_plugged=1;
					break;
				}
			}
			if(already_plugged)break;

			/*this elements is not yet plugged  into the connectivity for this node, do it, and increase counter: */
			connectivity[width*index+num_elements]=element;
			connectivity[width*index+maxels]=num_elements+1;

		}
	}

	/*Last check: is the number of elements on last column of the connectivity superior to maxels? If so, then error out and 
	 * warn the user to increase the connectivity width: */
	for(i=0;i<nods;i++){
		if (*(connectivity+width*i+maxels)>maxels)
		 _error_("max connectivity width reached (" << *(connectivity+width*i+maxels) << ")! increase width of connectivity table");
	}

	/*Assign output pointers: */
	*pconnectivity=connectivity;
	*pwidth=width;
}
