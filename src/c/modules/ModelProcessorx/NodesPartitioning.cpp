/*!\file:  NodesPartitioning.cpp
 * \brief: partition elements and nodes and vertices
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx/ModelProcessorx.h"

void  DiscontinuousGalerkinNodesPartitioning(bool** pmy_nodes,bool* my_elements,bool* my_vertices, IoModel* iomodel){

	/* Each element has it own nodes (as many as vertices) + additional nodes
	 * from neighboring elements for each face. This yields to a very different
	 * partition for the nodes and the vertices. The vertices are similar to
	 * continuous galerkin, but the nodes partitioning involves faces, which
	 * messes up sorting of ids. */

	/*Intermediaries*/
	int  i,i1,i2;
	int  e1,e2;
	int  pos;

	/*Get faces and elements*/
	CreateEdges(iomodel);

	/*Build discontinuous node partitioning
	 *  - there are three nodes per element (discontinous)
	 *  - for each element present of each partition, its three nodes will be in this partition
	 *  - the faces require the dofs of the 2 nodes of each elements sharing the face.
	 *    if the 2 elements sharing the face are on 2 different cpus, we must duplicate
	 *    the two nodes that are not on the cpus so that the face can access the dofs of
	 *    all its 4 nodes
	 */

	/*Allocate*/
	bool* my_nodes=xNewZeroInit<bool>(3*iomodel->numberofelements);

	/*First: add all the nodes of all the elements belonging to this cpu*/
	if(iomodel->domaintype==Domain2DhorizontalEnum || iomodel->domaintype==Domain2DverticalEnum){
		for (i=0;i<iomodel->numberofelements;i++){
			if (my_elements[i]){
				my_nodes[3*i+0]=true;
				my_nodes[3*i+1]=true;
				my_nodes[3*i+2]=true;
			}
		}
	}
	else{
		_error_("not implemented yet");
	}

	/*Second: add all missing nodes*/

	/*Get faces and elements*/
	CreateFaces(iomodel);

	if(iomodel->domaintype==Domain2DhorizontalEnum){
		/*!All elements have been partitioned above, only create elements for this CPU: */
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			e1=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]
			e2=iomodel->faces[4*i+3]-1; //faces are [node1 node2 elem1 elem2]

			/* 1) If the element e1 is in the current partition
			 * 2) and if the face of the element is shared by another element (internal face)
			 * 3) and if this element is not in the same partition:
			 * we must clone the nodes on this partition so that the loads (Numericalflux)
			 * will have access to their properties (dofs,...)*/
			if(my_elements[e1] && e2!=-2 && !my_elements[e2]){

				/*1: Get vertices ids*/
				i1=iomodel->faces[4*i+0];
				i2=iomodel->faces[4*i+1];

				/*2: Get the column where these ids are located in the index*/
				pos=UNDEF;
				for(int j=0;j<3;j++){
					if(iomodel->elements[3*e2+j]==i1) pos=j;
				}

				/*3: We have the id of the elements and the position of the vertices in the index
				 * we can now create the corresponding nodes:*/
				if(pos==0){
					my_nodes[e2*3+0]=true;
					my_nodes[e2*3+2]=true;
				}
				else if(pos==1){
					my_nodes[e2*3+1]=true;
					my_nodes[e2*3+0]=true;
				}
				else if(pos==2){
					my_nodes[e2*3+2]=true;
					my_nodes[e2*3+1]=true;
				}
				else{
					_error_("Problem in faces creation");
				}
			}
		}
	}

	/*Free data and assign output pointers */
	*pmy_nodes=my_nodes;
}
