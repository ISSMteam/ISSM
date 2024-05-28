/*!\file:  ElementsAndVerticesPartitioning.cpp
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

void  ElementsAndVerticesPartitioning(IoModel* iomodel){

	int numberofelements2d;
	int numberofvertices2d;
	int numlayers;

	/*intermediary: */
	int *epart          = NULL; //element partitioning.
	int *npart          = NULL; //node partitioning.
	int  elements_width;        //number of columns in elements (2d->3, 3d->6)
	int *elements2d     = NULL;

	/*Get my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*First, check that partitioning has not yet been carryed out. Just check whether my_elements pointers is not already assigned a value: */
	if(iomodel->my_elements) return;

	/*Number of vertices per elements, needed to correctly retrieve data: */
	/*Determine parallel partitioning of elements: we use Metis for now. First load the data, then partition*/
	switch(iomodel->meshelementtype){
		case TriaEnum:
			elements_width=3;
			numberofelements2d = 0;
			numberofvertices2d = 0;
			numlayers          = 0;
			break;
		case TetraEnum:
			elements_width=4;
			numberofelements2d = 0;
			numberofvertices2d = 0;
			numlayers          = 0;
			break;
		case PentaEnum:
			elements_width=6;
			iomodel->FetchData(&elements2d,NULL,NULL,"md.mesh.elements2d");
			iomodel->FindConstant(&numberofelements2d,"md.mesh.numberofelements2d");
			iomodel->FindConstant(&numberofvertices2d,"md.mesh.numberofvertices2d");
			iomodel->FindConstant(&numlayers,"md.mesh.numberoflayers");
			break;
		default:
			_error_("mesh elements "<< EnumToStringx(iomodel->meshelementtype) <<" not supported yet");
	}

	/*Use ice levelset for weights*/
	int fordan = 0;
	int* weights = NULL;
	if(fordan){
		IssmDouble* icelevelset = NULL;
		iomodel->FetchData(&icelevelset,NULL,NULL,"md.mask.ice_levelset");

		weights = xNew<int>(iomodel->numberofvertices);
		for(int i=0;i<iomodel->numberofvertices;i++){
			if(icelevelset[i]>=0) weights[i] = 1;
			if(icelevelset[i]<0)  weights[i] = 100;
		}
		xDelete<IssmDouble>(icelevelset);
	}

	/*Partition and free resouces*/
	MeshPartitionx(&epart,&npart,iomodel->numberofelements,iomodel->numberofvertices,iomodel->elements,numberofelements2d,numberofvertices2d,elements2d,weights,numlayers,elements_width,iomodel->meshelementtype,num_procs);
	xDelete<int>(elements2d);
	xDelete<int>(npart);
	xDelete<int>(weights);

	if(fordan){
		for(int i=0;i<IssmComm::GetSize();i++){
			if(i==IssmComm::GetRank()){
				int temp =0;
				for(int j=0;j<iomodel->numberofelements;j++) if(epart[j]==i) temp++;
				_printf_("Partition #"<<i<<" number of elements: "<<temp<<"\n");
			}
			ISSM_MPI_Barrier(IssmComm::GetComm());
		}
	}

	/*Deal with rifts, they have to be included into one partition only, not several: */
	int numrifts;
	iomodel->FindConstant(&numrifts,"md.rifts.numrifts");
	if(numrifts){
		IssmDouble *riftinfo = NULL;
		iomodel->FetchData(&riftinfo,&numrifts,NULL,"md.rifts.riftstruct");
		for(int i=0;i<numrifts;i++){
			const int RIFTINFOSIZE = 12;
			int el1=reCast<int>(*(riftinfo+RIFTINFOSIZE*i+2))-1; //matlab indexing to c indexing
			int el2=reCast<int>(*(riftinfo+RIFTINFOSIZE*i+3))-1; //matlab indexing to c indexing
			epart[el2]=epart[el1]; //ensures that this pair of elements will be in the same partition, as well as the corresponding vertices;
		}
		iomodel->DeleteData(riftinfo,"md.rifts.riftstruct");
	}

	/*Create my_elements, used by each partition */
	bool *my_elements = xNewZeroInit<bool>(iomodel->numberofelements);

	/*Start figuring out, out of the partition, which elements belong to this cpu: */
	bool check = false;
	for(int i=0;i<iomodel->numberofelements;i++){

		/*!All elements have been partitioned above, only deal with elements for this cpu: */
		if(my_rank==epart[i]){
			my_elements[i]=true;
			check = true;
		}
	}
	if(!check) _error_("partition "<<my_rank<<" does not have any element! Try reducing md.cluster.np");

	/*Assign pointers to iomodel*/
	iomodel->epart      =epart;
	iomodel->my_elements=my_elements;
}
