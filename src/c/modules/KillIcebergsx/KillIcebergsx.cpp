/*!\file KillIcebergsx
 * \brief: compute inverse method gradient
 */

#include "./KillIcebergsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

#include "../InputUpdateFromVectorx/InputUpdateFromVectorx.h"

int KillIcebergsx(FemModel* femmodel){

	/*Intermediaries*/
	int lid;

	/*retrieve vertex info and prepare element flag to speed up process*/
	int         nbv_local    = femmodel->vertices->Size();
	IssmDouble *local_mask   = xNewZeroInit<IssmDouble>(nbv_local);
	bool       *element_flag = xNewZeroInit<bool>(femmodel->elements->Size());
	IssmDouble ice;

	/*Step 1, go through all elements and put 1 in local_mask if the element is grounded*/
	int i=0;
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		if(!element->IsIceInElement()){
			/*Nothing to do, just flag element to speed up the computation*/
			element_flag[i] = true;
		}
		else {
			if(element->IsAllGrounded()){
				/* only look at element with ice but not fully grounded */
				int numvertices = element->GetNumberOfVertices();
				Gauss* gauss=element->NewGauss();
				Input* icelevel_input = element->GetInput(MaskIceLevelsetEnum);				_assert_(icelevel_input);

				for(int v=0;v<numvertices;v++) {
					gauss->GaussVertex(v);
					icelevel_input->GetInputValue(&ice,gauss);
					/* The initial mask is very strict, we look at all grounded elements and set the mask for ice nodes only. */
					if (ice < 0) local_mask[element->vertices[v]->Lid()] = 1.;
				}
				delete gauss;
			}
		}
		i++;
	}

	/*Now we have 2 loops, one across cpus, and one for each cpus: we are going
	 * to propagate the mask if an element is connected to a positive mask
	 * already.  We then communicate to the other partitions. We stop when the
	 * mask stops changing*/
	bool keepsyncing = true;
	while(keepsyncing){

		/*Get local mask from parallel vector*/
		femmodel->SyncLocalVectorWithClonesVerticesAdd(local_mask);

		/*Local iterations on partition*/
		bool keepgoing    = true;
		int  iter         = 1;
		while(keepgoing){
			//_printf0_("   -- Kill icebergs: local iteration "<<iter<<"\n");

			keepgoing    = false;
			int i=0;
			for(Object* & object : femmodel->elements->objects){
				Element* element=xDynamicCast<Element*>(object);

				if(!element_flag[i]){
					int numvertices = element->GetNumberOfVertices();
					bool found1 = false;
					IssmDouble sumlocalmask = 0.;

					for(int j=0;j<numvertices;j++){
						lid = element->vertices[j]->Lid();
						/*we need to have at least a sharing edge, to extend the mask*/
						sumlocalmask += local_mask[lid];
						if(sumlocalmask > 1.5){
							found1 = true;
							break;
						}
					}
					if(found1){
						element_flag[i] = true;
						for(int j=0;j<numvertices;j++){
							lid = element->vertices[j]->Lid();
							local_mask[lid]=1.;
						}
						keepgoing = true;
					}
				}
				i++;
			}
			iter++;
		}

		/*Check how many iterations all cpus did*/
		int iter_max;
		ISSM_MPI_Reduce(&iter,&iter_max,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(&iter_max,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		if(iter_max==2){
			/*If iter is only 2, nothing else was changed in the while loop above (iter is initialized as 1 and then ++)*/
			keepsyncing = false;
		}
	}

	/*Cleanup*/
	xDelete<bool>(element_flag);

	int killbergReinit = 0;
	/*OK, now deactivate iceberg and count the number of deactivated vertices*/
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		if(element->IsIceInElement()){
			int  numvertices = element->GetNumberOfVertices();
			bool deactivate = false;
			for(int j=0;j<numvertices;j++){
				lid = element->vertices[j]->Lid();
				if(local_mask[lid]==0.){
					deactivate = true;
					break;
				}
			}

			if(deactivate){
				int  numvertices = element->GetNumberOfVertices();
				IssmDouble* values = xNew<IssmDouble>(numvertices);
				for(int j=0;j<numvertices;j++) values[j] = 1.; /*Anything >0 = no ice*/
				element->AddInput(MaskIceLevelsetEnum,values,P1Enum);
				xDelete<IssmDouble>(values);
				killbergReinit += 1;
			}
		}
	}
	/*cleanup*/
	xDelete<IssmDouble>(local_mask);

	/*Recompute the sign distance for the levelset function*/
	return killbergReinit;
}
