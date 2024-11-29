/*!\file MmeToInputx
 * \brief: 
 */

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/TransientInput.h"
#include "../../classes/Inputs/DatasetInput.h"
#include "../../classes/Inputs/TriaInput.h"
#include "./MmeToInputx.h"

void  MmeToInputx(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* variable_partition,int npart,int rootenum, int interpolationenum){ 

	TransientInput* transientinput  = NULL;
	TransientInput* transientinput2 = NULL;
	Tria* element                    = NULL;
	IssmDouble value;
	IssmDouble* values               = NULL;
	IssmDouble* times                = NULL;
	int N;
	int id;

	/*find thickness dataset: */
	DatasetInput* datasetinput = femmodel->inputs->GetDatasetInput(rootenum);

	/*Initialize new transient input: */
	transientinput = datasetinput->GetTransientInputByOffset(0); _assert_(transientinput);
	transientinput->GetAllTimes(&times,&N);
	femmodel->inputs->SetTransientInput(DummyEnum,times,N);
	transientinput2 = femmodel->inputs->GetTransientInput(DummyEnum); transientinput2->Configure(femmodel->parameters);

	for(Object* & object : femmodel->elements->objects){
		Tria*   element=xDynamicCast<Tria*>(object);

		if(reCast<int>(variable_partition[element->Sid()])==-1)id=0; //grab background field
		else id=distributed_values[reCast<int>(variable_partition[element->Sid()])]-1; //grab partition field

		/*recover the right field from the mme: */
		transientinput = datasetinput->GetTransientInputByOffset(id); _assert_(transientinput);

		/*copy values from the transientinput to the final transientinput2: */
		for (int j=0;j<N;j++){
			TriaInput* tria_input=transientinput->GetTriaInput(j);
			element->InputServe(tria_input);
			if(interpolationenum==P0Enum){
				value=tria_input->element_values[0];
				transientinput2->AddTriaTimeInput( j,1,&(element->lid),&value,P0Enum);
			}
			else if(interpolationenum==P1Enum){

				/*Get values and lid list*/
				const int   numvertices     = element->GetNumberOfVertices();
				int        *vertexlids      = xNew<int>(numvertices);
				int        *vertexsids      = xNew<int>(numvertices);

				/*Recover vertices ids needed to initialize inputs*/
				element->GetVerticesLidList(&vertexlids[0]);
				element->GetVerticesSidList(&vertexsids[0]);
				values=tria_input->element_values;
				transientinput2->AddTriaTimeInput( j,numvertices,vertexlids,values,P1Enum);

				/*free memory: */
				xDelete<int>(vertexlids);
				xDelete<int>(vertexsids);
			}
		}
	}

	/*wipe out existing SurfaceloadIceThicknessRateEnum dataset:*/
	femmodel->inputs->ChangeEnum(DummyEnum,rootenum);

	/*free memory:*/
	xDelete<IssmDouble>(times);

}	
