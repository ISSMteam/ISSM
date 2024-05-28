/*!\file MmeToInputFromId
 * \brief: compute damage
 */
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../classes/Inputs/DatasetInput.h"
#include "../../classes/Inputs/TransientInput.h"
#include "../../classes/Inputs/TriaInput.h"
#include "./MmeToInputFromIdx.h"

void MmeToInputFromIdx(Inputs* inputs, Elements* elements, int id, int rootenum, int interpolationenum){

	TransientInput* transientinput  = NULL;
	TransientInput* transientinput2 = NULL;
	Tria* element                    = NULL;
	IssmDouble value;
	IssmDouble* values               = NULL;
	IssmDouble* times                = NULL;
	int N;

	/*find thickness dataset: */
	DatasetInput* datasetinput = inputs->GetDatasetInput(rootenum);

	/*Initialize new transient input: */
	transientinput = datasetinput->GetTransientInputByOffset(0); _assert_(transientinput);
	transientinput->GetAllTimes(&times,&N);
	inputs->SetTransientInput(DummyEnum,times,N);
	transientinput2 = inputs->GetTransientInput(DummyEnum);

	for(Object* & object : elements->objects){
		Tria*   element=xDynamicCast<Tria*>(object);

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
			}
		}
	}

	/*wipe out existing SurfaceloadIceThicknessChangeEnum dataset:*/
	inputs->ChangeEnum(DummyEnum,rootenum);

}
