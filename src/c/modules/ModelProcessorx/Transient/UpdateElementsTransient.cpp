/*
 * UpdateElementsTransient:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void	UpdateElementsTransient(Elements* elements, Parameters* parameters,Inputs* inputs,IoModel* iomodel){

	/*FIXME: this should go into parameterization update*/

	bool isgroundingline;
	parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);

	if(isgroundingline){
		iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
	}
}
