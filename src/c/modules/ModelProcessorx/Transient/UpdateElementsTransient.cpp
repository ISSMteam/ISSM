/*
 * UpdateElementsTransient:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"

void	UpdateElementsTransient(Elements* elements, Parameters* parameters,Inputs* inputs,IoModel* iomodel){

	/*FIXME: this should go into parameterization update*/

	bool isgroundingline;
	parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);

	if(isgroundingline){
		iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
	}

	/*If we run with adaptive time step, we need to make sure that Vz is provided*/
	int timestepping_type;
	iomodel->FindConstant(&timestepping_type,"md.timestepping.type");
	if(timestepping_type==AdaptiveTimesteppingEnum && iomodel->domaintype==Domain3DEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vz",VzEnum,0.);
	}
}
