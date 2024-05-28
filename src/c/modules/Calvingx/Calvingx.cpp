/*!\file Calvingx
 * \brief: compute inverse method gradient
 */

#include "./Calvingx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void Calvingx(FemModel* femmodel){

	/*Recover Calving law Enum*/
	int calvinglaw;
	femmodel->parameters->FindParam(&calvinglaw,CalvingLawEnum);

	/*Calculate calving rate*/
	switch(calvinglaw){
		case CalvingMinthicknessEnum:
		case CalvingHabEnum:
			femmodel->ElementOperationx(&Element::CalvingSetZeroRate);
			break;
		case DefaultCalvingEnum:
			femmodel->ElementOperationx(&Element::CalvingRateToVector);
			break;
		case CalvingCrevasseDepthEnum:
			femmodel->ElementOperationx(&Element::CalvingSetZeroRate);
			/*rate is 0 but we need to calculate a few things to determine where it will calve*/
			femmodel->StrainRateparallelx();
			femmodel->StrainRateeffectivex();
			femmodel->DeviatoricStressx();
			femmodel->ElementOperationx(&Element::CalvingCrevasseDepth);
			break;
		case CalvingLevermannEnum:
			femmodel->StrainRateparallelx();
			femmodel->StrainRateperpendicularx();
			femmodel->CalvingRateLevermannx();
			break;
		case CalvingVonmisesEnum:
		case CalvingDev2Enum:
			femmodel->ElementOperationx(&Element::CalvingRateVonmises);
			break;
		case CalvingVonmisesADEnum:
			femmodel->ElementOperationx(&Element::CalvingRateVonmisesAD);
			break;
		case CalvingTestEnum:
			femmodel->ElementOperationx(&Element::CalvingRateTest);
			break;
		case CalvingParameterizationEnum:
			femmodel->ElementOperationx(&Element::CalvingRateParameterization);
			break;
		case CalvingPollardEnum:
			femmodel->ElementOperationx(&Element::CalvingPollard);
			break;
		case CalvingCalvingMIPEnum:
			femmodel->ElementOperationx(&Element::CalvingRateCalvingMIP);
			break;
		default:
			_error_("Caving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}
}

void MovingFrontalVelx(FemModel* femmodel){
	femmodel->ElementOperationx(&Element::MovingFrontalVelocity);
}
