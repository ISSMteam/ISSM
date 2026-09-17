/*!\file Calvingx
 * \brief: compute inverse method gradient
 */

#include "./Calvingx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include <random>

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
		case CalvingStochasticEnum:
			CalvingStochasticx(femmodel);
			break;
		default:
			_error_("Caving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}
}

void MovingFrontalVelx(FemModel* femmodel){
	femmodel->ElementOperationx(&Element::MovingFrontalVelocity);
}
void CalvingStochasticx(FemModel* femmodel){

	/*Intermediaries*/
	IssmPDouble r;
	IssmDouble  delta_t,f,chi_max,k;

	/*1. generate random number*/
	static std::mt19937 gen(std::random_device{}());          // or gen(1234) for reproducibility
	static std::uniform_real_distribution<double> dis(0.0,1.0); // [0,1)
	r = dis(gen);

	/*2. Calculate random probability*/
	femmodel->parameters->FindParam(&delta_t, TimesteppingTimeStepEnum);
	femmodel->parameters->FindParam(&f, CalvingFEnum);
	femmodel->parameters->FindParam(&k, CalvingKEnum);
	femmodel->parameters->FindParam(&chi_max, CalvingChiMaxEnum);
	IssmDouble Pmax = 1. - exp(-f*delta_t);

	/*3. cap random probability*/
	IssmDouble P = min(r, Pmax);

	/*4. Stochastic waitime (charactiristic time)*/
	IssmDouble tau = -(f*delta_t)/log(1. - P);

	/*5. Define delta critical */
	IssmDouble chi_crit = chi_max - log(tau)/k;
	femmodel->parameters->SetParam(chi_crit, CalvingChiCritEnum);

	/*Loop over elements and compute crevasse depth*/
	femmodel->DeviatoricStressx();
	femmodel->ElementOperationx(&Element::CalvingCrevasseDepth);
}
