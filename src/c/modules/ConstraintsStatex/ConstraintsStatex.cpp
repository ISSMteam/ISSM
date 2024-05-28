/*!\file ConstraintsStatex
 * \brief: set up penalty constraints on loads
 */

#include "./ConstraintsStatex.h"
#include "./ConstraintsStateLocal.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ConstraintsStatex(int* pconverged, int* pnum_unstable_constraints,FemModel* femmodel){

	/*Early return if no rift and no penalties*/
	if(femmodel->loads->numrifts == 0 && femmodel->loads->numpenalties == 0){
		*pconverged                = 0;
		*pnum_unstable_constraints = 0;
		return;
	}

	/*output: */
	int converged                    = 1;
	int num_unstable_constraints     = 0;
	int min_mechanical_constraints   = 0;
	int unstable                     = 0;
	int sum_num_unstable_constraints = 0;
	int analysis_type;

	/*Display message*/
	if(VerboseModule()) _printf0_("   Constraining penalties\n");

	/*recover parameters: */
	femmodel->parameters->FindParam(&min_mechanical_constraints,StressbalanceRiftPenaltyThresholdEnum);
	femmodel->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Rift penalties first*/
	if(femmodel->loads->numrifts){
		RiftConstraintsState(&converged,&num_unstable_constraints,femmodel->loads,min_mechanical_constraints,analysis_type);
	}

	/*Deal with pengrid*/
	if(femmodel->loads->numpenalties){
		for(Object* & object : femmodel->loads->objects){
			Load* load=(Load*)object;
			if(load->ObjectEnum()==PengridEnum){
				Pengrid* pengrid=(Pengrid*)load;
				pengrid->ConstraintActivate(&unstable);
				num_unstable_constraints += unstable;
			}
		}
	}

	ISSM_MPI_Reduce(&num_unstable_constraints,&sum_num_unstable_constraints,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_num_unstable_constraints,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	num_unstable_constraints=sum_num_unstable_constraints;

	/*Have we converged? : */
	if(num_unstable_constraints) converged=0;

	/*Assign output pointers: */
	*pconverged                = converged;
	*pnum_unstable_constraints = num_unstable_constraints;
}
