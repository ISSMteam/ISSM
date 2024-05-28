/*!\file: ResetBoundaryConditions.cpp
 * \brief: change boundary conditions of a model, using a solution vector from another analysis
 */ 

#include "../classes/classes.h"
#include "../modules/modules.h"
#include "../shared/io/io.h"

void ResetBoundaryConditions(FemModel* femmodel, int analysis_type){

	/*variables: */
	Vector<IssmDouble>* yg = NULL;

	if(VerboseSolution()) _printf0_("   updating boundary conditions...\n");
	_assert_(femmodel->analysis_type_list[femmodel->analysis_counter]==analysis_type); 

	/*set current analysis: */
	femmodel->SetCurrentConfiguration(analysis_type);
	int index = femmodel->AnalysisIndex(analysis_type);

	/*retrieve boundary conditions from element inputs :*/
	GetSolutionFromInputsx(&yg,femmodel);

	/*update spcs using this new vector of constraints: */
	UpdateDynamicConstraintsx(femmodel->constraints,femmodel->nodes,femmodel->parameters,yg);

	/*Free resources:*/
	delete yg;
}
