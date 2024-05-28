/*!\file InputUpdateFromSolutionx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromSolutionx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void InputUpdateFromSolutionx(FemModel* femmodel,Vector<IssmDouble>* solution){

	/*GetAnalysis*/
	int analysisenum;
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);
	Analysis* analysis = EnumToAnalysis(analysisenum);

	/*Display message*/
	if(VerboseModule()) _printf0_("   Updating inputs from solution for " << EnumToStringx(analysisenum) << "\n");

	/*Get local vector with both masters and slaves:*/
	IssmDouble *local_ug = NULL;
	femmodel->GetLocalVectorWithClonesGset(&local_ug,solution);

	/*Now update inputs (analysis specific)*/
	for(Object* & object : femmodel->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		analysis->InputUpdateFromSolution(local_ug,element);
	}

	/*cleanup and return*/
	delete analysis;
	xDelete<IssmDouble>(local_ug);
}
