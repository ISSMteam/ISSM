/*!\file: depthaverage_core.cpp
 * \brief: core of the extrusion solution
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void depthaverage_core(FemModel* femmodel){

	/*Intermediaries*/
	int domaintype,elementtype;
	int inputenum,input_average_enum;

	/*Get parameters*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&elementtype,MeshElementtypeEnum);
	femmodel->parameters->FindParam(&inputenum,InputToDepthaverageInEnum);
	femmodel->parameters->FindParam(&input_average_enum,InputToDepthaverageOutEnum);

	/*If this is a 2D horizontal domain: no need to do anything, just copy input*/
	if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DsurfaceEnum){
		InputDuplicatex(femmodel,inputenum,input_average_enum);
		return;
	}

	if(VerboseSolution()) _printf0_("   depth averaging "<<EnumToStringx(inputenum)<<"\n");

	/*Special method for Penta, otherwise call solution sequence*/
	if(elementtype==PentaEnum){
		InputDepthAverageAtBasex(femmodel,inputenum,input_average_enum);
	}
	else{
		/*Call on core computations: */
		femmodel->SetCurrentConfiguration(DepthAverageAnalysisEnum);
		solutionsequence_linear(femmodel);
	}
}
