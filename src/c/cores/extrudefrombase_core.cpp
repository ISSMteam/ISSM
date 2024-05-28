/*!\file: extrudefrombase_core.cpp
 * \brief: core of the extrusion solution
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void extrudefrombase_core(FemModel* femmodel){

	/*Intermediaries*/
	int elementtype,domaintype;

	/*Get parameters*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&elementtype,MeshElementtypeEnum);

	/*If this is a 2D horizontal domain, return (no need to extrude)*/
	if(domaintype==Domain2DhorizontalEnum) return;
	if(domaintype==Domain3DsurfaceEnum) return;

	int inputenum;
	femmodel->parameters->FindParam(&inputenum,InputToExtrudeEnum);
	if(VerboseSolution()) _printf0_("   extruding "<<EnumToStringx(inputenum)<<" from base...\n");

	/*Special method for Penta, otherwise call solution sequence*/
	if(elementtype==PentaEnum){
		InputExtrudex(femmodel,inputenum,-1);
	}
	else{
		/*Call on core computations: */
		femmodel->SetCurrentConfiguration(ExtrudeFromBaseAnalysisEnum);
		femmodel->UpdateConstraintsExtrudeFromBasex();
		solutionsequence_linear(femmodel);
	}
}
