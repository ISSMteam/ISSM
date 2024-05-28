/*!\file GetSolutionFromInputsx
 * \brief: update datasets using  parameter inputs
 */

#include "./GetSolutionFromInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetSolutionFromInputsx(Vector<IssmDouble>** psolution,FemModel* femmodel){/*{{{*/

	if(VerboseModule()) _printf0_("   Get solution from inputs\n");

	/*retrieve parameters: */
	int analysisenum;
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);

	/*Get size of vector: */
	int gsize       = femmodel->nodes->NumberOfDofs(GsetEnum);
	int gsize_local = femmodel->nodes->NumberOfDofsLocal(GsetEnum);
	if(gsize==0) _error_("Allocating a Vec of size 0 as gsize=0 ");

	/*Initialize solution: */
	Vector<IssmDouble>* solution=new Vector<IssmDouble>(gsize_local,gsize);

	/*Go through elements and plug solution: */
	Analysis* analysis = EnumToAnalysis(analysisenum);
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		analysis->GetSolutionFromInputs(solution,element);
	}
	delete analysis;

	/*Assemble vector: */
	solution->Assemble();

	/*Assign output pointers:*/
	*psolution=solution;
}/*}}}*/
void GetBasalSolutionFromInputsx(Vector<IssmDouble>** psolution,FemModel* femmodel){ /*{{{*/

	if(VerboseModule()) _printf0_("   Get solution from inputs\n");

	/*retrieve parameters: */
	int analysisenum;
	int domaintype;
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	/*Get size of vector: */
	int gsize       = femmodel->nodes->NumberOfDofs(GsetEnum);
	int gsize_local = femmodel->nodes->NumberOfDofsLocal(GsetEnum);
	if(gsize==0) _error_("Allocating a Vec of size 0 as gsize=0 ");

	/*Initialize solution: */
	Vector<IssmDouble>* solution=new Vector<IssmDouble>(gsize_local,gsize);

	/*Go through elements and plug solution: */
	Analysis* analysis = EnumToAnalysis(analysisenum);
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		switch(domaintype){
			case Domain2DhorizontalEnum:
				analysis->GetSolutionFromInputs(solution,element);
			break;
		case Domain3DEnum:
			if(!element->IsOnBase()) continue;
			analysis->GetSolutionFromInputs(solution,element);
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
		}
	}
	delete analysis;

	/*Assemble vector: */
	solution->Assemble();

	/*Assign output pointers:*/
	*psolution=solution;
}/*}}}*/
