/*!\file CreateJacobianMatrixx
 * \brief retrieve vector from inputs in elements
 */

#include "./CreateJacobianMatrixx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../AllocateSystemMatricesx/AllocateSystemMatricesx.h"

void CreateJacobianMatrixx(Matrix<IssmDouble>** pJff,FemModel* femmodel,IssmDouble kmax){

	int      configuration_type,analysisenum;
	Element *element = NULL;
	Load    *load    = NULL;
	Matrix<IssmDouble>* Jff = NULL;

	/*Checks*/
	_assert_(femmodel && femmodel->nodes && femmodel->elements);

	/*Recover some parameters*/
	femmodel->parameters->FindParam(&configuration_type,ConfigurationTypeEnum);
	femmodel->parameters->FindParam(&analysisenum,AnalysisTypeEnum);
	Analysis* analysis = EnumToAnalysis(analysisenum);

	/*Initialize Jacobian Matrix*/
	AllocateSystemMatricesx(&Jff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		ElementMatrix* Je = analysis->CreateJacobianMatrix(element);
		if(Je) Je->AddToGlobal(Jff);
		delete Je;
	}
	for (Object* & object : femmodel->loads->objects){
		load=(Load*)object;
		load->CreateJacobianMatrix(Jff);
		load->PenaltyCreateJacobianMatrix(Jff,kmax);
	}
	Jff->Assemble();

	/*Assign output pointer*/
	delete analysis;
	*pJff=Jff;

}
