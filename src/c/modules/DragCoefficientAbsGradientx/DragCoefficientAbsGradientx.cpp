/*!\file DragCoefficientAbsGradientx
 * \brief: compute misfit between observations and model
 */

#include "./DragCoefficientAbsGradientx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/DatasetInput.h"

void DragCoefficientAbsGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		J+=DragCoefficientAbsGradient(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble DragCoefficientAbsGradient(Element* element){

	int         domaintype,numcomponents;
	int frictionlaw;
	IssmDouble  Jelem=0.;
	IssmDouble  misfit,Jdet;
	IssmDouble  dp[2],weight;
	IssmDouble* xyz_list      = NULL;

	/*Get basal element*/
	if(!element->IsOnBase()) return 0.;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   numcomponents   = 1; break;
		case Domain3DEnum:           numcomponents   = 2; break;
		case Domain2DhorizontalEnum: numcomponents   = 2; break;
		default: _error_("not supported yet");
	}

	/*Spawn basal element*/
	Element* basalelement = element->SpawnBasalElement();

	/* Get node coordinates*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	DatasetInput* weights_input=basalelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum);   _assert_(weights_input);

	/* get the friction law: if 2-Weertman, 11-Schoof or 14-RegularizedCoulomb, which has a different names of C */
	element->FindParam(&frictionlaw, FrictionLawEnum);
	Input* drag_input = NULL;
	switch(frictionlaw) {
		case 2:
		case 11:
		case 13:
		case 14:
			drag_input = basalelement->GetInput(FrictionCEnum); _assert_(drag_input);
			break;
		default:
			drag_input = basalelement->GetInput(FrictionCoefficientEnum); _assert_(drag_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,DragCoefficientAbsGradientEnum);
		drag_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

		/*Compute Tikhonov regularization J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */
		Jelem+=weight*.5*dp[0]*dp[0]*Jdet*gauss->weight;
		if(numcomponents==2) Jelem+=weight*.5*dp[1]*dp[1]*Jdet*gauss->weight;

	}

	/*clean up and Return: */
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
