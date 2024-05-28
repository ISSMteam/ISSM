/*!\file ThicknessAbsMisfitx
 * \brief: compute misfit between observations and model
 */

#include "./ThicknessAbsMisfitx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/DatasetInput.h"

void ThicknessAbsMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		J+=ThicknessAbsMisfit(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble ThicknessAbsMisfit(Element* element){

	IssmDouble  thickness,thicknessobs,weight;
	IssmDouble  Jelem=0.;
	IssmDouble  misfit,Jdet;
	IssmDouble* xyz_list = NULL;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/* Get node coordinates*/
	element->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	DatasetInput* weights_input     =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* thickness_input   =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);
	Input* thicknessobs_input=element->GetInput(InversionThicknessObsEnum);              _assert_(thicknessobs_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,ThicknessAbsMisfitEnum);
		thickness_input->GetInputValue(&thickness,gauss);
		thicknessobs_input->GetInputValue(&thicknessobs,gauss);

		/*Compute ThicknessAbsMisfitEnum*/
		Jelem+=0.5*(thickness-thicknessobs)*(thickness-thicknessobs)*weight*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
