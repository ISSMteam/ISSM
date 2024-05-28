/*!\file ThicknessAcrossGradientx
 * \brief: compute misfit between observations and model
 */

#include "./ThicknessAcrossGradientx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/DatasetInput.h"

void ThicknessAcrossGradientx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(Object* & object : elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		J+=ThicknessAcrossGradient(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble ThicknessAcrossGradient(Element* element){

	IssmDouble  thickness,thicknessobs,weight;
	IssmDouble  Jelem=0.;
	IssmDouble  misfit,Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];
	IssmDouble  vx,vy,vel;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/* Get node coordinates*/
	element->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	DatasetInput* weights_input   =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
	Input* thickness_input =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);
	Input* vx_input        =element->GetInput(VxEnum);                                 _assert_(vx_input);
	Input* vy_input        =element->GetInput(VyEnum);                                 _assert_(vy_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		element->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,ThicknessAcrossGradientEnum);
		thickness_input->GetInputValue(&thickness,gauss);
		thickness_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel = sqrt(vx*vx+vy*vy);
		vx  = vx/(vel+1.e-9);
		vy  = vy/(vel+1.e-9);

		/*J = 1/2 ( -vy*dH/dx + vx*dH/dy )^2 */
		Jelem+=weight*1/2*(-vy*dp[0] + vx*dp[1])*(-vy*dp[0] + vx*dp[1])*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
