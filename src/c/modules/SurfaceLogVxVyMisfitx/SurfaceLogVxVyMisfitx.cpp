/*!\file SurfaceLogVxVyMisfitx
 * \brief: compute misfit between observations and model
 */

#include <float.h> /*defines DBL_EPSILON*/
#include "./SurfaceLogVxVyMisfitx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/DatasetInput.h"

void SurfaceLogVxVyMisfitx( IssmDouble* pJ, Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials,Parameters* parameters){

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	/*Compute Misfit: */
	for(Object* & object : elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		J+=SurfaceLogVxVyMisfit(element);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}

IssmDouble SurfaceLogVxVyMisfit(Element* element){

	int        domaintype,numcomponents;
	IssmDouble Jelem=0.;
	IssmDouble epsvel=DBL_EPSILON;
	IssmDouble meanvel=3.170979198376458e-05; /*1000 m/yr*/
	IssmDouble misfit,Jdet;
	IssmDouble vx,vy,vxobs,vyobs,weight;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnSurface()) return 0.;

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

	/*Spawn surface element*/
	Element* topelement = element->SpawnTopElement();

	/* Get node coordinates*/
	topelement->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	DatasetInput* weights_input=topelement->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
   Input* vx_input     = NULL;
   Input* vxobs_input  = topelement->GetInput(InversionVxObsEnum);   _assert_(vxobs_input);
   Input* vy_input     = NULL;
   Input* vyobs_input  = NULL;

   /*Read SurfaceEnum from 2D models:SSA, L1L2, MOLHO*/
   if (domaintype == Domain2DhorizontalEnum) {
      vx_input = topelement->GetInput(VxSurfaceEnum);             _assert_(vx_input);
      if(numcomponents==2){
         vy_input    =topelement->GetInput(VySurfaceEnum);        _assert_(vy_input);
         vyobs_input =topelement->GetInput(InversionVyObsEnum);   _assert_(vyobs_input);
      }
	}
   else {
      vx_input = topelement->GetInput(VxEnum);                    _assert_(vx_input);
      if(numcomponents==2){
         vy_input    =topelement->GetInput(VyEnum);               _assert_(vy_input);
         vyobs_input =topelement->GetInput(InversionVyObsEnum);   _assert_(vyobs_input);
      }
   }

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(4);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		weights_input->GetInputValue(&weight,gauss,SurfaceLogVxVyMisfitEnum);
		vx_input->GetInputValue(&vx,gauss);
		vxobs_input->GetInputValue(&vxobs,gauss);
		if(numcomponents==2){
			vy_input->GetInputValue(&vy,gauss);
			vyobs_input->GetInputValue(&vyobs,gauss);
		}

		/*Compute SurfaceRelVelMisfit:
		 *
		 *      1            [        |u| + eps     2          |v| + eps     2  ]
		 * J = --- \bar{v}^2 | log ( -----------  )   +  log ( -----------  )   |
		 *      2            [       |u    |+ eps              |v    |+ eps     ]
		 *                              obs                       obs
		 */

		if(numcomponents==1){
			misfit=0.5*meanvel*meanvel*pow(log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)),2);
		}
		else{
			misfit=0.5*meanvel*meanvel*(
						pow(log((fabs(vx)+epsvel)/(fabs(vxobs)+epsvel)),2) +
						pow(log((fabs(vy)+epsvel)/(fabs(vyobs)+epsvel)),2) );
		}

		/*Add to cost function*/
		Jelem+=misfit*weight*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}
