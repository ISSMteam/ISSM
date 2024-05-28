/*!\file SurfaceAreax
 * \brief: compute Surface area
 */

#include "./SurfaceAreax.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"

void SurfaceAreax(IssmDouble* pS,FemModel* femmodel){

	/*Intermediary*/
	Element* element=NULL;

	/*output: */
	IssmDouble S = 0.;
	IssmDouble S_sum;

	/*Compute gradients: */
	for(Object* & object : femmodel->elements->objects){
		element=xDynamicCast<Element*>(object);
		S+=element->SurfaceArea();
	}

	/*Sum all J from all cpus of the cluster:*/
 	ISSM_MPI_Reduce (&S,&S_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&S_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm()); 
	S=S_sum;

	/*add surface area to element inputs:*/
	InputUpdateFromConstantx(femmodel,S,SurfaceAreaEnum);

	/*Assign output pointers: */
	if(pS) *pS=S;
}
