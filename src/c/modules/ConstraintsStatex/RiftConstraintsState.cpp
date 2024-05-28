/*!\file RiftConstraintsState.cpp
 * \brief: manage penalties for rifts
 */
#include "./ConstraintsStateLocal.h"
#include "../../shared/shared.h"

/*current module: */
/*RiftConstraintsState(int* pconverged, int* pnum_unstable_constraints,Loads* loads,int min_mechanical_constraints,int configuration_type){{{*/
void RiftConstraintsState(int* pconverged, int* pnum_unstable_constraints,Loads* loads,int min_mechanical_constraints,int configuration_type){

	int num_unstable_constraints=0;
	int converged=0;

	RiftConstrain(&num_unstable_constraints,loads,configuration_type);
	if(num_unstable_constraints==0)converged=1;

	if(RiftIsFrozen(loads,configuration_type)){
		converged=1;
		num_unstable_constraints=0;
	}
	else if(num_unstable_constraints<=min_mechanical_constraints){
		if(VerboseModule()) _printf0_("   freezing constraints\n");
		RiftFreezeConstraints(loads,configuration_type);
	}

	/*Assign output pointers: */
	*pconverged=converged;
	*pnum_unstable_constraints=num_unstable_constraints;
}
/*}}}*/
/*RiftConstrain(int* pnum_unstable_constraints,Loads* loads,int configuration_type){{{*/
void RiftConstrain(int* pnum_unstable_constraints,Loads* loads,int configuration_type){

	/* generic object pointer: */
	Riftfront* riftfront=NULL;
	Load*      load=NULL;

	int unstable;
	int sum_num_unstable_constraints;
	int num_unstable_constraints=0;

	/*Enforce constraints: */
	for (Object* & object : loads->objects){
		if (RiftfrontEnum==object->ObjectEnum()){
			load=(Load*)object;
			riftfront=(Riftfront*)load;
			riftfront->Constrain(&unstable);
			num_unstable_constraints+=unstable;
		}
	}

	ISSM_MPI_Reduce (&num_unstable_constraints,&sum_num_unstable_constraints,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_num_unstable_constraints,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	num_unstable_constraints=sum_num_unstable_constraints;

	/*Assign output pointers: */
	*pnum_unstable_constraints=num_unstable_constraints;

}
/*}}}*/
/*RiftIsFrozen(Loads* loads,int configuration_type){{{*/
int RiftIsFrozen(Loads* loads,int configuration_type){

	int			i;

	/* generic object pointer: */
	Load*      load=NULL;
	Riftfront* riftfront=NULL;
	int found=0;
	int mpi_found=0;

	/*Enforce constraints: */
	for (Object* & object : loads->objects){
		if (RiftfrontEnum==object->ObjectEnum()){
			load=(Load*)object;
			riftfront=(Riftfront*)load;
			if (riftfront->IsFrozen()){
				found=1;
				break;
			}
		}
	}

	/*Is there just one found? that would mean we have frozen! : */
	ISSM_MPI_Reduce (&found,&mpi_found,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&mpi_found,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	found=mpi_found;

	return found;
}
/*}}}*/
/*RiftFreezeConstraints(Loads* loads,int configuration_type){{{*/
void RiftFreezeConstraints(Loads* loads,int configuration_type){

	int			i;

	/* generic object pointer: */
	Load*      load=NULL;
	Riftfront* riftfront=NULL;

	/*Enforce constraints: */
	for (Object* & object : loads->objects){
		if (RiftfrontEnum==object->ObjectEnum()){
			load=(Load*)object;
			riftfront=(Riftfront*)load;
			riftfront->FreezeConstraints();
		}
	}

}
/*}}}*/
