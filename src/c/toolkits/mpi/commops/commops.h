/*! \file mpipatches.h
 *  \brief: prototype header for all ISSM add-ons to MPI
 */

#ifndef MPI_PATCHES_H_ 
#define MPI_PATCHES_H_

#include "../../../shared/Numerics/types.h" 
#include "../../../shared/io/Comm/IssmComm.h"
#include "../../mpi/issmmpi.h"

int DetermineLocalSize(int global_size,ISSM_MPI_Comm comm);
int* DetermineRowRankFromLocalSize(int global_size,int localsize,ISSM_MPI_Comm comm);
void GetOwnershipBoundariesFromRange(int* plower_row,int* pupper_row,int range,ISSM_MPI_Comm comm);
int DetermineGlobalSize(int local_size,ISSM_MPI_Comm comm);

#endif
