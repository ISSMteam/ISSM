/*! \file IssmComm.cpp
 * \brief  file containing the methods for IssmComm.h
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./IssmComm.h"
#include "../../Numerics/types.h"
#include "../../Exceptions/exceptions.h"

#ifdef _DO_NOT_LOAD_GLOBALS_
ISSM_MPI_Comm	IssmComm::comm;
bool			IssmComm::parallel;
#endif

void IssmComm::SetComm(ISSM_MPI_Comm incomm){ /*{{{*/

	/*A comm is provided, we are running in parallel (this is not a module)*/
	parallel = true;
	comm     = incomm;

}/*}}}*/
void IssmComm::SetComm(void){ /*{{{*/

	/*no comm provided, This is a matlab/python module*/
	parallel = false;

	/*No need to initialise comm*/

}/*}}}*/
ISSM_MPI_Comm IssmComm::GetComm(){  /*{{{*/
	if(!parallel) _error_("Cannot return comm in serial mode");
	return comm;
}/*}}}*/
int IssmComm::GetRank(){  /*{{{*/

	int my_rank = 0;

	/*for matlab and python modules*/
	if(!parallel) return my_rank;

	ISSM_MPI_Comm_rank(comm,&my_rank);

	return my_rank;

}/*}}}*/
int IssmComm::GetSize(){  /*{{{*/

	int size = 1;

	/*for matlab and python modules*/
	if(!parallel) return size;

	ISSM_MPI_Comm_size(comm,&size);

	return size;

}/*}}}*/
