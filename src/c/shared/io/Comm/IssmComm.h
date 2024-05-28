/* \file IssmComm.h
 * \brief  create a class with a static comm, and static methods to access it
 * This is a way of protecting access to the communicator.
 */

#ifndef _ISSM_COMM_H
#define _ISSM_COMM_H

/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../../toolkits/mpi/issmmpi.h"

/*}}}*/

class IssmComm {

	private:
		static ISSM_MPI_Comm comm;
		static bool parallel;

	public:
		static void SetComm(ISSM_MPI_Comm incomm);
		static void SetComm(void);
		static ISSM_MPI_Comm GetComm(void);
		static int GetRank(void);
		static int GetSize(void);
};

#endif  /* _ISSM_COMM_H */
