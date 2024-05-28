/*!\file global.h:
 * \brief: these are the global variables always needed. 
 */

#ifndef _GLOBALS_H_
#define _GLOBALS_H_
#include "../shared/io/Comm/IssmComm.h"
#include "../toolkits/ToolkitOptions.h"

/*Communicators: */
#ifndef _DO_NOT_LOAD_GLOBALS_ 
ISSM_MPI_Comm IssmComm::comm;
bool IssmComm::parallel;

/*String that is used to characterize our toolkits, ends up in Petsc Options
 * database if we use Petsc. Can also be used to characterize the ISSM toolkit,
 * often used when Petsc is not allowed*/
char* ToolkitOptions::toolkittype;
char* ToolkitOptions::toolkitoptions;
#endif

#endif
