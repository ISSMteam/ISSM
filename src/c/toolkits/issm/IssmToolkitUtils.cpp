/*!\file:  IssmToolkitUtils.cpp
 * \brief utilities used throughout our ISSM toolkit
 */ 

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/MemOps/MemOps.h"
#include "../../shared/io/Comm/IssmComm.h"
#include "../../shared/Enum/Enum.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../ToolkitOptions.h"
#include "./IssmToolkitUtils.h"
#include <string.h>
/*}}}*/

/*Routines: */
int IssmMatTypeFromToolkitOptions(void){ /*{{{*/

	char *mat_type      = NULL;
	int   mat_type_enum;
	int   num_procs     = 0;
	bool  isparallel    = false;

	/*first, figure out if we are running in parallel: */
	num_procs=IssmComm::GetSize();
	if(num_procs>1)isparallel=true;

	/*retrieve matrix type as a string, from the Toolkits Options database, similar to what Petsc does. Actually, 
	 *we try and stick with the Petsc matrix types: */
	mat_type=ToolkitOptions::GetToolkitOptionValue("mat_type");

	if (strcmp(mat_type,"mpidense")==0){
		mat_type_enum=MpiDenseEnum;
	}
	else if (strcmp(mat_type,"mpisparse")==0){
		mat_type_enum=MpiSparseEnum;
	}
	else if (strcmp(mat_type,"dense")==0){
		if (isparallel) _error_("Dense matrix type not supported for parallel runs with num_procs>1");
		else mat_type_enum=DenseEnum;
	}
	else _error_("matrix type not supported yet!");

	/*Free resources: */
	xDelete<char>(mat_type);

	/*return: */
	return mat_type_enum;
} /*}}}*/
int IssmVecTypeFromToolkitOptions(void){ /*{{{*/

	char* vec_type=NULL;
	int   vec_type_enum;
	int   num_procs=0;
	bool  isparallel=false;

	/*first, figure out if we are running in parallel: */
	num_procs=IssmComm::GetSize();
	if(num_procs>1)isparallel=true;

	/*retrieve vector type as a string, from the Toolkits Options database, similar to what Petsc does. Actually, 
	 *we try and stick with the Petsc vector types: */
	vec_type=ToolkitOptions::GetToolkitOptionValue("vec_type");

	if (strcmp(vec_type,"mpi")==0){
		vec_type_enum=MpiEnum;
	}
	else if (strcmp(vec_type,"seq")==0){
		if (isparallel) _error_("Dense vector type not supported for parallel runs with num_procs>1");
		else vec_type_enum=SeqEnum;
	}
	else _error_("vector type not supported yet!");

	/*Free resources: */
	xDelete<char>(vec_type);

	/*return: */
	return vec_type_enum;
} /*}}}*/  
int IssmSolverTypeFromToolkitOptions(void){ /*{{{*/

	int   solver_type_enum;
	bool  isparallel=false;

	/*first, figure out if we are running in parallel: */
	int num_procs=IssmComm::GetSize();
	if(num_procs>1)isparallel=true;

	/*retrieve solver type as a string, from the Toolkits Options database, similar to what Petsc does. Actually, 
	 *we try and stick with the Petsc vector types: */
	char* solver_type=ToolkitOptions::GetToolkitOptionValue("solver_type");
	if(!solver_type) _error_("Solver not set");

	if (strcmp(solver_type,"mumps")==0){
		solver_type_enum=MumpsEnum;
	}
	else if (strcmp(solver_type,"gsl")==0){
		if (isparallel) _error_("Gsl solver type not supported for parallel runs with num_procs>1");
		else solver_type_enum=GslEnum;
	}
	else _error_("solver type not supported yet!");

	/*Free resources: */
	xDelete<char>(solver_type);

	/*return: */
	return solver_type_enum;
} /*}}}*/  
