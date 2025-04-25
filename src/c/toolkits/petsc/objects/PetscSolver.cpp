/*!\file SolverxPetsc
 * \brief Petsc implementation of solver
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./PetscSolver.h"
#include "../../../shared/Numerics/Verbosity.h"
#include "../../../shared/MemOps/MemOps.h"
#include "../../../shared/Exceptions/exceptions.h"
#include "../../../shared/io/Comm/IssmComm.h"
#include "../../../shared/Enum/Enum.h"
#include "../../../shared/io/Print/Print.h"
#include "../../../shared/Numerics/recast.h"

void	PetscSolve(PetscVec<IssmDouble>** puf, PetscMat<IssmDouble>* Kff, PetscVec<IssmDouble>* pf, PetscVec<IssmDouble>* uf0,PetscVec<IssmDouble>* df, Parameters* parameters) { /*{{{*/

	PetscVec<IssmDouble>* uf=new PetscVec<IssmDouble>();

	PVec uf0_vector = NULL;
	PVec df_vector  = NULL;

	if(uf0) uf0_vector = uf0->vector;
	if(df)  df_vector  = df->vector;

	SolverxPetsc(&uf->vector, Kff->matrix, pf->vector, uf0_vector, df_vector, parameters);

	*puf=uf;

}
/*}}}*/
void	SolverxPetsc(PVec* puf, PMat Kff, PVec pf, PVec uf0,PVec df, Parameters* parameters){ /*{{{*/

	/*Output: */
	PVec        uf = NULL;

	/*Intermediary: */
	int        local_m,local_n,global_m,global_n;

	/*Solver */
	PKSP        ksp              = NULL;
	PC         pc               = NULL;
	int        iteration_number;
	int        solver_type;
	bool       fromlocalsize    = true;
	#if PETSC_VERSION_LT(3,2,0)
	PetscTruth flag,flg;
	#else
	PetscBool flag,flg;
	#endif

	/*FS: */
	IS         isv=NULL;
	IS         isp=NULL;
	char ksp_type[50];

	/*Display message*/
	#if PETSC_VERSION_LT(3,2,0)
	if(VerboseSolver())PetscOptionsPrint(stdout);
	#else
		#if PETSC_VERSION_LT(3,7,0)
		if(VerboseSolver())PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
		#else
		if(VerboseSolver())PetscOptionsView(NULL,PETSC_VIEWER_STDOUT_WORLD);
		#endif
	#endif

	/*First, check that f-set is not NULL, i.e. model is fully constrained:*/ 
	_assert_(Kff);
	MatGetSize(Kff,&global_m,&global_n); _assert_(global_m==global_n);
	if(!global_n){
		*puf=NewVec<PVec>(0,IssmComm::GetComm()); return;
	}

	/*Initial guess */
	/*Now, check that we are not giving an initial guess to the solver, if we are running a direct solver: */
	#if PETSC_VERSION_LT(3,7,0)
	PetscOptionsGetString(PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
	#elif PETSC_VERSION_LT(3,19,0)
	PetscOptionsGetString(NULL,PETSC_NULL,"-ksp_type",ksp_type,49,&flg);
	#else
	PetscOptionsGetString(NULL,PETSC_NULLPTR,"-ksp_type",ksp_type,49,&flg);
	#endif
	if(flg!=PETSC_TRUE) _error_("could not find option -ksp_type, maybe you are not using the right toolkit?");
	if (strcmp(ksp_type,"preonly")==0) uf0=NULL;

	/*If initial guess for the solution exists, use it to create uf, otherwise, 
	 * duplicate the right hand side so that the solution vector has the same structure*/
	if(uf0){
		VecDuplicate(uf0,&uf); VecCopy(uf0,uf);
	}
	else{
		MatGetLocalSize(Kff,&local_m,&local_n);uf=NewVec<PVec>(local_n,IssmComm::GetComm(),fromlocalsize);
	}

	/*Process petsc options to see if we are using special types of external solvers*/
	PetscOptionsDetermineSolverType(&solver_type);

	/*Check the solver is available*/
	if(solver_type==MUMPSPACKAGE_LU || solver_type==MUMPSPACKAGE_CHOL){
		#ifndef _HAVE_MUMPS_
		_error_("requested MUMPS solver, which was not compiled into ISSM!\n");
		#endif
	}

	/*Prepare solver*/
	KSPCreate(IssmComm::GetComm(),&ksp);
	#if PETSC_VERSION_GE(3,5,0)
		KSPSetOperators(ksp,Kff,Kff);
	#else
		KSPSetOperators(ksp,Kff,Kff,DIFFERENT_NONZERO_PATTERN);
	#endif
	KSPSetFromOptions(ksp);

	/*Specific solver?: */
	KSPGetPC(ksp,&pc);
	if (solver_type==MUMPSPACKAGE_LU){
		#if PETSC_VERSION_GE(3,9,0)
		PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);
		#else
		PCFactorSetMatSolverPackage(pc,MATSOLVERMUMPS);
		#endif
	}

	/*FS: */
	if (solver_type==FSSolverEnum){
		/*Make indices out of doftypes: */
		if(!df)_error_("need doftypes for FS solver!\n");
		DofTypesToIndexSet(&isv,&isp,df,FSSolverEnum);

		/*Set field splits: */
		KSPGetPC(ksp,&pc);

		#if PETSC_VERSION_LT(3,1,0)
		PCFieldSplitSetIS(pc,isv);
		PCFieldSplitSetIS(pc,isp);
		#elif PETSC_VERSION_LT(3,19,0)
		PCFieldSplitSetIS(pc,PETSC_NULL,isv);
		PCFieldSplitSetIS(pc,PETSC_NULL,isp);
		#else
		PCFieldSplitSetIS(pc,PETSC_NULLPTR,isv);
		PCFieldSplitSetIS(pc,PETSC_NULLPTR,isp);
		#endif

	}

	/*If there is an initial guess for the solution, use it
	 * except if we are using the MUMPS direct solver
	 * where any initial guess will crash Petsc*/
	if (uf0){
		if((solver_type!=MUMPSPACKAGE_LU) && (solver_type!=MUMPSPACKAGE_CHOL) && (solver_type!=SPOOLESPACKAGE_LU) && (solver_type!=SPOOLESPACKAGE_CHOL) && (solver_type!=SUPERLUDISTPACKAGE)){
			KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
		}
	}

	/*Solve: */
	if(VerboseSolver())KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	KSPSolve(ksp,pf,uf);

	/*Check convergence*/
	KSPGetIterationNumber(ksp,&iteration_number);
	if (iteration_number<0) _error_("Solver diverged at iteration number: " << -iteration_number);
	if (VerboseSolver())  _printf0_("Petsc: "<< iteration_number << " KSP iterations\n"); 

	/*Free resources:*/
	KSPFree(&ksp);

	/*Assign output pointers:*/
	*puf=uf;
} 
/*}}}*/
void DofTypesToIndexSet(IS* pisv, IS* pisp, PVec df,int typeenum){ /*{{{*/

	/*output: */
	IS          isv=NULL;
	IS          isp=NULL;

	int         start,end;
	PArray			df_local;
	int         df_local_size;

	int*     pressure_indices=NULL;
	int*     velocity_indices=NULL;
	int      pressure_num=0;
	int      velocity_num=0;
	int      pressure_count=0;
	int      velocity_count=0;

	if(typeenum==FSSolverEnum){

		/*Ok, recover doftypes vector values and indices: */
		VecGetOwnershipRange(df,&start,&end);
		VecGetLocalSize(df,&df_local_size);
		VecGetArray(df,&df_local);

		pressure_num=0;
		velocity_num=0;
		for(int i=0;i<df_local_size;i++){
			if (reCast<int, IssmDouble>(df_local[i])==PressureEnum)pressure_num++;
			else velocity_num++;
		}

		/*Allocate indices: */
		if(pressure_num)pressure_indices=xNew<int>(pressure_num);
		if(velocity_num)velocity_indices=xNew<int>(velocity_num);

		pressure_count=0;
		velocity_count=0;
		for(int i=0;i<df_local_size;i++){
			if(reCast<int, IssmDouble>(df_local[i])==PressureEnum){
				pressure_indices[pressure_count]=start+i;
				pressure_count++;
			}
			if(reCast<int, IssmDouble>(df_local[i])==VelocityEnum){
				velocity_indices[velocity_count]=start+i;
				velocity_count++;
			}
		}
		VecRestoreArray(df,&df_local);

		/*Create indices sets: */
		#if PETSC_VERSION_LT(3,2,0)
		ISCreateGeneral(IssmComm::GetComm(),pressure_num,pressure_indices,&isp);
		ISCreateGeneral(IssmComm::GetComm(),velocity_num,velocity_indices,&isv);
		#else
		ISCreateGeneral(IssmComm::GetComm(),pressure_num,pressure_indices,PETSC_COPY_VALUES,&isp);
		ISCreateGeneral(IssmComm::GetComm(),velocity_num,velocity_indices,PETSC_COPY_VALUES,&isv);
		#endif
	}

	/*Free resources:*/
	xDelete<int>(pressure_indices);
	xDelete<int>(velocity_indices);

	/*Assign output pointers:*/
	*pisv=isv;
	*pisp=isp;
}
/*}}}*/
