/*!\file MpiDenseMumpsSolve.cpp
 * \brief: solve dense matrix system with MUMPS
 */

/*Header files: */
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/Numerics/types.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/io/Comm/IssmComm.h"
#include "../../classes/Params/GenericParam.h"
#include "../../classes/Params/Parameters.h"
#include "../mpi/issmmpi.h"
#include "../adolc/adolcincludes.h"
#include "./mumpsincludes.h"

/*Mumps header files: */
#include <dmumps_c.h>

void MumpsInit(DMUMPS_STRUC_C &theMumpsStruc){ 
	theMumpsStruc.par          = 1;  
	theMumpsStruc.sym          = 0;
	theMumpsStruc.comm_fortran = MPI_Comm_c2f(IssmComm::GetComm());
	theMumpsStruc.job          = -1;
	dmumps_c(&theMumpsStruc);
}

// must be preceded by a call to MumpsInit
void MumpsSettings(DMUMPS_STRUC_C &theMumpsStruc) { 
	/*Control statements:{{{ */
	theMumpsStruc.icntl[1-1] = 6; //error verbose: default 6, 0 or negative -> suppressed
	theMumpsStruc.icntl[2-1] = 0; //std verbose: default 1, 0 or negative -> suppressed
	theMumpsStruc.icntl[3-1] = 0; //global information verbose: default 6, 0 or negative -> suppressed
	theMumpsStruc.icntl[4-1] = 0; //verbose everything: default is 4
	theMumpsStruc.icntl[5-1] = 0;
	theMumpsStruc.icntl[18-1] = 3;

	theMumpsStruc.icntl[20-1] = 0;
	theMumpsStruc.icntl[21-1] = 0;
	theMumpsStruc.icntl[30-1] = 0;
	/*}}}*/
}

// must be preceded by a call to MumpsInit
void MumpsAnalyze(DMUMPS_STRUC_C &theMumpsStruc) { 
	theMumpsStruc.job          = 1;
	dmumps_c(&theMumpsStruc);
}

// must be preceded by a call to MumpsAnalyze
void MumpsFactorize(DMUMPS_STRUC_C &theMumpsStruc) { 
	theMumpsStruc.job          = 2;
	dmumps_c(&theMumpsStruc);
}

// must be preceded by a call to MumpsFactorize
void MumpsBacksubstitute(DMUMPS_STRUC_C &theMumpsStruc) { 
	theMumpsStruc.job          = 3;
	dmumps_c(&theMumpsStruc);
}

// must be preceded at least  by a call to MumpsInit
void MumpsFinalize(DMUMPS_STRUC_C &theMumpsStruc) { 
	theMumpsStruc.job          = -2;
	dmumps_c(&theMumpsStruc);
}

void MumpsSolve(int n,int nnz,int local_nnz,int* irn_loc,int* jcn_loc, IssmPDouble *a_loc, IssmPDouble *rhs, Parameters* parameters=0 /*unused here*/){
	/*Initialize mumps*/
	DMUMPS_STRUC_C theMumpsStruc;
	MumpsInit(theMumpsStruc);
	MumpsSettings(theMumpsStruc);

	/*now setup the rest of theMumpsStruc */
	theMumpsStruc.n=n;
	theMumpsStruc.nz=nnz;
	theMumpsStruc.nz_loc=local_nnz;
	theMumpsStruc.irn_loc=irn_loc;
	theMumpsStruc.jcn_loc=jcn_loc;
	theMumpsStruc.a_loc=a_loc;
	theMumpsStruc.rhs=rhs;
	theMumpsStruc.nrhs=1;
	theMumpsStruc.lrhs=1;

	/*Solve system*/
	MumpsAnalyze(theMumpsStruc);
	MumpsFactorize(theMumpsStruc);
	MumpsBacksubstitute(theMumpsStruc);
	MumpsFinalize(theMumpsStruc);
}

#ifdef _HAVE_ADOLC_
// prototype for active variant
void MumpsSolve(int n,
		int nnz,
		int local_nnz,
		int* irn_loc,
		int* jcn_loc,
		IssmDouble *a_loc,
		IssmDouble *rhs,
		Parameters* parameters);
#endif 

void MpiDenseMumpsSolve( /*output: */ IssmDouble* uf, int uf_M, int uf_m, /*matrix input: */ IssmDouble* Kff, int Kff_M, int Kff_N, int Kff_m, /*right hand side vector: */ IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters){ /*{{{*/

	/*Variables*/
	ISSM_MPI_Comm  comm;
	int            num_procs;
	int            i,j;
	int            nnz,local_nnz;
	int           *irn_loc    = NULL;
	int           *jcn_loc    = NULL;
	IssmDouble    *a_loc      = NULL;
	int            count;
	int            lower_row;
	int            upper_row;
	IssmDouble    *rhs        = NULL;
	int           *recvcounts = NULL;
	int           *displs     = NULL;

	/*Communicator info */
	num_procs = IssmComm::GetSize();
	comm      = IssmComm::GetComm();

	/*First, some checks*/
	if (Kff_M!=Kff_N)_error_("stiffness matrix Kff should be square");
	if (uf_M!=Kff_M || uf_M!=pf_M)_error_("solution vector should be the same size as stiffness matrix Kff and load vector pf");
	if (uf_m!=Kff_m || uf_m!=pf_m)_error_("solution vector should be locally the same size as stiffness matrix Kff and load vector pf");

	/*Initialize matrix */
	/*figure out number of non-zero entries: */
	local_nnz=0;
	for(i=0;i<Kff_m;i++){
		for(j=0;j<Kff_N;j++){
			if (Kff[i*Kff_N+j]!=0)local_nnz++;
		}
	}

	ISSM_MPI_Reduce(&local_nnz,&nnz,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,comm);
	ISSM_MPI_Bcast(&nnz,1,ISSM_MPI_INT,0,comm);

	/*Allocate: */
	if(local_nnz){
		irn_loc=xNew<int>(local_nnz);
		jcn_loc=xNew<int>(local_nnz);
		a_loc=xNew<IssmDouble>(local_nnz);
	}

	/*Populate the triplets: */
	GetOwnershipBoundariesFromRange(&lower_row,&upper_row,Kff_m,comm);
	count=0;
	for(i=0;i<Kff_m;i++){
		for(j=0;j<Kff_N;j++){
			if (Kff[i*Kff_N+j]!=0){
				irn_loc[count]=lower_row+i+1; //fortran indexing
				jcn_loc[count]=j+1; //fortran indexing
				a_loc[count]=Kff[i*Kff_N+j];
				count++;
			}
		}
	}
	/*Deal with right hand side. We need to ISSM_MPI_Gather it onto cpu 0: */
	rhs=xNew<IssmDouble>(pf_M);

	recvcounts=xNew<int>(num_procs);
	displs=xNew<int>(num_procs);

	/*recvcounts:*/
	ISSM_MPI_Allgather(&pf_m,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,comm);

	/*displs: */
	ISSM_MPI_Allgather(&lower_row,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,comm);

	/*Gather:*/
	ISSM_MPI_Gatherv(pf, pf_m, ISSM_MPI_DOUBLE, rhs, recvcounts, displs, ISSM_MPI_DOUBLE,0,comm);

	MumpsSolve(Kff_M,nnz,local_nnz,irn_loc,jcn_loc,a_loc,rhs,parameters);

	/*Now scatter from cpu 0 to all other cpus*/
	ISSM_MPI_Scatterv( rhs, recvcounts, displs, ISSM_MPI_DOUBLE, uf, uf_m, ISSM_MPI_DOUBLE, 0, comm); 

	/*Cleanup*/
	xDelete<int>(irn_loc);
	xDelete<int>(jcn_loc);
	xDelete<IssmDouble>(a_loc);
	xDelete<IssmDouble>(rhs);
	xDelete<int>(recvcounts);
	xDelete<int>(displs);
} /*}}}*/

#ifdef _HAVE_ADOLC_

int mumpsSolveEDF(int iArrLength, int* iArr, int /* ignored */, IssmPDouble* dp_x, int /* ignored */, IssmPDouble* dp_y) {
  // unpack parameters
  int n=iArr[0];
  int nnz=iArr[1];
  int local_nnz=iArr[2];
  int *local_irn=xNew<int>(local_nnz);
  int *local_jcn=xNew<int>(local_nnz);
  IssmPDouble *A=xNew<IssmPDouble>(local_nnz);
  for (int i=0;i<local_nnz;++i) {
    local_irn[i]=iArr[3+i];
    local_jcn[i]=iArr[3+local_nnz+i];
    A[i]=dp_x[i];
  }
  IssmPDouble *rhs_sol=xNew<IssmPDouble>(n);
  for (int i=0;i<n;++i) { 
    rhs_sol[i]=dp_x[local_nnz+i];
  }
  MumpsSolve(n,nnz,local_nnz,local_irn,local_jcn,A,rhs_sol);
  for (int i=0;i<n;++i) {
    dp_y[i]=rhs_sol[i];
  }
  xDelete(rhs_sol);
  xDelete(A);
  xDelete(local_jcn);
  xDelete(local_irn);
  return 0;
}

void MumpsSolve(int n,int nnz,int local_nnz,int* irn_loc,int* jcn_loc,IssmDouble *a_loc,IssmDouble *rhs,Parameters* parameters){
  int packedDimsSparseArrLength=1+1+1+local_nnz+local_nnz;
  int *packedDimsSparseArr=xNew<int>(packedDimsSparseArrLength);
  packedDimsSparseArr[0]=n;
  packedDimsSparseArr[1]=nnz;
  packedDimsSparseArr[2]=local_nnz;
  for (int i=0;i<local_nnz;++i) {
    packedDimsSparseArr[3+i]=irn_loc[i];
    packedDimsSparseArr[3+local_nnz+i]=jcn_loc[i];
  }
  IssmDouble *pack_A_rhs=xNew<IssmDouble>(local_nnz+n);
  for (int i=0;i<local_nnz;++i) { 
    pack_A_rhs[i]=a_loc[i];
  }
  for (int i=0;i<n;++i) { 
    pack_A_rhs[local_nnz+i]=rhs[i];
  }
  IssmPDouble *passivePack_A_rhs=xNew<IssmPDouble>(local_nnz+n);
  IssmPDouble *passiveSol=xNew<IssmPDouble>(n);
  IssmDouble *sol=xNew<IssmDouble>(n);
  call_ext_fct(xDynamicCast<GenericParam<Adolc_edf> * >(parameters->FindParamObject(AdolcParamEnum))->GetParameterValue().myEDF_for_solverx_p,
	       packedDimsSparseArrLength, packedDimsSparseArr,
	       local_nnz+n, passivePack_A_rhs, pack_A_rhs, 
	       n, passiveSol,sol);
  for (int i=0;i<n;++i) { 
    rhs[i]=sol[i];
  }
  xDelete(sol);
  xDelete(passiveSol);
  xDelete(passivePack_A_rhs);
  xDelete(pack_A_rhs);
  xDelete(packedDimsSparseArr);
}

int fos_reverse_mumpsSolveEDF(int iArrLength, int* iArr, 
			      int m, IssmPDouble *dp_U,
			      int nPlusNz, IssmPDouble *dp_Z,
			      IssmPDouble *dp_x, IssmPDouble *dp_y) {
  // unpack parameters
  int n=iArr[0];
  int nnz=iArr[1];
  int local_nnz=iArr[2];
  int *local_irn=xNew<int>(local_nnz);
  int *local_jcn=xNew<int>(local_nnz);
  IssmPDouble *a_loc=xNew<IssmPDouble>(local_nnz);
  for (int i=0;i<local_nnz;++i) {
	  local_irn[i]=iArr[3+i];
	  local_jcn[i]=iArr[3+local_nnz+i];
	  a_loc[i]=dp_x[i];
  }
  IssmPDouble *rhs_sol=xNew<IssmPDouble>(n);
  for (int i=0;i<n;++i) { 
    rhs_sol[i]=dp_U[i];
  }
  DMUMPS_STRUC_C theMumpsStruc;
  MumpsInit(theMumpsStruc);
  MumpsSettings(theMumpsStruc);
  theMumpsStruc.n=n;
  theMumpsStruc.nz=nnz;
  theMumpsStruc.nz_loc=local_nnz;
  theMumpsStruc.irn_loc=local_irn;
  theMumpsStruc.jcn_loc=local_jcn;
  theMumpsStruc.a_loc=a_loc;
  theMumpsStruc.rhs=rhs_sol;
  theMumpsStruc.nrhs=1;
  theMumpsStruc.lrhs=1;
  theMumpsStruc.icntl[9-1] = 0; //solve for the transpose
  MumpsAnalyze(theMumpsStruc);
  MumpsFactorize(theMumpsStruc);
  MumpsBacksubstitute(theMumpsStruc);
  MumpsFinalize(theMumpsStruc);
  // update the adjoint of the rhs:
  for (int i=0;i<n;++i) {
    dp_Z[local_nnz+i]+=rhs_sol[i];
  }
  // bcast dp_y (the solution of the forward system)
  ISSM_MPI_Bcast(dp_y,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
  // bcast the adjoint of the right-hand-side, i.e. this solution
  ISSM_MPI_Bcast(rhs_sol,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
  // update the adjoint of the matrix with the outer product of right-hand-side adjoint and the original solution
  for (int i=0;i<local_nnz;++i) {
    dp_Z[i]-=rhs_sol[iArr[3+i]-1]*dp_y[iArr[3+local_nnz+i]-1];
  }
  xDelete(rhs_sol);
  xDelete(a_loc);
  xDelete(local_jcn);
  xDelete(local_irn);
  return 3;
}

#endif
