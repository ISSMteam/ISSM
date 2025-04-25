/*!\file MumpsSolve.cpp
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
#include "../codipack/codipackincludes.h"
#include "../issm/SparseRow.h"
#include "./mumpsincludes.h"

/*Mumps header files: */
#include <dmumps_c.h>

void MumpsInit(DMUMPS_STRUC_C &theMumpsStruc){ 
	theMumpsStruc.n = 0;
	theMumpsStruc.nz = 0;
	theMumpsStruc.a = NULL;
	theMumpsStruc.jcn = NULL;
	theMumpsStruc.irn = NULL;
	theMumpsStruc.par          = 1;  
	theMumpsStruc.sym          = 0;
	#ifdef _HAVE_MPI_
	theMumpsStruc.comm_fortran = MPI_Comm_c2f(IssmComm::GetComm());
	#endif
	theMumpsStruc.job          = -1;
	dmumps_c(&theMumpsStruc);
}

// must be preceded by a call to MumpsInit
void MumpsSettings(DMUMPS_STRUC_C &theMumpsStruc) { 
	/*Control statements:{{{ */
	theMumpsStruc.icntl[1-1] = 0; //error verbose: default 6, 0 or negative -> suppressed
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

#ifdef _HAVE_AD_
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

#ifdef _HAVE_MPI_
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
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
	rhs=xNew<IssmDouble>(pf_M,"t");
#else
	rhs=xNew<IssmDouble>(pf_M);
#endif

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
void MpiSparseMumpsSolve( /*output: */ IssmDouble* uf, int uf_M, int uf_m, /*matrix input: */ SparseRow<IssmDouble>** Kff, int Kff_M, int Kff_N, int Kff_m, /*right hand side vector: */ IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters){ /*{{{*/

	/*Variables*/
	ISSM_MPI_Comm  comm;
	int            num_procs;
	int            i;
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
	if(Kff_M!=Kff_N)_error_("stiffness matrix Kff should be square");
	if(uf_M!=Kff_M || uf_M!=pf_M)_error_("solution vector should be the same size as stiffness matrix Kff and load vector pf");
	if(uf_m!=Kff_m || uf_m!=pf_m)_error_("solution vector should be locally the same size as stiffness matrix Kff and load vector pf");

	/*Initialize matrix */
	/*figure out number of non-zero entries: */
	local_nnz=0;
	for(i=0;i<Kff_m;i++){
		local_nnz+=Kff[i]->Nnz();
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
		Kff[i]->SetIrnJcnA(irn_loc,jcn_loc,a_loc,lower_row+i+1,count);
		count+=Kff[i]->Nnz();
	}
	/*Deal with right hand side. We need to ISSM_MPI_Gather it onto cpu 0: */
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
	rhs=xNew<IssmDouble>(pf_M,"t");
#else
	rhs=xNew<IssmDouble>(pf_M);
#endif

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
#endif
void SeqDenseMumpsSolve(IssmDouble* uf,int uf_M,int uf_m, IssmDouble* Kff,int Kff_M, int Kff_N, int Kff_m, IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters){/*{{{*/

	/*Variables*/
	int        *irn_loc = NULL;
	int        *jcn_loc = NULL;
	IssmDouble *a_loc   = NULL;
	IssmDouble *rhs     = NULL;

	/*First, some checks*/
	if (Kff_M!=Kff_N)_error_("stiffness matrix Kff should be square");
	if (Kff_M!=Kff_m)_error_("stiffness matrix Kff is not serial");
	if (uf_M!=Kff_M || uf_M!=pf_M)_error_("solution vector should be the same size as stiffness matrix Kff and load vector pf");
	if (uf_m!=Kff_m || uf_m!=pf_m)_error_("solution vector should be locally the same size as stiffness matrix Kff and load vector pf");

	/*Initialize matrix */
	/*figure out number of non-zero entries: */
	int nnz = 0;
	for(int i=0;i<Kff_M;i++){
		for(int j=0;j<Kff_N;j++){
			if(Kff[i*Kff_N+j]!=0)nnz++;
		}
	}

	/*Allocate: */
	if(nnz){
		irn_loc = xNew<int>(nnz);
		jcn_loc = xNew<int>(nnz);
		a_loc   = xNew<IssmDouble>(nnz);
	}

	/*Populate the triplets: */
	int count=0;
	for(int i=0;i<Kff_M;i++){
		for(int j=0;j<Kff_N;j++){
			if(Kff[i*Kff_N+j]!=0){
				irn_loc[count] = i+1; //fortran indexing
				jcn_loc[count] = j+1; //fortran indexing
				a_loc[count]   = Kff[i*Kff_N+j];
				count++;
			}
		}
	}

	/*Deal with right hand side. We need to ISSM_MPI_Gather it onto cpu 0: */
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
	rhs=xNew<IssmDouble>(pf_M,"t");
#else
	rhs=xNew<IssmDouble>(pf_M);
#endif
	xMemCpy<IssmDouble>(rhs,pf,Kff_M);

	MumpsSolve(Kff_M,nnz,nnz,irn_loc,jcn_loc,a_loc,rhs,parameters);

	/*Now scatter from cpu 0 to all other cpus*/
	xMemCpy<IssmDouble>(uf,rhs,Kff_M);

	/*Cleanup*/
	xDelete<int>(irn_loc);
	xDelete<int>(jcn_loc);
	xDelete<IssmDouble>(a_loc);
	xDelete<IssmDouble>(rhs);
}/*}}}*/

#ifdef _HAVE_ADOLC_
int mumpsSolveEDF(int iArrLength, int* iArr, int /* ignored */, IssmPDouble* dp_x, int /* ignored */, IssmPDouble* dp_y){/*{{{*/
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
}/*}}}*/
void MumpsSolve(int n,int nnz,int local_nnz,int* irn_loc,int* jcn_loc,IssmDouble *a_loc,IssmDouble *rhs,Parameters* parameters){/*{{{*/
  int packedDimsSparseArrLength=1+1+1+local_nnz+local_nnz;
  int *packedDimsSparseArr=xNew<int>(packedDimsSparseArrLength);
  packedDimsSparseArr[0]=n;
  packedDimsSparseArr[1]=nnz;
  packedDimsSparseArr[2]=local_nnz;
  for (int i=0;i<local_nnz;++i) {
    packedDimsSparseArr[3+i]=irn_loc[i];
    packedDimsSparseArr[3+local_nnz+i]=jcn_loc[i];
  }
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
  IssmDouble *pack_A_rhs=xNew<IssmDouble>(local_nnz+n,"t");
#else
  IssmDouble *pack_A_rhs=xNew<IssmDouble>(local_nnz+n);
#endif

  for (int i=0;i<local_nnz;++i) { 
    pack_A_rhs[i]=a_loc[i];
  }
  for (int i=0;i<n;++i) { 
    pack_A_rhs[local_nnz+i]=rhs[i];
  }
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
  IssmDouble *sol=xNew<IssmDouble>(n,"t");
#else
  IssmDouble *sol=xNew<IssmDouble>(n);
#endif

  call_ext_fct(xDynamicCast<GenericParam<Adolc_edf> * >(parameters->FindParamObject(AdolcParamEnum))->GetParameterValue().myEDF_for_solverx_p,
	       packedDimsSparseArrLength, packedDimsSparseArr,
	       local_nnz+n, pack_A_rhs, 
	       n,sol);
  for (int i=0;i<n;++i) { 
    rhs[i]=sol[i];
  }
  xDelete(sol);
  xDelete(pack_A_rhs);
  xDelete(packedDimsSparseArr);
}/*}}}*/
int fos_reverse_mumpsSolveEDF(int iArrLength, int* iArr, int m, IssmPDouble *dp_U, int nPlusNz, IssmPDouble *dp_Z, IssmPDouble *dp_x, IssmPDouble *dp_y) {/*{{{*/
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
}/*}}}*/
#endif

#ifdef _HAVE_CODIPACK_
#if _CODIPACK_MAJOR_==2
using Tape = typename IssmDouble::Tape;
using AccessInterface = codi::VectorAccessInterface<typename Tape::Real, typename Tape::Identifier>;
void MumpsSolve_codi_b(Tape* tape, void* data_in, AccessInterface* ra) {/*{{{*/

	/*recast data_in and tape*/
  codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)data_in;

  IssmDouble::Real* valueATrans;
  IssmDouble::Identifier* indexATrans;
  int* irnATrans;
  int* jcnATrans;
  IssmDouble::Identifier* indexB;
  IssmDouble::Real* valueX;
  IssmDouble::Identifier* indexX;
  int n;
  int nnz;
  int local_nnz;
  Parameters* parameters;

  data->getData(valueATrans);
  data->getData(indexATrans);
  data->getData(irnATrans);
  data->getData(jcnATrans);
  data->getData(indexB);
  data->getData(valueX);
  data->getData(indexX);
  data->getData(n);
  data->getData(nnz);
  data->getData(local_nnz);
  data->getData(parameters);

  // create the adjoint vector for x and reset the adjoint values on the tape
  IssmDouble::Gradient* adjX = xNew<IssmDouble::Gradient>(n);
  getVectorAdjoint(*tape, indexX, adjX, n);

  MumpsSolve(n, nnz, local_nnz, irnATrans, jcnATrans, valueATrans, adjX, parameters);
  // adjX contains now the solution

  updateVectorAdjoint(*tape, indexB, adjX, n);

  // bcast dp_y (the solution of the forward system)
  ISSM_MPI_Bcast(valueX,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
  // bcast the adjoint of the right-hand-side, i.e. this solution
  ISSM_MPI_Bcast(adjX,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

  for(int i=0; i<local_nnz; ++i) {
    // we access the transposed matrix here because we stored the indices in a transposed way
    // -1 is substracted because jcn and irn are stored with fortran indexing
    if(indexATrans[i] != 0) {
      updateAdjoint(*tape, indexATrans[i], -adjX[jcnATrans[i]-1]*valueX[irnATrans[i]-1]);
    }
  }

  xDelete(adjX);
}
/*}}}*/
void MumpsSolve_codi_delete(Tape* tape,void* data_in) {/*{{{*/

	/*recast data_in*/
	codi::ExternalFunctionUserData* data = (codi::ExternalFunctionUserData*)data_in;

  IssmDouble::Real* valueATrans;
  IssmDouble::Identifier* indexATrans;
  int* irnATrans;
  int* jcnATrans;
  IssmDouble::Identifier* indexB;
  IssmDouble::Real* valueX;
  IssmDouble::Identifier* indexX;
  int n;
  int nnz;
  int local_nnz;
  Parameters* parameters;

  data->getData(valueATrans);
  data->getData(indexATrans);
  data->getData(irnATrans);
  data->getData(jcnATrans);
  data->getData(indexB);
  data->getData(valueX);
  data->getData(indexX);
  data->getData(n);
  data->getData(nnz);
  data->getData(local_nnz);
  data->getData(parameters);

  xDelete(valueATrans);
  xDelete(indexATrans);
  xDelete(irnATrans);
  xDelete(jcnATrans);
  xDelete(indexB);
  xDelete(valueX);
  xDelete(indexX);

  delete data;
}
/*}}}*/
void MumpsSolve(int n,int nnz,int local_nnz,int* irn_loc,int* jcn_loc,IssmDouble *a_loc,IssmDouble *rhs,Parameters* parameters){/*{{{*/
  IssmDouble::Tape& tape = IssmDouble::getTape();
  codi::ExternalFunctionUserData* dataHandler = NULL;

  if(tape.isActive()) {
    dataHandler = new codi::ExternalFunctionUserData();

    // create the index and double vector for the matrix
    IssmDouble::Real* valueATrans = xNew<IssmDouble::Real>(local_nnz);
    IssmDouble::Identifier* indexATrans = xNew<IssmDouble::Identifier>(local_nnz);
    int* irnATrans = xNew<int>(local_nnz);
    int* jcnATrans = xNew<int>(local_nnz);

    // read the data for the matrix A in a transposed fashion
	  for (int i=0; i<local_nnz; ++i) {
        getPrimalAndGradData(a_loc[i], valueATrans[i], indexATrans[i]);
        irnATrans[i]=jcn_loc[i];  // transposed store
        jcnATrans[i]=irn_loc[i];  // transposed store
    }

    // create the index vector for a (primal values are not needed for a)
    IssmDouble::Identifier* indexB = xNew<IssmDouble::Identifier>(n);
    getVectorGradData(rhs, indexB, n);

    dataHandler->addData(valueATrans);
    dataHandler->addData(indexATrans);
    dataHandler->addData(irnATrans);
    dataHandler->addData(jcnATrans);
    dataHandler->addData(indexB);
  }

  // unpack the primal values from the matrix and the vector
  IssmDouble::Real* valueA = xNew<IssmDouble::Real>(local_nnz);
  IssmDouble::Real* valueB = xNew<IssmDouble::Real>(n);
  // read the data from A and B
  getVectorPrimal(a_loc, valueA, local_nnz);
  getVectorPrimal(rhs, valueB, n);

  MumpsSolve(n, nnz, local_nnz, irn_loc, jcn_loc, valueA, valueB, parameters);
  // valueB contains now the solution

  // pack the values into rhs
  setVectorPrimal(rhs, valueB, n);

  if(tape.isActive()) {
    // create the index vector X and register x as active variables
    IssmDouble::Identifier* indexX = xNew<IssmDouble::Identifier>(n);
    registerVector(rhs, indexX, n);

    dataHandler->addData(valueB); // contains the values from x
    dataHandler->addData(indexX);

    // store other arguments
    dataHandler->addData(n);
    dataHandler->addData(nnz);
    dataHandler->addData(local_nnz);
    dataHandler->addData(parameters); // we assume here that parameters is still intact when the reverse run is called

	 tape.pushExternalFunction(codi::ExternalFunction<Tape>::create(&MumpsSolve_codi_b,(void*)dataHandler, &MumpsSolve_codi_delete));
  }
  else{
    // if the tape is active valueB is stored in the dataHandler and deleted in the reverse sweep
    xDelete(valueB);
  }

  xDelete(valueA);
}
/*}}}*/
#elif _CODIPACK_MAJOR_==1
void MumpsSolve_codi_b(void* tape_in,void* data_in,void* ra) {/*{{{*/

	/*recast data_in and tape*/
	codi::DataStore* data = (codi::DataStore*)data_in;
	//IssmDouble::TapeType& tape = (IssmDouble::TapeType&)tape_in;
	IssmDouble::TapeType& tape = IssmDouble::getGlobalTape();

	IssmDouble::Real* valueATrans;
	IssmDouble::GradientData* indexATrans;
	int* irnATrans;
	int* jcnATrans;
	IssmDouble::GradientData* indexB;
	IssmDouble::Real* valueX;
	IssmDouble::GradientData* indexX;
	int n;
	int nnz;
	int local_nnz;
	Parameters* parameters;

	data->getData(valueATrans);
	data->getData(indexATrans);
	data->getData(irnATrans);
	data->getData(jcnATrans);
	data->getData(indexB);
	data->getData(valueX);
	data->getData(indexX);
	data->getData(n);
	data->getData(nnz);
	data->getData(local_nnz);
	data->getData(parameters);

	// create the adjoint vector for x and reset the adjoint values on the tape
	IssmDouble::GradientValue* adjX = xNew<IssmDouble::GradientValue>(n);
	getVectorAdjoint(tape, indexX, adjX, n);

	MumpsSolve(n, nnz, local_nnz, irnATrans, jcnATrans, valueATrans, adjX, parameters);
	// adjX contains now the solution

	updateVectorAdjoint(tape, indexB, adjX, n);

	// bcast dp_y (the solution of the forward system)
	ISSM_MPI_Bcast(valueX,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
	// bcast the adjoint of the right-hand-side, i.e. this solution
	ISSM_MPI_Bcast(adjX,n,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

	for(int i=0; i<local_nnz; ++i) {
		// we access the transposed matrix here because we stored the indices in a transposed way
		// -1 is substracted because jcn and irn are stored with fortran indexing
		updateAdjoint(tape, indexATrans[i], -adjX[jcnATrans[i]-1]*valueX[irnATrans[i]-1]);
	}

	xDelete(adjX);
}
/*}}}*/
void MumpsSolve_codi_delete(void* tape_in,void* data_in) {/*{{{*/

	/*recast data_in*/
	codi::DataStore* data = (codi::DataStore*)data_in;

	IssmDouble::Real* valueATrans;
	IssmDouble::GradientData* indexATrans;
	int* irnATrans;
	int* jcnATrans;
	IssmDouble::GradientData* indexB;
	IssmDouble::Real* valueX;
	IssmDouble::GradientData* indexX;
	int n;
	int nnz;
	int local_nnz;
	Parameters* parameters;

	data->getData(valueATrans);
	data->getData(indexATrans);
	data->getData(irnATrans);
	data->getData(jcnATrans);
	data->getData(indexB);
	data->getData(valueX);
	data->getData(indexX);
	data->getData(n);
	data->getData(nnz);
	data->getData(local_nnz);
	data->getData(parameters);

	xDelete(valueATrans);
	xDelete(indexATrans);
	xDelete(irnATrans);
	xDelete(jcnATrans);
	xDelete(indexB);
	xDelete(valueX);
	xDelete(indexX);

	delete data;
}
/*}}}*/
void MumpsSolve(int n,int nnz,int local_nnz,int* irn_loc,int* jcn_loc,IssmDouble *a_loc,IssmDouble *rhs,Parameters* parameters){/*{{{*/
	IssmDouble::TapeType& tape = IssmDouble::getGlobalTape();
	codi::DataStore* dataHandler = NULL;

	if(tape.isActive()) {
		dataHandler = new codi::DataStore();

		// create the index and double vector for the matrix
		IssmDouble::Real* valueATrans = xNew<IssmDouble::Real>(local_nnz);
		IssmDouble::GradientData* indexATrans = xNew<IssmDouble::GradientData>(local_nnz);
		int* irnATrans = xNew<int>(local_nnz);
		int* jcnATrans = xNew<int>(local_nnz);

		// read the data for the matrix A in a transposed fashion
		for (int i=0; i<local_nnz; ++i) {
			getPrimalAndGradData(a_loc[i], valueATrans[i], indexATrans[i]);
			irnATrans[i]=jcn_loc[i];  // transposed store
			jcnATrans[i]=irn_loc[i];  // transposed store
		}

		// create the index vector for a (primal values are not needed for a)
		IssmDouble::GradientData* indexB = xNew<IssmDouble::GradientData>(n);
		getVectorGradData(rhs, indexB, n);

		dataHandler->addData(valueATrans);
		dataHandler->addData(indexATrans);
		dataHandler->addData(irnATrans);
		dataHandler->addData(jcnATrans);
		dataHandler->addData(indexB);
	}

	// unpack the primal values from the matrix and the vector
	IssmDouble::Real* valueA = xNew<IssmDouble::Real>(local_nnz);
	IssmDouble::Real* valueB = xNew<IssmDouble::Real>(n);
	// read the data from A and B
	getVectorPrimal(a_loc, valueA, local_nnz);
	getVectorPrimal(rhs, valueB, n);

	MumpsSolve(n, nnz, local_nnz, irn_loc, jcn_loc, valueA, valueB, parameters);
	// valueB contains now the solution

	// pack the values into rhs
	setVectorPrimal(rhs, valueB, n);

	if(tape.isActive()) {
		// create the index vector X and register x as active variables
		IssmDouble::GradientData* indexX = xNew<IssmDouble::GradientData>(n);
		registerVector(rhs, indexX, n);

		dataHandler->addData(valueB); // contains the values from x
		dataHandler->addData(indexX);

		// store other arguments
		dataHandler->addData(n);
		dataHandler->addData(nnz);
		dataHandler->addData(local_nnz);
		dataHandler->addData(parameters); // we assume here that parameters is still intact when the reverse run is called

		//tape.pushExternalFunction(&MumpsSolve_codi_b, dataHandler, &MumpsSolve_codi_delete);
		tape.pushExternalFunctionHandle(&MumpsSolve_codi_b,(void*)dataHandler, &MumpsSolve_codi_delete);
	}
	else{
		// if the tape is active valueB is stored in the dataHandler and deleted in the reverse sweep
		xDelete(valueB);
	}

	xDelete(valueA);
}
/*}}}*/
#else
#error "_CODIPACK_MAJOR_ not supported"
#endif
#endif
