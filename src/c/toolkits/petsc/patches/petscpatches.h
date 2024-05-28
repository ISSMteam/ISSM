/*\file petscpatches.h
 * \brief: our own patches for petsc use
 */

#ifndef _PETSC_PATCHES_H_
#define _PETSC_PATCHES_H_

#include <petscksp.h>

#include "./SolverEnum.h"
#include "../../toolkitsenums.h"
#include "../../../shared/io/Comm/IssmComm.h"

class Parameters;

Vec NewVec(int size,ISSM_MPI_Comm comm,bool fromlocalsize=false);
Mat NewMat(int M,int N,ISSM_MPI_Comm comm);
Mat NewMat(int M,int N,double sparsity,ISSM_MPI_Comm comm);
Mat NewMat(int M,int N,int connectivity,int numberofdofspernode, ISSM_MPI_Comm comm);

int VecToMPISerial(double** pgathered_vector, Vec vector,ISSM_MPI_Comm comm,bool broadcast=true);
void MatFree(Mat* pmat);
void ISFree(IS* pis);
void VecFree(Vec* pvec);
void KSPFree(KSP* pksp);
int MatPartition(Mat* poutmatrix,Mat matrixA,double* row_partition_vector,int row_partition_vector_size, double* col_partition_vector,int col_partition_vector_size);
void PetscOptionsDetermineSolverType(int* psolver_type);
void MatMultPatch(Mat A,Vec X, Vec AX,ISSM_MPI_Comm comm);
void MatToMPISerial(double** poutmatrix,Mat matrix,ISSM_MPI_Comm comm,bool broadcast=true);
Vec  SerialToVec(double* vector,int vector_size);
InsertMode ISSMToPetscInsertMode(InsMode mode);
NormType ISSMToPetscNormMode(NormMode mode);
MatType ISSMToPetscMatrixType(MatrixType type);

#endif
