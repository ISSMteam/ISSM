/*!\file: matrix.h
 * \brief prototypes for matrix.h
 */ 

#ifndef _MATRIXUTILS_H_
#define _MATRIXUTILS_H_

#include "../Numerics/types.h"

int  TripleMultiply( IssmDouble* a, int nrowa, int ncola, int itrna, IssmDouble* b, int nrowb, int ncolb, int itrnb, IssmDouble* c, int nrowc, int ncolc, int itrnc, IssmDouble* d, int iaddd);
int  MatrixMultiply( IssmDouble* a, int nrowa, int ncola, int itrna, IssmDouble* b, int nrowb, int ncolb, int itrnb, IssmDouble* c, int iaddc );
int  MatrixInverse( IssmDouble* a, int ndim, int nrow, IssmDouble* b, int nvec, IssmDouble* pdet );

void Matrix2x2Invert(IssmDouble* Ainv, IssmDouble* A);
void Matrix2x2Determinant(IssmDouble* Adet,IssmDouble* A);
void Matrix2x2Eigen(IssmDouble* plambda1,IssmDouble* plambda2,IssmDouble* pvx, IssmDouble* pvy,IssmDouble a11, IssmDouble a21,IssmDouble a22);

void Matrix3x3Invert(IssmDouble* Ainv, IssmDouble* A);
void Matrix3x3Determinant(IssmDouble* Adet,IssmDouble* A);
IssmDouble Matrix3x3Determinant( IssmDouble a1,IssmDouble a2,IssmDouble a3, IssmDouble b1,IssmDouble b2,IssmDouble b3, IssmDouble c1,IssmDouble c2,IssmDouble c3);
void Matrix3x3Solve(IssmDouble* X,IssmDouble* A,IssmDouble* B);

void Matrix4x4Adjoint(IssmDouble* Aadj, IssmDouble* A);
void Matrix4x4Invert(IssmDouble* Ainv, IssmDouble* A);
void Matrix4x4Determinant(IssmDouble* Adet,IssmDouble* A);
void Matrix4x4Solve(IssmDouble* X,IssmDouble* A,IssmDouble* B);

void newcell(IssmDouble** pcell, IssmDouble newvalue, bool top, int m);
IssmDouble  cellsum(IssmDouble* cell, int m);
void celldelete(IssmDouble** pcell, int m, int* indices, int nind);
void cellsplit(IssmDouble** pcell, int m, int i,IssmDouble scale);
void cellecho(int numcells, int m, ...);
void CholeskyRealPositiveDefinite(IssmDouble* Lchol, IssmDouble* A, int ndim);
#endif //ifndef _MATRIXUTILS_H_
