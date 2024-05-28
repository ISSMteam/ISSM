/*!\file:  AllocateSystemMatricesx.h
*/ 

#ifndef _ALLOCATESYSTEMMATRICESX_H
#define _ALLOCATESYSTEMMATRICESX_H

#include "../../classes/classes.h"

/* local prototypes: */
void AllocateSystemMatricesx(Matrix<IssmDouble>** pKff,Matrix<IssmDouble>** pKfs,Vector<IssmDouble>** pdf,Vector<IssmDouble>** ppf,FemModel* femmodel);
void MatrixNonzeros(int** pd_nnz,int** po_nnz,FemModel* femmodel,int set1enum,int set2enum);

#endif
