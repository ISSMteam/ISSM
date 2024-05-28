/*!\file:  PetscMat.h
 */ 

#ifndef _PETSC_SOLVER_H_
#define _PETSC_SOLVER_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../petscincludes.h"
class Parameters;

/*}}}*/

void	PetscSolve(PetscVec** puf, PetscMat* Kff, PetscVec* pf, PetscVec* uf0,PetscVec* df, Parameters* parameters);
void	SolverxPetsc(Vec* puf, Mat Kff, Vec pf, Vec uf0,Vec df, Parameters* parameters);
void    DofTypesToIndexSet(IS* pisv, IS* pisp, Vec df,int typeenum);

#endif //#ifndef _PETSCSOLVER_H_
