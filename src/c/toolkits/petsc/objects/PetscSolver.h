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

void	PetscSolve(PetscVec<IssmDouble>** puf, PetscMat<IssmDouble>* Kff, PetscVec<IssmDouble>* pf, PetscVec<IssmDouble>* uf0,PetscVec<IssmDouble>* df, Parameters* parameters);
void	SolverxPetsc(PVec* puf, PMat Kff, PVec pf, PVec uf0,PVec df, Parameters* parameters);
void  DofTypesToIndexSet(IS* pisv, IS* pisp, PVec df,int typeenum);

#endif //#ifndef _PETSCSOLVER_H_
