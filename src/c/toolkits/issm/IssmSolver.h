/*!\file:  IssmSolver.h
 * \brief main hook up from Solver toolkit object to the ISSM toolkit
 */ 

#ifndef _ISSM_SOLVER_H_
#define _ISSM_SOLVER_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/Numerics/types.h"

/*}}}*/

template <class doubletype> class IssmVec;
template <class doubletype> class IssmMat;
class Parameters;

void IssmSolve(IssmVec<IssmDouble>** puf,IssmMat<IssmDouble>* Kff, IssmVec<IssmDouble>* pf,Parameters* parameters);

#endif 
