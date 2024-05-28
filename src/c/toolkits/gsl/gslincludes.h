/* \file gslsincludes.h
 * \brief all includes for GSL layer
 */

#ifndef _GSL_INCLUDES_H_
#define _GSL_INCLUDES_H_

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

void DenseGslSolve(IssmPDouble** pX,IssmPDouble* A,IssmPDouble* B, int n);
void DenseGslSolve(IssmDouble** pX,IssmDouble* Kff,int Kff_M,int Kff_N,IssmDouble* pf,int pf_M,Parameters* parameters);

void SolverxSeq(IssmPDouble *X, IssmPDouble *A, IssmPDouble *B,int n);

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
void SolverxSeq(IssmDouble *X,IssmDouble *A,IssmDouble *B,int n, Parameters* parameters);
// call back functions:
ADOLC_ext_fct EDF_for_solverx;
ADOLC_ext_fct_fos_forward EDF_fos_forward_for_solverx;
ADOLC_ext_fct_fos_reverse EDF_fos_reverse_for_solverx;
ADOLC_ext_fct_fov_forward EDF_fov_forward_for_solverx;
ADOLC_ext_fct_fov_reverse EDF_fov_reverse_for_solverx;
#endif

#endif
