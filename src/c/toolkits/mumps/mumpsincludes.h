/* \file mumpsincludes.h
 * \brief all includes for MUMPS layer
 */

#ifndef _MUMPS_INCLUDES_H_
#define _MUMPS_INCLUDES_H_

/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/Numerics/types.h"
/*}}}*/

class Parameters;
template <class doubletype> class SparseRow;

#ifdef _HAVE_MPI_
void MpiDenseMumpsSolve(IssmDouble* uf,int uf_M,int uf_n, IssmDouble* Kff,int Kff_M, int Kff_N, int Kff_m, IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters);
void MpiSparseMumpsSolve(IssmDouble* uf,int uf_M,int uf_n, SparseRow<IssmDouble>** Kff,int Kff_M, int Kff_N, int Kff_m, IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters);
#endif
void SeqDenseMumpsSolve(IssmDouble* uf,int uf_M,int uf_n, IssmDouble* Kff,int Kff_M, int Kff_N, int Kff_m, IssmDouble* pf, int pf_M, int pf_m, Parameters* parameters);

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
// call back functions:
ADOLC_ext_fct_iArr mumpsSolveEDF;
ADOLC_ext_fct_iArr_fos_reverse fos_reverse_mumpsSolveEDF;
ADOLC_ext_fct_iArr_fov_reverse fov_reverse_mumpsSolveEDF;
#endif

#endif
