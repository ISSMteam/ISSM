/* \file issmtoolkit.h
 * \brief all includes for MPI layer
 */

#ifndef _ISSM_TOOLKIT_H_
#define _ISSM_TOOLKIT_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./IssmAbsMat.h"
#include "./IssmAbsVec.h"
#include "./IssmDenseMat.h"
#include "./IssmMat.h"
#include "./IssmSeqVec.h"
#include "./IssmVec.h"
#include "./IssmSolver.h"

#ifdef _HAVE_MPI_
#include "./IssmMpiDenseMat.h"
#include "./IssmMpiSparseMat.h"
#include "./IssmMpiVec.h"
#endif

#endif
