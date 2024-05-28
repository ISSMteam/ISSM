/* \file toolkits.h
 * \brief: this API allows use of external packages, provides patches, etc ...
 */

#ifndef _TOOLKITS_H_
#define _TOOLKITS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#ifdef _HAVE_PETSC_
#include "./petsc/petscincludes.h"
#endif

#include "./mpi/issmmpi.h"

#ifdef _HAVE_METIS_
#include "./metis/metisincludes.h"
#endif

#ifdef _HAVE_GSL_
#include "./gsl/gslincludes.h"
#endif

#ifdef _HAVE_ADOLC_
#include "./adolc/adolcincludes.h"
#endif

#ifdef _HAVE_CODIPACK_
#include "./codipack/codipackincludes.h"
#endif

#ifdef _HAVE_TRIANGLE_
#include "./triangle/triangleincludes.h"
#endif

#include "./objects/toolkitobjects.h"
#include "./toolkitsenums.h"
#include "./issm/issmtoolkit.h"
#include "./ToolkitOptions.h"
#endif
