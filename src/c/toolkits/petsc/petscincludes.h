/* \file petscincludes.h
 * \brief all includes from Petsc + our own patches
 */

#ifndef _PETSC_INCLUDES_H_
#define _PETSC_INCLUDES_H_

#ifdef _HAVE_PETSC_

/*Petsc includes: */
#include <petscksp.h>

/*our own patches: */
#include "patches/petscpatches.h"
#include "objects/petscobjects.h"

#endif
#endif
