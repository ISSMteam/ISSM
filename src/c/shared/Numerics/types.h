/*!\file: types.h
 * \brief prototypes for types.h
 */ 

#ifndef _TYPES_H_
#define  _TYPES_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <complex>

/*here are our abstracted types: inspired on petsc */
#if ISSM_USE_64BIT_INDICES == 1
typedef long long IssmInt;
#else
typedef int IssmInt;
#endif  

#if defined(_HAVE_ADOLC_) &&  !defined(_WRAPPERS_)
/*ADOLC typedefs*/
#include "adolc/adolc.h"
typedef adouble              IssmDouble;  /*for active variables*/
//typedef acomplex             IssmComplex; /*for active variables*/ /*FIXME!*/
typedef std::complex<double> IssmComplex; /*for active variables*/ /*FIXME!*/
typedef double               IssmPDouble; /*for passive variables*/
typedef std::complex<double> IssmPComplex;/*for passive variables*/

#elif defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
/*CoDiPack typedefs*/
#include <codi.hpp>
//typedef codi::RealReverseIndex          IssmDouble;
typedef codi::RealReverse               IssmDouble;
typedef std::complex<codi::RealReverse> IssmComplex;
typedef double                          IssmPDouble;
typedef IssmComplex                     IssmPComplex;

#else 
/*Non-AD typedefs*/
typedef double               IssmDouble; 
typedef std::complex<double> IssmComplex; 
typedef IssmDouble           IssmPDouble;
typedef IssmComplex          IssmPComplex;
#endif

#endif //ifndef _TYPES_H_
