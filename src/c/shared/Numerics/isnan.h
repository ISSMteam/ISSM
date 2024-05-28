#ifndef _XISNAN_H_
#define _XISNAN_H_

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*cmath defines isnan and isinf*/
#include <cmath>

#ifdef _INTEL_WIN_
template <class T> int xIsNan(const T& X){return (X!=X)?1:0;}
#else
template <class T> int xIsNan(const T& X){return std::isnan(X); }
#endif
template <class T> int xIsInf(const T& X){return std::isinf(X); }

/*Special overloading definitions for AD*/

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
#include "./types.h"
template <> int xIsNan<adouble> (const adouble& X);
template <> int xIsInf<adouble> (const adouble& X);
#endif

#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
#include "./types.h"
template <> int xIsNan<IssmDouble> (const IssmDouble& X);
template <> int xIsInf<IssmDouble> (const IssmDouble& X);
#endif

#endif
