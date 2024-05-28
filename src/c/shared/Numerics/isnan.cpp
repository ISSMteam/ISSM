#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Special overloading definitions for AD*/
#include "./isnan.h"

#if defined(_HAVE_ADOLC_) && !defined(_WRAPPERS_)
template <> int xIsNan<adouble> (const adouble& X){ return std::isnan(X.getValue()); }
template <> int xIsInf<adouble> (const adouble& X){ return std::isinf(X.getValue()); }
#endif

#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
template <> int xIsNan<IssmDouble> (const IssmDouble& X){ return std::isnan(X.getValue()); }
template <> int xIsInf<IssmDouble> (const IssmDouble& X){ return std::isinf(X.getValue()); }
#endif
