#ifndef SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_
#define SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_CODIPACK_)

#include <codi.hpp>

#ifndef CODIPACK_TYPE
#define CODIPACK_TYPE codi::RealReverse
#endif

using CoDiReal = CODIPACK_TYPE;

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_ */
