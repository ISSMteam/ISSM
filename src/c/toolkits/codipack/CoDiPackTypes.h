/*
 * CoDiPackTypes.h
 *
 *  Created on: Jul 20, 2016
 *      Author: a_h_ck
 */

#ifndef SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_
#define SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#if defined(_HAVE_CODIPACK_)

#include <codi.hpp>

using CoDiReal = codi::RealReverse;

#endif /* _HAVE_CODIPACK_ */
#endif /* SRC_C_TOOLKITS_CODIPACK_CODIPACKTYPES_H_ */
