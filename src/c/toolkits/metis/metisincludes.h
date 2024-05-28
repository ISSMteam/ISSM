/* \file metisincludes.h
 * \brief all includes from metis + our own patches
 */

#ifndef _METIS_INCLUDES_H_
#define _METIS_INCLUDES_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

/*Metis includes: */
#include <metis.h>

/*our own patches: */
#include "patches/metispatches.h"

#endif
