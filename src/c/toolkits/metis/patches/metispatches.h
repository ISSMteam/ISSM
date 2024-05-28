/*\file metispatches.h
 * \brief: our own patches for metis. Mainly to work through new apis from 4.0 to 5.0 version.
 */

#ifndef _METIS_PATCHES_H_
#define _METIS_PATCHES_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

void METIS_PartMeshNodalPatch(int numberofelements,int numberofnodes,int* index,int* vweights,int num_procs,int* epart,int* npart);

#endif
