/*!\file:  BamgConvertMeshx.h
 * \brief header file for Bamg module
 */ 

#ifndef _BAMGCONVERTMESHX_H
#define _BAMGCONVERTMESHX_H

#include "../../classes/classes.h"
#include "../../bamg/bamgobjects.h"

/* local prototypes: */
int BamgConvertMeshx(BamgMesh* bamgmesh,BamgGeom* bamggeom,int* index,double* x,double* y,int nods,int nels);

#endif
