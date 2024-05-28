/*!\file:  KMLMeshWritex.h
 * \brief header file for kml mesh writer routines.
 */ 

#ifndef _KMLMESHWRITEX_H
#define _KMLMESHWRITEX_H

#include <float.h>    /*  DBL_MAX  */
#include "../../kml/kmlobjects.h"

/* local prototypes: */
void KMLMeshWritex(int* ierror,
				   char* name,
				   char* notes,
				   int* elem,int melem,int nelem,
				   int* nodecon,int mncon,int nncon,
				   double* lat, double* lng,
				   int* part,
				   double* data, int mdata, int ndata,
				   double* cmap, int mcmap, int ncmap,
				   FILE* fid);

KML_Folder* KMLMeshElem(int* elem,int melem,int nelem,
						int* nodecon,int mncon,int nncon,
						double* lat, double* lng,
						double* edata,
						double* cmap, int mcmap, int ncmap);

#endif  /* _KMLMESHWRITEX_H */
