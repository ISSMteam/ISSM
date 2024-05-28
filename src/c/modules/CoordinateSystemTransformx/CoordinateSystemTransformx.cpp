/*!\file CoordinateSystemTransformx
 * \brief: x code for CoordinateSystemTransformx
 */

/*Header files*/
#include "./CoordinateSystemTransformx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include <proj.h>

void CoordinateSystemTransformx(double** px_dest,double** py_dest,double* x_src,double* y_src,int size,const char* str_src,const char* str_dst){

/*
NOTE:
- Compilation of this module is fenced in Makefile.am, so we do not have to 
check again if _HAVE_PROJ_ is defined
*/
#if defined(PROJ_VERSION_MAJOR) && PROJ_VERSION_MAJOR >= 6
	/*
	Converts an array of longitude (x_src) and array of latitude (y_src) in the 
	projection described by str_src to the projection described by str_dest.

	NOTE:
	- API types and calls have been migrated from PROJ 4 to PROJ 6. See SVN 
	revision history for changes.

	Sources:
	- https://proj.org/development/migration.html
	- https://www.gaia-gis.it/fossil/libspatialite/wiki?name=PROJ.6

	TODO:
	- Because transformations are done in place, we could save memory and time 
	by simply passing in x_src and y_src and removing the px_dest and py_dest 
	parameters entirely?
	*/

	/*Allocate output and initialize values as src*/
	_assert_(size>0);
	double* x_dest = xNew<double>(size);
	double* y_dest = xNew<double>(size);

	for(int i=0;i<size;i++){
		x_dest[i] = x_src[i];
		y_dest[i] = y_src[i];
	}

	PJ *P;
	size_t sx = sizeof(double);
	size_t sy = sizeof(double);
	size_t nx = size;
	size_t ny = size;

	P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,str_src,str_dst,NULL);

	if ( 0 == P ) {
		proj_destroy(P);
		_error_("failed to initialize CRS transformation object");
	}

	int p = proj_trans_generic(P, PJ_FWD, x_dest, sx, nx, y_dest, sy, ny, 0, 0, 0, 0, 0, 0);

	if ( 0 == p ){
		proj_destroy(P);
		_error_("no successful transformations");
	}

	proj_destroy(P);

	/*Output : */
	*px_dest=x_dest;
	*py_dest=y_dest;
#else
	_error_("PROJ version >= 6 required");
#endif
}
