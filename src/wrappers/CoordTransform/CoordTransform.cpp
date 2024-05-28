/*
 * CoordTransform: mesh a domain using an .exp file
 */

/*Header files*/
#include "./CoordTransform.h"
#include <proj.h>

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
	*/

	void CoordTransformUsage(void){/*{{{*/
		_printf_(" type help CoordTransform\n");
	}/*}}}*/
	WRAPPER(CoordTransform_python){
		/*intermediary: */
		double *xin     = NULL;
		double *yin     = NULL;
		char   *projin  = NULL;
		char   *projout = NULL;
		int     M,N;
		int     test1,test2;

		/*Boot module: */
		MODULEBOOT();

		/*checks on arguments: */
		CHECKARGUMENTS(NLHS,NRHS,&CoordTransformUsage);

		/*Fetch data needed for meshing: */
		FetchData(&xin,&M,&N,XIN);
		FetchData(&yin,&test1,&test2,YIN);
		if(!M*N)     _error_("no coordinate provided");
		if(test1!=M) _error_("x and y do not have the same size");
		if(test2!=N) _error_("x and y do not have the same size");
		FetchData(&projin,PROJIN);
		FetchData(&projout,PROJOUT);

		/*Calculate size once*/
		int size = M*N;

		/*Initialize output*/
		double* xout = xNew<double>(size);
		double* yout = xNew<double>(size);

		for(int i=0;i<size;i++){
			xout[i] = xin[i];
			yout[i] = yin[i];
		}

		size_t sx = sizeof(double);
		size_t sy = sizeof(double);
		size_t nx = size;
		size_t ny = size;

		PJ* P = proj_create_crs_to_crs(PJ_DEFAULT_CTX,projin,projout,NULL);

		if(P==0){
			proj_destroy(P);
			_error_("Projection string not recognized");
		}

		int p = proj_trans_generic(P, PJ_FWD, xout, sx, nx, yout, sy, ny, 0, 0, 0, 0, 0, 0);

		if(p==0){
			proj_destroy(P);
			_error_("projection failed");
		}

		/*Cleanup*/
		proj_destroy(P);

		/*write outputs: */
		WriteData(XOUT,xout,M,N);
		WriteData(YOUT,yout,M,N);

		/*Clean-up and return*/
		xDelete<double>(xin);
		xDelete<double>(yin);
		xDelete<double>(xout);
		xDelete<double>(yout);
		xDelete<char>(projin);
		xDelete<char>(projout);

		/*end module: */
		MODULEEND();
	}
#else
	_error_("PROJ version >= 6 required");
#endif
