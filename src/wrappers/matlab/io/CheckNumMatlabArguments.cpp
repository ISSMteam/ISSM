/*!\file CheckNumMatlabArguments.cpp:
 * \brief: check number of arguments and report an usage error message.
 */

#include "./matlabio.h"

int CheckNumMatlabArguments(int nlhs,int NLHS, int nrhs,int NRHS, const char* THISFUNCTION, void (*function)( void )){

	/*checks on arguments on the matlab side: */
	if (nrhs==0 && nlhs==0) {
		/*unless NLHS=0 and NRHS=0, we are just asking for documentation: */
		if (NRHS==0 && NLHS==0)return 1;
		/* special case: */
		function();
		_error_("usage: see above");
	}
	else if (nlhs!=NLHS || nrhs!=NRHS ) {
		function(); 
		_error_("usage error.");
	}
	return 1;
}
