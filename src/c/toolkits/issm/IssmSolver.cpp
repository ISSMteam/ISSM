/*!\file IssmSolver
 * \brief implementation of solver 
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <cstring>

#include "../../shared/shared.h"
#include "./issmtoolkit.h"

void IssmSolve(IssmVec<IssmDouble>** pout,IssmMat<IssmDouble>* Kff, IssmVec<IssmDouble>* pf, Parameters* parameters){/*{{{*/

	/*Let matrix decide, to retain object orientation: */
	IssmVec<IssmDouble>* outvector=NULL;

	outvector=Kff->Solve(pf,parameters);

	*pout=outvector;
}
/*}}}*/
