/*!\file Zgesvx
 * \brief: solves double complex linear equation system A*x=b 
 */

#include "./Zgesvx.h"

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"

extern "C" { 
	int zgesv_(int* pnyi, int* pnrhs, IssmComplex* yilocal, int* plda, int* ipiv, IssmComplex* rhslocal, int* pldb, int* pinfo);
}

void Zgesvx(int* pnyi, int* pnrhs, IssmComplex* yilocal, int* plda, int* ipiv, IssmComplex* rhslocal, int* pldb, int* pinfo){

	/*Call fortran driver: */
	zgesv_(pnyi, pnrhs, yilocal, plda, ipiv, rhslocal, pldb, pinfo);

}
