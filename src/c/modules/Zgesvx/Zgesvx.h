/*!\file:  Zgesvx.h
 * \brief header file for ...
 */ 

#ifndef _ZGESVX_H
#define _ZGESVX_H

#include "../../classes/classes.h"

/* local prototypes: */
void Zgesvx(int* pnyi, int* pnrhs, IssmComplex* yilocal, int* plda, int* ipiv, IssmComplex* rhslocal, int* pldb, int* pinfo);

#endif  /* _ZGESVX_H */
