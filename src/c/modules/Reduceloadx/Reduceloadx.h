/*!\file:  Reduceloadx.h
 * \brief reduce loads (wring out boundary conditions)
 */ 

#ifndef _REDUCELOADX_H
#define _REDUCELOADX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	Reduceloadx( Vector<IssmDouble>* pf, Matrix<IssmDouble>* Kfs, Vector<IssmDouble>* ys,bool flag_ys0=false);

#endif  /* _REDUCELOADX_H */
