/*!\file:  CreateJacobianMatrixx.h
*/ 

#ifndef _CREATEJACOBIANMATRIXX_H
#define _CREATEJACOBIANMATRIXX_H

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

/* local prototypes: */
void CreateJacobianMatrixx(Matrix<IssmDouble>** pJff,FemModel* femmodel,IssmDouble kmax);

#endif
