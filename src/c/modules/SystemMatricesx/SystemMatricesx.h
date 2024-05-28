/*!\file:  SystemMatricesx.h
*/ 

#ifndef _SYSTEMMATRICESX_H
#define _SYSTEMMATRICESX_H

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

/* local prototypes: */
void SystemMatricesx(Matrix<IssmDouble>** pKff, Matrix<IssmDouble>** pKfs, Vector<IssmDouble>** ppf, Vector<IssmDouble>** pdf, IssmDouble* pkmax,FemModel* femmodel, bool isAllocated=false);

#endif  /* _SYSTEMMATRICESX_H */
