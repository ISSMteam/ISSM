/*!\file:  GetSolutionFromInputsx.h
 * \brief header file for updating datasets from inputs
 */

#ifndef _GETSOLUTIONFROMINPUTSXX_H
#define _GETSOLUTIONFROMINPUTSXX_H

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

/* local prototypes: */
void GetSolutionFromInputsx(Vector<IssmDouble>** psolution,FemModel* femmodel);
void GetBasalSolutionFromInputsx(Vector<IssmDouble>** psolution,FemModel* femmodel);

#endif  /* _GETSOLUTIONFROMINPUTSXX_H */
