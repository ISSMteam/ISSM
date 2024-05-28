/*!\file:  FloatingiceMeltingRatex.h
 * \brief header file for Floatingice melting rate
 */ 

#ifndef _FloatingiceMeltingRatex_H
#define _FloatingiceMeltingRatex_H

#include "../../classes/classes.h"
#include "../FloatingiceMeltingRatePicox/FloatingiceMeltingRatePicox.h"

/* local prototypes: */
void FloatingiceMeltingRatex(FemModel* femmodel);

void LinearFloatingiceMeltingRatex(FemModel* femmodel);
void SpatialLinearFloatingiceMeltingRatex(FemModel* femmodel);
void MismipFloatingiceMeltingRatex(FemModel* femmodel);
void FloatingiceMeltingRateIsmip6x(FemModel* femmodel);
void BeckmannGoosseFloatingiceMeltingRatex(FemModel* femmodel);
void LinearFloatingiceMeltingRatearmax(FemModel* femmodel);

#endif  /* _FloatingiceMeltingRatex_H*/
