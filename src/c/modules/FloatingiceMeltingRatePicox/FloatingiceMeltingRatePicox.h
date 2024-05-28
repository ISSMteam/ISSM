/*!\file:  FloatingiceMeltingRatePicox.h
 * \brief header file for Floatingice melting rate
 */ 

#ifndef _FloatingiceMeltingRatePicox_H
#define _FloatingiceMeltingRatePicox_H

#include "../../classes/classes.h"

/* local prototypes: */
void FloatingiceMeltingRatePicox(FemModel* femmodel);

void UpdateBoxIdsPico(FemModel* femmodel);
void ComputeBoxAreasPico(FemModel* femmodel);
void UpdateBoxPico(FemModel* femmodel, int loopboxid);
void ComputeAverageOceanvarsPico(FemModel* femmodel, int boxid);
void ComputeBasalMeltPlume(FemModel* femmodel);

#endif  /* _FloatingiceMeltingRatePicox_H*/
