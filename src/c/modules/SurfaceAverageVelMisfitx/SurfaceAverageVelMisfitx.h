/*!\file:  SurfaceAverageVelMisfitx.h
 * \brief header file for inverse methods misfit computation
 */ 

#ifndef _SURFACEAVERAGEVELMISFITX_H
#define _SURFACEAVERAGEVELMISFITX_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceAverageVelMisfitx(IssmDouble* pJ,FemModel* femmodel);
IssmDouble SurfaceAverageVelMisfit(Element* element);

#endif
