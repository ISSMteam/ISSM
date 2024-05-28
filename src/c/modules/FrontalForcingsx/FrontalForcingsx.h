#ifndef _FRONTALFORCINGSX_H
#define _FRONTALFORCINGSX_H

#include "../../classes/classes.h"
#include "../../analyses/analyses.h"

/* local prototypes: */
void FrontalForcingsx(FemModel* femmodel);
void Subglacialdischargearmax(FemModel* femmodel);
void Thermalforcingarmax(FemModel* femmodel);

#endif
