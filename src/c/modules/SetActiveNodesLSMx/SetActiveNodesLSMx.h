/*!\file:  SetActiveNodesLSMx.h
 * \brief header file for updating single point constraints  for next time step
 */ 

#ifndef _SETACTIVENODESLSMX_H
#define _SETACTIVENODESLSMX_H

#include "../../classes/classes.h"

void SetActiveNodesLSMx(FemModel* femmodel,bool ishydrology=false,bool isdebris=false);
void GetMaskOfIceVerticesLSMx0(FemModel* femmodel,bool ishydrology=false,bool isdebris=false);
void GetMaskOfIceVerticesLSMx(FemModel* femmodel,bool ishydrology=false,bool isdebris=false);
#endif  /* _SETACTIVENODESLSMX_H*/
