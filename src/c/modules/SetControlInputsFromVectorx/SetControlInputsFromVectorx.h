/*!\file:  SetControlInputsFromVectorx.h
 */ 

#ifndef _SETCONTROLINPUTSXFROMVECTOR_H
#define _SETCONTROLINPUTSXFROMVECTOR_H

#include "../../classes/classes.h"

/* local prototypes: */
void SetControlInputsFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector);
void SetControlInputsFromVectorx(FemModel* femmodel,IssmDouble* vector);

#ifdef _HAVE_AD_
void SetControlInputsFromVectorx(FemModel* femmodel,IssmPDouble* vector);
#endif

#endif 
