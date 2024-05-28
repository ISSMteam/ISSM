/*!\file:  InputUpdateFromConstantx.h
 * \brief header file for updating datasets from inputs
 */ 

#ifndef _UPDATEINPUTSFROMCONSTANTXX_H
#define _UPDATEINPUTSFROMCONSTANTXX_H

#include "../../classes/classes.h"
class Inputs;

/* local prototypes: */
void InputUpdateFromConstantx(FemModel* femmodel,bool       constant,int name);
void InputUpdateFromConstantx(FemModel* femmodel,int        constant,int name);
void InputUpdateFromConstantx(FemModel* femmodel,int        constant,int name, int type);
void InputUpdateFromConstantx(FemModel* femmodel,IssmDouble constant,int name);
#ifdef _HAVE_AD_
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmPDouble constant,int name);
#endif
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmDouble constant,int name);
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmDouble constant,int name, int type);
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,bool       constant,int name);

#endif  /* _UPDATEINPUTSFROMCONSTANTXX_H */
