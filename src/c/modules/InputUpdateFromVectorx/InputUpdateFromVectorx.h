/*!\file:  InputUpdateFromVectorx.h
 * \brief header file for updating datasets from inputs
 */ 

#ifndef _UPDATEINPUTSFROMVECTORXX_H
#define _UPDATEINPUTSFROMVECTORXX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	InputUpdateFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector, int name,int type);
void	InputUpdateFromVectorx(FemModel* femmodel,IssmDouble* vector, int name,int type);

#endif  /* _UPDATEINPUTSFROMVECTORXX_H */
