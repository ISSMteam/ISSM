/*!\file:  InputUpdateFromVectorDakotax.h
 * \brief header file for updating datasets from inputs
 */ 

#ifndef _UPDATEINPUTSFROMVECTORDAKOTAXX_H
#define _UPDATEINPUTSFROMVECTORDAKOTAXX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	InputUpdateFromVectorDakotax(FemModel* femmodel,Vector<IssmDouble>* vector, int name,int type);
void	InputUpdateFromVectorDakotax(FemModel* femmodel,IssmDouble* vector, int name,int type);

#endif  /* _UPDATEINPUTSFROMVECTORDAKOTAXX_H */
