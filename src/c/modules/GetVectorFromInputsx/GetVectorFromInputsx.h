/*!\file:  GetVectorFromInputsx.h
 */

#ifndef _GETVECTORFROMINPUTSXX_H
#define _GETVECTORFROMINPUTSXX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	GetVectorFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type);
void  GetVectorFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type,IssmDouble time);
void	GetVectorFromInputsx(IssmDouble** pvector,FemModel* femmodel,int name,int type);
void  GetVectorFromInputsx(IssmDouble** pvector,int* pvector_size, FemModel* femmodel,int name);

void	GetVectoronBaseFromInputsx(IssmDouble** pvector,FemModel* femmodel,int name,int type);
void	GetVectoronBaseFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type);

#endif  /* _GETVECTORFROMINPUTSXX_H */
