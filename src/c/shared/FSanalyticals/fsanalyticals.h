/*!\file: elements.h
 * \brief prototypes for elements.h
 */ 

#ifndef _SHARED_ANALYTICALS_H_
#define _SHARED_ANALYTICALS_H_

#include "../Numerics/types.h"

IssmDouble fx(IssmDouble x_coord, IssmDouble y_coord, IssmDouble z_coord, int testid);
IssmDouble fy(IssmDouble x_coord, IssmDouble y_coord, IssmDouble z_coord, int testid);
IssmDouble fz(IssmDouble x_coord, IssmDouble y_coord, IssmDouble z_coord, int testid);
IssmDouble alpha(IssmDouble x_coord, IssmDouble y_coord, IssmDouble z_coord, int testid);

#endif //ifndef _SHARED_ANALYTICALS_H_
