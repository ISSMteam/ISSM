/*!\file:  UpdateDynamicConstraintsx.h
 * \brief header file for updating single point constraints  for next time step
 */ 

#ifndef _UPDATEDYNAMICCONSTRAINTSXX_H
#define _UPDATEDYNAMICCONSTRAINTSXX_H

#include "../../classes/classes.h"

void UpdateDynamicConstraintsx(Constraints* constraints,Nodes* nodes,Parameters* parameters,Vector<IssmDouble>* yg);

#endif  /* _UPDATESPCSX_H */
