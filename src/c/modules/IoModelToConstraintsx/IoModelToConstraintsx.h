/*!\file:  IoModelToConstraintsx.h
 */ 

#ifndef _IOMODEL_TO_CONSTRAINTS_H_
#define _IOMODEL_TO_CONSTRAINTS_H_

#include "../../classes/classes.h"

/* local prototypes: */
void IoModelToConstraintsx(Constraints* constraints,IoModel* iomodel,const char* spc_name,int analysis_type,int finite_element,int dof=0);
void IoModelToConstraintsx(Constraints* constraints,IoModel* iomodel,IssmDouble* spcdata,int M,int N,int analysis_type,int finite_element,int dof=0);
void IoModelToDynamicConstraintsx(Constraints* constraints,IoModel* iomodel,const char* spc_name,int analysis_type,int finite_element,int dof=0);
void IoModelToDynamicConstraintsx(Constraints* constraints,IoModel* iomodel,IssmDouble* spcdata,int M,int N,int analysis_type,int finite_element,int dof=0);

#endif  /* _IOMODELTOELEMENTINPUTX_H */
