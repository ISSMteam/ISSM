/*!\file:  InputUpdateFromDakotax.h
 * \brief header file for updating datasets from inputs
 */ 

#ifndef _INPUTUPDATEFROMDAKOTAXX_H
#define _INPUTUPDATEFROMDAKOTAXX_H

#include "../../classes/classes.h"

void  InputUpdateFromDakotax(FemModel* femmodel,double* variables,char* *variables_descriptors,int numvariables);
void  InputUpdateSpecialtyCode(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* variable_partition,int npart,char* root);
void  MmeToInput(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* variable_partition,int npart,int rootenum, int interpolationenum);
void InputScaleFromDakotax(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* partition, int npart, int nt, int name);

#endif  /* _INPUTUPDATEFROMDAKOTAXX_H */
