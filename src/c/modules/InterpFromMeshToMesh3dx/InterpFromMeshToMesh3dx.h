/*!\file InterpFromMeshToMesh3dx.h
 * \brief: header file for Data interpolation routines.
 */

#ifndef _INTERPFROMMESHTOMESH3DX_H
#define _INTERPFROMMESHTOMESH3DX_H

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"

int InterpFromMeshToMesh3dx(IssmSeqVec<IssmPDouble>** pdata_prime,double* index_data, double* x_data, double* y_data, double* z_data, int nods_data,int nels_data, double* data, int data_length, double* x_prime, double* y_prime, double* z_prime, int nods_prime,double default_value);

#endif /* _INTERPFROMMESHTOMESH3DX_H */
