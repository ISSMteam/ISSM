/*!\file:  InterpFromMeshToMesh2dx.h
 * \brief header file for Bamg module
 */ 

#ifndef _INTERPFROMMESHTOMESH2DX_H
#define _INTERPFROMMESHTOMESH2DX_H

#include "../../classes/classes.h"

int InterpFromMeshToMesh2dx(double** pdata_interp,int* index_data,double* x_data,double* y_data,int nods_data,int nels_data,
			double* data,int M_data,int N_data,double* x_interp,double* y_interp,int N_interp,Options* options=NULL);

#endif
