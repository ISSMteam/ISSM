/*!\file InterpFromMesh2dx.h
 * \brief: header file for Data interpolation routines.
 */

#ifndef _INTERPFROMMESH2DX_H
#define _INTERPFROMMESH2DX_H

#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

/*threading: */
typedef struct{

	int                 interpolation_type;
	bool                debug;
	int                 nels_data;
	int                *index_data;
	double              *x_data;
	double              *y_data;
	double              *data;
	double              xmin,xmax;
	double              ymin,ymax;
	int                 nods_prime;
	IssmSeqVec<IssmPDouble> *data_prime;
	double              *x_prime;
	double              *y_prime;
	double              *default_values;
	int                 num_default_values;
	double              *incontour;

} InterpFromMesh2dxThreadStruct;

int InterpFromMesh2dx(IssmSeqVec<IssmPDouble>** pdata_prime,int* index_data, double* x_data, double* y_data, int nods_data,int nels_data, double* data, int data_length, double* x_prime, double* y_prime, int nods_prime,
		double* default_values,int num_default_values,Contour<IssmPDouble>** contours,int numcontours);

void* InterpFromMesh2dxt(void* vInterpFromMesh2dxThreadStruct);

#endif /* _INTERPFROMMESH2DX_H */
