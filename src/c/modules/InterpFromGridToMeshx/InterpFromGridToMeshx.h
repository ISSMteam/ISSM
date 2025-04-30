/*!\file InterpFromGridToMeshx.h
 * \brief: header file for Data interpolation routines.
 */

#ifndef _INTERPFROMGRIDTOMESHX_H
#define _INTERPFROMGRIDTOMESHX_H

#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"
#include "../../shared/shared.h"

/*threading: */
typedef struct{
	double*             x;
	int                 x_rows;
	double*             y;
	int                 y_rows;
	double*             data;
	double              default_value;
	const char*         interp;
	int                 M;
	int                 N;
	int                 nods;
	double*             x_mesh;
	double*             y_mesh;
	IssmSeqVec<IssmPDouble>* data_mesh;
} InterpFromGridToMeshxThreadStruct;

int    InterpFromGridToMeshx(IssmSeqVec<IssmPDouble>** pdata_mesh,double* x, int x_rows, double* y, int y_rows, double* data, int M, int N, double* x_mesh, double* y_mesh, int nods, double default_value, const char* interptype);
void*  InterpFromGridToMeshxt(void* vInterpFromGridToMeshxThreadStruct);
#endif /* _INTERPFROMGRIDTOMESHX_H */
