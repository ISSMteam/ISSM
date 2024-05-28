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
bool   findindices(int* pn,int* pm,double* x,int x_rows, double* y,int y_rows, double xgrid,double ygrid);
double triangleinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);
double bilinearinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);
double nearestinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);

#endif /* _INTERPFROMGRIDTOMESHX_H */
