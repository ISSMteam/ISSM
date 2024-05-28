/*!\file Kriging.h
 * \brief: header file for Kriging
 */

#ifndef _KRIGINGX_H
#define _KRIGINGX_H

#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

class Observations;
class Variogram;

int  Krigingx(double** ppredictions,double **perror,double* x, double* y, double* observations, int n_obs,double* x_interp,double* y_interp,int n_interp,Options* options);
int  pKrigingx(double** ppredictions,double **perror,double* x, double* y, double* observations, int n_obs,double* x_interp,double* y_interp,int n_interp,Options* options);
void ProcessVariogram(Variogram **pvariogram,Options* options);
void ProcessVariogram2(Variogram **pvariogram,Options* options);

/*threading: */
typedef struct{
	int           n_interp;
	double       *x_interp;
	double       *y_interp;
	double        radius;
	int           mindata;
	int           maxdata;
	Variogram    *variogram;
	Observations *observations;
	double       *predictions;
	double       *error;
	int          *numdone;
	double        power;//for idw
}KrigingxThreadStruct;

void* Krigingxt(void*);
void* NearestNeighbort(void*);
void* idwt(void*);
void* v4t(void*);
void* Distancest(void*);
#endif /* _KRIGINGX_H */
