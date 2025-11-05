/*!\file: numerics.h
 * \brief prototypes for numerics.h
 */ 

#ifndef _NUMERICS_H_
#define  _NUMERICS_H_

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Verbosity.h"
#include "./GaussPoints.h"
#include "./isnan.h"
#include "./recast.h"
#include "./types.h"
#include "./constants.h"
#include "./OptPars.h"
#include "./interpolation.h"

#if !defined(_HAVE_CODIPACK_)
// already defined in codipack headers
IssmDouble  min(IssmDouble a,IssmDouble b);
IssmDouble  max(IssmDouble a,IssmDouble b);
#endif

#ifdef _HAVE_AD_
IssmPDouble  min(IssmPDouble a,IssmPDouble b);
IssmPDouble  max(IssmPDouble a,IssmPDouble b);
#endif

int         min(int a,int b);
int         max(int a,int b);
void        BrentSearch(IssmDouble** pJ,OptPars optpars,IssmDouble* X0,IssmDouble (*f)(IssmDouble*,void*),IssmDouble (*g)(IssmDouble**,IssmDouble*,void*),void* usr);
void        cross(IssmDouble *result,IssmDouble*vector1,IssmDouble*vector2);
bool        XZvectorsToCoordinateSystem(IssmDouble *T,IssmDouble*xzvectors);
int         cubic(IssmDouble a, IssmDouble b, IssmDouble c, IssmDouble d,IssmDouble X[3], int *num);
IssmDouble  legendre(IssmDouble Pn1, IssmDouble Pn2, IssmDouble x, int n);
IssmDouble*  p_polynomial_value ( int m, int n, IssmDouble* x);

int         NewtonSolveDnorm(IssmDouble* pdnorm,IssmDouble c1,IssmDouble c2,IssmDouble c3,IssmDouble n,IssmDouble dnorm);
IssmDouble  ODE1(IssmDouble alpha,IssmDouble beta,IssmDouble Si, IssmDouble dt,int method);

void LineSectionNormal(IssmDouble* result, IssmDouble* xyz_section);
void TriangleFacetNormal(IssmDouble* normal, IssmDouble* xyz_facet);

/*Interpolation*/
IssmDouble bilinearinterp(IssmDouble* x_grid,IssmDouble* y_grid,IssmDouble* data,IssmDouble x,IssmDouble y,int m,int n,int Nx);
/*templates*/
double triangleinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);
double nearestinterp(double x1,double x2,double y1,double y2,double Q11,double Q12,double Q21,double Q22,double x,double y);

#endif //ifndef _NUMERICS_H_
