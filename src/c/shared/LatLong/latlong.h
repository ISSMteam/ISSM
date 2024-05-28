/*!\file: latlong.h
 * \brief prototypes for latlong.h
 */ 

#ifndef _SHARED_LATLONG_H_
#define _SHARED_LATLONG_H_

int Ll2xyx(double* x, double* y, double* lat, double* lon, int ncoord, int sgn);
int Ll2xyx(double* x, double* y, double* lat, double* lon, int ncoord, int sgn, double central_meridian, double standard_parallel);
void Ll2xydef(double* pdelta, double* pslat, int sgn);

int Xy2llx(double* lat, double* lon, double* x, double* y, int ncoord, int sgn);
int Xy2llx(double* lat, double* lon, double* x, double* y, int ncoord, int sgn, double central_meridian, double standard_parallel);
void Xy2lldef(double* pdelta, double* pslat, int sgn);

#endif //ifndef _SHARED_LATLONG_H_
