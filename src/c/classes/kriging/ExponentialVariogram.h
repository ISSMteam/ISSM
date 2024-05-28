/*! \file ExponentialVariogram.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _EXPONENTIALVARIOGRAM_H_
#define _EXPONENTIALVARIOGRAM_H_

/*Headers:*/
#include "./Variogram.h"

class ExponentialVariogram: public Variogram{

	public:
		double nugget; //The height of the jump of the semivariogram at the discontinuity at the origin
		double sill;   //Limit of the variogram tending to infinity lag distances
		double range;  //The distance in which the difference of the variogram from the sill becomes negligible

		/*ExponentialVariogram constructors, destructors*/
		ExponentialVariogram();
		ExponentialVariogram(Options* options);
		~ExponentialVariogram();

		/*Object virtual functions definitions*/
		Object* copy();
		void  DeepEcho(){_error_("Not implemented yet");};
		void  Echo();
		int   Id(){_error_("Not implemented yet");}; 
		void Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!"); };
		int   ObjectEnum(){_error_("Not implemented yet");};

		/*Variogram functions*/
		double Covariance(double deltax,double deltay);
		double SemiVariogram(double deltax,double deltay);
};
#endif  /* _EXPONENTIALVARIOGRAM_H */
