/*! \file PowerVariogram.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _POWERVARIOGRAM_H_
#define _POWERVARIOGRAM_H_

/*Headers:*/
#include "./Variogram.h"

class PowerVariogram: public Variogram{

	public:
		double nugget; //The height of the jump of the semivariogram at the discontinuity at the origin
		double slope;  
		double power; 

		/*PowerVariogram constructors, destructors*/
		PowerVariogram();
		PowerVariogram(Options* options);
		~PowerVariogram();

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
#endif  /* _POWERVARIOGRAM_H */
