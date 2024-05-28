/*!\file:  Variogram.h
 * \brief abstract class for Variogram object
 */ 

#ifndef _VARIOGRAM_H_
#define _VARIOGRAM_H_

#include "../../datastructures/datastructures.h"

class Variogram: public Object{

	public: 
		virtual ~Variogram(){};
		virtual double Covariance(double deltax,double deltay)=0;
		virtual double SemiVariogram(double deltax,double deltay)=0;

};
#endif
