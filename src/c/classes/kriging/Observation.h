/*! \file Observation.h 
 *  \brief: header file for Observation object
 */

#ifndef _OBSERVATION_H_
#define _OBSERVATION_H_

#include "../../datastructures/datastructures.h"

class Observation: public Object{

	public:
		double x,y;
		int    xi,yi;
		int    index;
		double weight;
		double value;

		/*Observation constructors, destructors*/
		Observation();
		Observation(double x_in,double y_in,int xi_in,int yi_in,int index_in,double value_in);
		Observation(double x_in,double y_in,double value_in);
		~Observation();

		/*Object virtual functions definitions*/
		bool operator==(const Observation& ob) const;
		Object *copy();
		void    DeepEcho()  {_error_("Not implemented yet"); };
		double  distance(const Observation& ob) const;
		void    Echo();
		int     Id()        {_error_("Not implemented yet"); };
		void    print() const;
		void    Marshall(MarshallHandle* marshallhandle){ _error_("not implemented yet!");};
		int     ObjectEnum(){_error_("Not implemented yet"); };

		/*Management*/
		void WriteXYObs(const Observation& ob, double* px, double* py, double* pobs);
		void WriteXYObs(double* px,double* py,double* pobs);
};
#endif  /* _OBSERVATION_*/
