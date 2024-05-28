/*!\file Misfit.h
 * \brief: header file for Misfit object
 */

#ifndef _MISFIT_H_
#define _MISFIT_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Misfit: public Object, public Definition{

	public: 

		int         definitionenum;
		int         local;     
		int         model_enum;
		char*       name;
		int         observation_enum;
		char*       timeinterpolation;
		int         weights_enum;

		int         lock; // if lock is on, we just return the value stored in "misfit".  this is used so we don't compute misfit past the final_time
		IssmDouble  misfit; //value carried over in time.

		/*Misfit constructors, destructors :*/
		Misfit();
		Misfit(char* in_name, int in_definitionenum, int in_model_enum, int in_observation_enum, char* in_timeinterpolation, int in_local, int in_weights_enum);
		~Misfit();

		/*Object virtual function resolutoin: */
		Object* copy();
		void DeepEcho(void);
		void Echo(void);
		int Id(void);
		void Marshall(MarshallHandle* marshallhandle);
		int ObjectEnum(void);

		/*Definition virtual function resolutoin: */
		int DefinitionEnum();
		char* Name();
		IssmDouble Response(FemModel* femmodel);
};
#endif  /* _MISFIT_H_ */
