/*!\file Cfdragcoeffabsgradtransient.h
 * \brief: header file for Cfdragcoeffabsgradtransient object
 */

#ifndef _CFDRAGCOEFFABSGRADTRANSIENT_H_
#define _CFDRAGCOEFFABSGRADTRANSIENT_H_

/*Headers:*/
#include "./Definition.h"
class FemModel;

class Cfdragcoeffabsgradtransient: public Object, public Definition{

	public: 

		int         definitionenum;
		char       *name;
		int         num_datatimes;
		IssmDouble *datatimes;
		bool       *passedflags;
		IssmDouble  J;

		/*Cfdragcoeffabsgradtransient constructors, destructors :*/
		Cfdragcoeffabsgradtransient();
		Cfdragcoeffabsgradtransient(char* in_name, int in_definitionenum, int num_datatimes, IssmDouble* in_datatime);
		Cfdragcoeffabsgradtransient(char* in_name, int in_definitionenum, int num_datatimes, IssmDouble* in_datatime, bool* in_timepassedflag, IssmDouble in_J);
		~Cfdragcoeffabsgradtransient();

		/*Object virtual function resolutoin: */
		Object *copy();
		void    DeepEcho(void);
		void    Echo(void);
		int     Id(void);
		void    Marshall(MarshallHandle  *marshallhandle);
		int     ObjectEnum(void);

		/*Definition virtual function resolutoin: */
		int         DefinitionEnum();
		char       *Name();
		IssmDouble  Response(FemModel                       *femmodel);
		IssmDouble  Cfdragcoeffabsgradtransient_Calculation(Element  *element);
};
#endif  /* _CFDRAGCOEFFABSGRADTRANSIENT_H_ */
