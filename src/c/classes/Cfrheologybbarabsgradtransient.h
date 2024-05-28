/*!\file Cfrheologybbarabsgradtransient.h
 * \brief: header file for Cfrheologybbarabsgradtransient object
 */

#ifndef _CFRHEOLOGYBBARABSGRADTRANSIENT_H_
#define _CFRHEOLOGYBBARABSGRADTRANSIENT_H_

/*Headers:*/
#include "./Definition.h"
class FemModel;

class Cfrheologybbarabsgradtransient: public Object, public Definition{

	public: 

		int         definitionenum;
		char       *name;
		int         num_datatimes;
		IssmDouble *datatimes;
		bool       *passedflags;
		IssmDouble  J;

		/*Cfrheologybbarabsgradtransient constructors, destructors :*/
		Cfrheologybbarabsgradtransient();
		Cfrheologybbarabsgradtransient(char* in_name, int in_definitionenum, int num_datatimes, IssmDouble* in_datatime);
		Cfrheologybbarabsgradtransient(char* in_name, int in_definitionenum, int num_datatimes, IssmDouble* in_datatime, bool* in_timepassedflag, IssmDouble in_J);
		~Cfrheologybbarabsgradtransient();

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
		IssmDouble Cfrheologybbarabsgradtransient_Calculation(Element* element);
};
#endif  /* _CFRHEOLOGYBBARABSGRADTRANSIENT_H_ */
