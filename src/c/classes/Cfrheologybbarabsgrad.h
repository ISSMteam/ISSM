/*!\file Cfrheologybbarabsgrad.h
 * \brief: header file for Cfrheologybbarabsgrad object
 */

#ifndef _CFRHEOLOGYBBARABSGRAD_H_
#define _CFRHEOLOGYBBARABSGRAD_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Cfrheologybbarabsgrad: public Object, public Definition{

	public: 

		int         definitionenum;
		char*       name;
		bool			firsttimepassed;
		IssmDouble  J;

		/*Cfrheologybbarabsgrad constructors, destructors :*/
		Cfrheologybbarabsgrad();
		Cfrheologybbarabsgrad(char* in_name, int in_definitionenum);
		Cfrheologybbarabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J);
		Cfrheologybbarabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J, bool in_firsttimepassed);
		~Cfrheologybbarabsgrad();

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
		IssmDouble Cfrheologybbarabsgrad_Calculation(Element* element);
};
#endif  /* _CFRHEOLOGYBBARABSGRAD_H_ */
