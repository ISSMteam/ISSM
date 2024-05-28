/*!\file Cfdragcoeffabsgrad.h
 * \brief: header file for Cfdragcoeffabsgrad object
 */

#ifndef _CFDRAGCOEFFABSGRAD_H_
#define _CFDRAGCOEFFABSGRAD_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Cfdragcoeffabsgrad: public Object, public Definition{

	public: 

		int         definitionenum;
		char       *name;
		bool			firsttimepassed;
		IssmDouble  J;

		/*Cfdragcoeffabsgrad constructors, destructors :*/
		Cfdragcoeffabsgrad();
		Cfdragcoeffabsgrad(char* in_name, int in_definitionenum);
		Cfdragcoeffabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J);
		Cfdragcoeffabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J, bool in_firsttimepassed);
		~Cfdragcoeffabsgrad();

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
		IssmDouble  Cfdragcoeffabsgrad_Calculation(Element  *element);
};
#endif  /* _CFDRAGCOEFFABSGRAD_H_ */
