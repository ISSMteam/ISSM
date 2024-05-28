/*!\file Cflevelsetmisfit.h
 * \brief: header file for Cflevelsetmisfit object
 */

#ifndef _CFLEVELSETMISFIT_H_
#define _CFLEVELSETMISFIT_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Cflevelsetmisfit: public Object, public Definition{

	public: 

		int         definitionenum;
		int         model_enum;
		char       *name;
		IssmDouble  datatime;
		bool        timepassedflag;
		IssmDouble  J;

		/*Cflevelsetmisfit constructors, destructors :*/
		Cflevelsetmisfit();
		Cflevelsetmisfit(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime);
		Cflevelsetmisfit(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime, bool timepassedflag, IssmDouble in_J);
		~Cflevelsetmisfit();

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
		IssmDouble Cflevelsetmisfit_Calculation(Element* element, int model_enum);
		IssmDouble Heaviside(IssmDouble x);
};
#endif  /* _CFLEVELSETMISFIT_H_ */
