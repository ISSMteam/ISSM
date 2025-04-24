/*!\file Cfsurfacesquare.h
 * \brief: header file for Cfsurfacesquare object
 */

#ifndef _CFSURFACESQUARE_H_
#define _CFSURFACESQUARE_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Cfsurfacesquare: public Object, public Definition{

	public: 

		int         definitionenum;
		int         surfaceid;
		int         model_enum;
		char       *name;
		IssmDouble  datatime;
		bool        timepassedflag;
		IssmDouble  J;

		/*Cfsurfacesquare constructors, destructors :*/
		Cfsurfacesquare();
		Cfsurfacesquare(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime, int in_surfaceid);
		Cfsurfacesquare(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime, bool timepassedflag, IssmDouble in_J,int in_surfaceid);
		~Cfsurfacesquare();

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
		IssmDouble  Response(FemModel                    *femmodel);
		IssmDouble  Cfsurfacesquare_Calculation(Element  *element, int model_enum);
};
#endif  /* _CFSURFACESQUARE_H_ */
