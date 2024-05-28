/*!\file Cfsurfacesquaretransient.h
 * \brief: header file for Cfsurfacesquaretransient object
 */

#ifndef _CFSURFACESQUARETRANSIENT_H_
#define _CFSURFACESQUARETRANSIENT_H_

/*Headers:*/
#include "./Definition.h"
class FemModel;

class Cfsurfacesquaretransient: public Object, public Definition{

	public: 

		int         definitionenum;
		int         model_enum;
		char       *name;
		int         num_datatimes;
		IssmDouble *datatimes;
		bool       *passedflags;
		IssmDouble  J;

		/*Cfsurfacesquaretransient constructors, destructors :*/
		Cfsurfacesquaretransient();
		Cfsurfacesquaretransient(char* in_name, int in_definitionenum, int in_model_enum,int num_datatimes, IssmDouble* in_datatime);
		Cfsurfacesquaretransient(char* in_name, int in_definitionenum, int in_model_enum,int num_datatimes, IssmDouble* in_datatime, bool* in_timepassedflag, IssmDouble in_J);
		~Cfsurfacesquaretransient();

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
		IssmDouble  Response(FemModel *femmodel);
		IssmDouble  Cfsurfacesquaretransient_Calculation(Element  *element, int model_enum);
};
#endif  /* _CFSURFACESQUARE_H_ */
