/*!\file Cfsurfacelogvel.h
 * \brief: header file for Cfsurfacelogvel object
 */

#ifndef _CFSURFACELOGVEL_H_
#define _CFSURFACELOGVEL_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Cfsurfacelogvel: public Object, public Definition{

	public: 

		int         definitionenum;
		char*       name;
		IssmDouble	datatime;
		bool			timepassedflag;
		IssmDouble  J;

		/*Cfsurfacelogvel constructors, destructors :*/
		Cfsurfacelogvel();
		Cfsurfacelogvel(char* in_name, int in_definitionenum, IssmDouble in_datatime);
		Cfsurfacelogvel(char* in_name, int in_definitionenum, IssmDouble in_datatime, bool timepassedflag, IssmDouble in_J);
		~Cfsurfacelogvel();

		/*Object virtual function resolutoin: */
		Object *copy();
		void    DeepEcho(void);
		void    Echo(void);
		int     Id(void);
		void    Marshall(MarshallHandle  *marshallhandle);
		int     ObjectEnum(void);

		/*Definition virtual function resolutoin: */
		int DefinitionEnum();
		char* Name();
		IssmDouble Response(FemModel* femmodel);
		IssmDouble Cfsurfacelogvel_Calculation(Element* element, int definitionenum);
};
#endif  /* _CFSURFACELOGVEL_H_ */
