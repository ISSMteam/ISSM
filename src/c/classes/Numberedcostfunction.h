/*!\file Numberedcostfunction.h
 * \brief: header file for Numberedcostfunction object
 */

#ifndef _NUMBEREDCOSTFUNCTION_H_
#define _NUMBEREDCOSTFUNCTION_H_

/*Headers:*/
#include "./Definition.h"
#include "./FemModel.h"

class Numberedcostfunction: public Object, public Definition{

	public: 

		int   definitionenum;
		char* name;
		int   number_cost_functions;
		int*  cost_functions_list;

		/*Numberedcostfunction constructors, destructors :*/
		Numberedcostfunction();
		Numberedcostfunction(char* in_name, int in_definitionenum,int number_cost_functions_in,int* cost_functions_list_in);
		~Numberedcostfunction();

		/*Object virtual function resolutoin: */
		Object*	copy();
		void		DeepEcho(void);
		void		Echo(void);
		int		Id(void);
		void		Marshall(MarshallHandle* marshallhandle);
		int		ObjectEnum(void);

		/*Definition virtual function resolutoin: */
		int		DefinitionEnum();
		char*		Name();
		IssmDouble Response(FemModel* femmodel);
};

#endif  /* _NUMBEREDCOSTFUNCTION_H_ */
