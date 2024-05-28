/*!\file Nodalvalue.h
 * \brief: header file for Nodalvalue object
 */

#ifndef _NODALVALUE_H_
#define _NODALVALUE_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "./FemModel.h"
/*}}}*/

void NodalValuex( IssmDouble* pnodalvalue, int natureofdataenum,Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters);

class Nodalvalue: public Object, public Definition{

	public: 

		int         definitionenum;
		int         model_enum;
		char*       name;
		int         node;

		/*Nodalvalue constructors, destructors :*/
		Nodalvalue();
		Nodalvalue(char* in_name, int in_definitionenum, int in_model_enum, int in_node);
		~Nodalvalue();

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
};

#endif  /* _NODALVALUE_H_ */
