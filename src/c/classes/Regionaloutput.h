/*!\file Regionaloutput.h
 * \brief: header file for Regionaloutput object
 */

#ifndef _REGIONALOUTPUT_H_
#define _REGIONALOUTPUT_H_

/*Headers:*/
/*{{{*/
#include "./Definition.h"
#include "./FemModel.h"

/*}}}*/

class Regionaloutput: public Object, public Definition{

public: 

	int         definitionenum;
	char*       outputname;
	char*       name;
	IssmDouble* mask;
	int         M;

	/*Regionalicevolume: constructors, destructors :*/
	Regionaloutput();
	Regionaloutput(char* in_name, int in_definitionenum, char* in_outputname, IssmDouble* maskin, int Min);
	~Regionaloutput();

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
};
#endif  /* _REGIONALOUTPUT_H_ */
