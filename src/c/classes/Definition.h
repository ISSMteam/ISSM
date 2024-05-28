/*!\file: Definition.h: abstract class used by some objects to behave like response objects 
 * that can be called according to a string (the output definition)
 */ 

#ifndef _DEFINITION_H_
#define  _DEFINITION_H_

/*Headers:*/
class FemModel;
class Definition {

	public:
		virtual       ~Definition(){};
		virtual char*  Name()=0;
		virtual int    DefinitionEnum()=0;
		virtual IssmDouble  Response(FemModel*)=0;

};

#endif //ifndef _DEFINITION_H_
