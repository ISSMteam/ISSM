/*!\file: DependentObject.h
 * \brief prototype for DependentObject.h
 */ 

#ifndef _DEPENDENTOBJECT_H_
#define  _DEPENDENTOBJECT_H_

/*{{{*/
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
/*}}}*/

class FemModel;

class DependentObject: public Object{

	public:

		char* name;
		IssmDouble response_value;

		/*DependentObject constructors, destructors */
		DependentObject();
		DependentObject(char* name);
		DependentObject(char* name, IssmDouble in_response);
		~DependentObject();

		/*Object virtual functions definitions*/
		Object *copy(void);
		void    DeepEcho();
		void    Echo();
		int     Id();
		int     ObjectEnum();
		void    Marshall(MarshallHandle  *marshallhandle);

		/*DependentObject methods: */
		void       RecordResponsex(FemModel*femmodel);
		IssmDouble GetValue(void);
		void       ResetResponseValue(void);

};
#endif //ifndef _DEPENDENTOBJECT_H_
