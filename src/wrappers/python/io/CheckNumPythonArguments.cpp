/*!\file CheckNumPythonArguments.cpp:
 * \brief: check number of arguments and report an usage error message.
 */

#define PY_ARRAY_UNIQUE_SYMBOL PythonIOSymbol
#define NO_IMPORT

#include "./pythonio.h"
#include "../../c/shared/Exceptions/exceptions.h"

int CheckNumPythonArguments(PyObject* inputs,int NRHS, void (*function)( void )){

	Py_ssize_t size=0;

	/*figure out size of tuple in input: */
	size=PyTuple_Size(inputs);

	/*check on requested size: */
	if (size==0){
		function();
		_error_("usage: see above");
	}
	else if (size!=NRHS ) {
		function(); 
		_error_("usage error.");
	}
	return 1;
}
