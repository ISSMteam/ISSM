#ifndef _CONTAINER_RESULTS_H_
#define _CONTAINER_RESULTS_H_

#include "../../datastructures/datastructures.h"

/*forward declarations */
class Parameters;
class ExternalResult;

/*!\brief Declaration of Results class.
 *
 * Declaration of Results class.  Results are vector lists (Containers) of Result objects.
 */ 
class Results: public DataSet{

	public:

		/*constructors, destructors*/
		Results();
		~Results();

		/*Mehthos*/
		int AddResult(ExternalResult* result);
		int DeleteResult(int result_enum,int result_step);
		ExternalResult* FindResult(int result_enum);
		void Write(Parameters* parameters);
};
#endif //ifndef _RESULTS_H_
