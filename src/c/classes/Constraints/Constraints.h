#ifndef _CONTAINER_CONSTRAINTS_H_
#define  _CONTAINER_CONSTRAINTS_H_

/*forward declarations */
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"

/*! \brief Declaration of Constraints class. 
 *
 * Declaration of Constraints class for handling Single Point Constraints (SPCs).
 * Constraints are vector lists (Containers) of Constraint objects.
 */ 
class Constraints: public DataSet{

	public:

		/*Object constructors and destructor*/
		Constraints(){/*{{{*/
			enum_type=ConstraintsEnum;
			return;
		}
		/*}}}*/
		~Constraints(){/*{{{*/
			return;
		}
		/*}}}*/

		/*numerics*/
		void ActivatePenaltyMethod(int in_analysis);
};

#endif //ifndef _CONSTRAINTS_H_
