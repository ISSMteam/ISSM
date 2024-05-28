#ifndef _CONTAINER_MATERIALS_H_
#define  _CONTAINER_MATERIALS_H_

/*forward declarations */
#include "../../datastructures/datastructures.h"
class Parameters;
class Elements;
class Vertices;
class Loads;
class Nodes;

/*! \brief Declaration of Materials class.
 *
 * Declaration of Materials class.  Materials are vector lists (Containers) of Material objects.
 */ 
class Materials: public DataSet{

	public:

		/*constructors, destructors*/
		Materials();
		~Materials();

		/*numerics*/
		void  Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
		void  ResetHooks();

};

#endif //ifndef _MATERIALS_H_
