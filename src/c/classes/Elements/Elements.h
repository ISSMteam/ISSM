#ifndef _CONTAINER_ELEMENTS_H_
#define  _CONTAINER_ELEMENTS_H_

/*forward declarations */
#include "../../datastructures/datastructures.h"
class Materials;
class Parameters;
class Vertices;
class Loads;
class Nodes;
class Results;

/*! \brief Declaration of Elements class 
 *
 * Declaration of Elements class.  Elements are vector lists (Containers) of Element objects.
 */ 
class Elements: public DataSet{

	public:

		/*constructors, destructors*/
		Elements();
		~Elements();

		/*numerics*/
		void   Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters,Inputs* inputs);
		int    MaxNumNodes(void);
		int    NumberOfElements(void);
		void   SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
		void   ResetHooks();
};

#endif //ifndef _ELEMENTS_H_
