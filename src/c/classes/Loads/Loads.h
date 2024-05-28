#ifndef _CONTAINER_LOADS_H_
#define  _CONTAINER_LOADS_H_

/*forward declarations */
#include "../../datastructures/datastructures.h"
class Materials;
class Parameters;
class Elements;
class Vertices;
class Nodes;

/*!\brief Declaration of Loads class.
 *
 * Declaration of Loads class.  Loads are vector lists (Containers) of Load objects.
 */ 
class Loads: public DataSet{

	public:

		int numrifts;
		int numpenalties;

		/*constructors, destructors*/
		Loads();
		~Loads();

		/*Objects virtual functions*/
		Loads* Copy();
		void   Marshall(MarshallHandle* marshallhandle);

		/*numerics*/
		void  Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
		bool  IsPenalty();
		void  Finalize();
		int   MaxNumNodes();
		int   NumberOfLoads();
		void  ResetHooks();
		void  SetCurrentConfiguration(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters);
};

#endif //ifndef _LOADS_H_
