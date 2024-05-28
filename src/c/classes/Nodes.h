#ifndef _CONTAINER_NODES_H_
#define  _CONTAINER_NODES_H_

#include "../datastructures/datastructures.h"
#include "../toolkits/toolkits.h"
class Parameters;
class Elements;
class Vertices;
class Loads;
class Nodes;
class Materials;
class MarshallHandle;

/*!\brief Declaration of Nodes class.
 *
 * Declaration of Nodes class.  Nodes are vector lists of objects (Containers) of Node objects.
 * Node objects are the degrees of freedom (DOFs) for a particular analysis type (not to be 
 * confused with a vertex, which defines the (x,y,z) location of a point).
 */ 
class Nodes: public DataSet{

	private:
		int numberofnodes;
		int numberofnodes_local;
		int numberofmasters_local;
	public:
		int*  common_recv;
		int** common_recv_ids;
		int*  common_send;
		int** common_send_ids;

		/*constructors, destructors*/
		Nodes();
		~Nodes();

		/*Objects virtual functions*/
		Nodes* Copy();
		void   Marshall(MarshallHandle* marshallhandle);

		/*numerics*/
		void  DistributeDofs(int SETENUM);
		void  Finalize(void);
		int   MaxNumDofs(int setenum);
		int   NumberOfDofs(int setenum);
		int   NumberOfDofsLocal(int setenum);
		int   NumberOfDofsLocalAll(int setenum);
		int   NumberOfNodes(void);
		int   NumberOfNodesLocal(void);
		int   NumberOfNodesLocalAll(void);
		bool  RequiresDofReindexing(void);
		void  CheckDofListAcrossPartitions(void);
		void  GetLocalVectorWithClonesGset(IssmDouble** plocal_vector,Vector<IssmDouble> *vector);
};

#endif //ifndef _NODES_H_
