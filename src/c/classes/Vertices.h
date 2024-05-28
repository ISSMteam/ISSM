#ifndef _CONTAINER_VERTICES_H_
#define  _CONTAINER_VERTICES_H_

/*forward declarations */
#include "../datastructures/datastructures.h"
#include "../shared/shared.h"
class IoModel;

/*!\brief Declaration of Vertices class.
 *
 * Declaration of Vertices class.  Vertices are vector lists (Containers) of Vertex objects.
 * A vertex is a set of (x,y,z) coordinates defining the location of points in the mesh (not
 * to be confused with a node, which is a degree of freedom (DOF) for a particular analysis).
 */ 
class Vertices: public DataSet{

	private:
		int numberofvertices;
		int numberofvertices_local;
		int numberofmasters_local;
	public:
		int*  common_recv;
		int** common_recv_ids;
		int*  common_send;
		int** common_send_ids;

		/*constructors, destructors:*/ 
		Vertices();
		~Vertices();

		/*Objects virtual functions*/
		Vertices* Copy();
		void      Marshall(MarshallHandle* marshallhandle);

		/*numerics:*/
		void  Finalize(IoModel* iomodel);
		int   NumberOfVertices(void);
		int   NumberOfVerticesLocal(void);
		int   NumberOfVerticesLocalAll(void);
		void  LatLonList(IssmDouble** lat,IssmDouble** lon);
		void  XYList(IssmDouble** pxcoords,IssmDouble** pycoords);
};

#endif //ifndef _VERTICES_H_
