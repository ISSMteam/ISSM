/*! \file Vertex.h 
 *  \brief: header file for vertex object
 */

#ifndef _VERTEX_H_
#define _VERTEX_H_

/*Headers:*/
/*{{{*/
#include "./classes.h"
#include "../shared/Exceptions/exceptions.h"
#include "../toolkits/toolkits.h"
template <class doubletype> class Vector;
template <class doubletype> class Matrix;
class Parameters;
class IoModel;
/*}}}*/

class Vertex: public Object{

	public: 
		bool       clone;
		int        domaintype;
		int        id;           // random index
		int        sid;          // "serial" id (rank of this vertex if the dataset was on 1 cpu)
		int        pid;          // "parallel" id
		int        lid;          // "local" id
		IssmDouble x;
		IssmDouble y;
		IssmDouble z;
		IssmDouble latitute;
		IssmDouble longitude;
		IssmDouble R;
		IssmDouble sigma;        //sigma coordinate: (z-bed)/thickness
		int        connectivity; //number of vertices connected to this vertex

		/*Vertex constructors, destructors {{{*/
		Vertex();
		Vertex(int id, int sid,bool clone, IoModel* iomodel,bool isamr);
		~Vertex();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		void  Echo();
		void  DeepEcho();
		int   Id(); 
		int   ObjectEnum();
		Object* copy();
		void Marshall(MarshallHandle* marshallhandle);

		/*}}}*/
		/*Vertex management:*/ 
		int        Connectivity(void); 
		IssmDouble GetLatitude(void); 
		IssmDouble GetLongitude(void); 
		IssmDouble GetRadius(void); 
		IssmDouble GetX(void); 
		IssmDouble GetY(void); 
		IssmDouble GetZ(void); 
		int        Pid(void); 
		int        Lid(void); 
		int        Sid(void); 
		void       UpdatePosition(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz,Parameters* parameters,IssmDouble* thickness,IssmDouble* bed);
		void       VertexCoordinates(Vector<IssmDouble>* vx,Vector<IssmDouble>* vy,Vector<IssmDouble>* vz,bool spherical=false);
};

/*Methods relating to Vertex object: */
void GetVerticesCoordinates(IssmDouble* xyz,Vertex** vertices, int numvertices,bool spherical=false);

#endif  /* _VERTEX_H */
