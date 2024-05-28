#ifndef _VERTEXONEDGE_H_
#define _VERTEXONEDGE_H_

#include "./include.h"
#include "./Edge.h"

namespace bamg {

	class Mesh;
	class BamgVertex;

	class VertexOnEdge {

		public:
			BamgVertex* v;
			Edge*   be;
			double abcisse;

			//Constructors
			VertexOnEdge(BamgVertex * w, Edge *bw,double s) :v(w),be(bw),abcisse(s) {};
			VertexOnEdge(){};

			//Operators
			operator double () const { return abcisse;}
			operator BamgVertex* () const { return v;}  
			operator Edge* () const { return be;}  
			BamgVertex & operator[](int i) const { return (*be)[i];}

			//Methods
			void SetOnBTh();
			void Set(const Mesh &,long,Mesh &);  
	};

}
#endif
