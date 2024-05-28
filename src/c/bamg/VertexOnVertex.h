#ifndef _VERTEXONVERTEX_H_
#define _VERTEXONVERTEX_H_

#include "./include.h"
#include "./BamgVertex.h"

namespace bamg {

	class Mesh;

	class VertexOnVertex {

		public:
			BamgVertex* v;
			BamgVertex* bv;

			//Constructors
			VertexOnVertex();
			VertexOnVertex(BamgVertex * w,BamgVertex *bw);

			//Methods
			void SetOnBTh();
			void Set(const Mesh &,long,Mesh &);
	};

}
#endif
