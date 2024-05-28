#ifndef _VERTEXONGEOM_H_
#define _VERTEXONGEOM_H_

#include "./include.h"
#include "./GeomVertex.h"

namespace bamg {

	class Mesh;
	class BamgVertex;
	class GeomEdge;

	class VertexOnGeom{

		public:

			BamgVertex* meshvertex;
			double curvilincoord;  
			union{ 
				GeomVertex* gv; // if curvilincoord <0; 
				GeomEdge*   ge; // if curvilincoord in [0..1]
			};

			//Constructors/Destructors
			VertexOnGeom();
			VertexOnGeom(BamgVertex & m,GeomVertex &g);
			VertexOnGeom(BamgVertex & m,GeomEdge &g,double s);

			//Operators
			operator BamgVertex*() const  {return meshvertex;}
			operator GeomVertex * () const  {return gv;}
			operator GeomEdge * () const  {return ge;}
			operator const double & () const {return curvilincoord;}

			//Methods
			int  OnGeomVertex()const;
			int  OnGeomEdge() const;
			int  IsRequiredVertex();
			void SetOn();

			//Inline methods
			void Set(const VertexOnGeom&,const Mesh &,Mesh &);  
	};
}
#endif
