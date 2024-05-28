#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "./include.h"
#include "./BamgGeom.h"
#include "./BamgOpts.h"
#include "./GeomVertex.h"
#include "./GeomEdge.h"
#include "./Curve.h"

namespace bamg {

	class Triangle;
	class BamgQuadtree;
	class GeomSubDomain;
	class Edge;

	class Geometry { 

		public:

			long           NbRef;                 // counter of ref on the this class if 0 we can delete
			long           nbv;                   // number of vertices
			long           nbe;                   // number of edges
			long           nbsubdomains;
			long           nbcurves;
			GeomVertex    *vertices;
			GeomEdge      *edges;
			BamgQuadtree  *quadtree;
			GeomSubDomain *subdomains;
			Curve         *curves;
			R2             pmin,pmax;             // domain extrema coordinates
			double         coefIcoor;             // coef to integer coordinates;

			//Constructor/Destructors
			~Geometry(); 
			Geometry();
			Geometry(const Geometry & Gh);
			Geometry(BamgGeom* bamggeom, BamgOpts* bamgopts);

			//Operators
			const GeomVertex &operator[](long i) const { return vertices[i]; };
			GeomVertex       &operator[](long i) { return vertices[i];       };
			const GeomEdge   &operator()(long i) const { return edges[i];    };
			GeomEdge         &operator()(long  i) { return edges[i];         };

			//Methods
			void             Echo();
			I2               R2ToI2(const R2 &P) const;
			double           MinimalHmin();
			double           MaximalHmax();
			void             ReadGeometry(BamgGeom *bamggeom, BamgOpts*bamgopts);
			void             Init(void);
			void             PostRead(bool checkcurve=false);
			long             GetId(const GeomVertex &t) const;
			long             GetId(const GeomVertex *t) const;
			long             GetId(const GeomEdge &t) const;
			long             GetId(const GeomEdge *t) const;
			long             GetId(const Curve *c) const;
			void             UnMarkEdges();
			GeomEdge        *ProjectOnCurve(const Edge &,double,BamgVertex &,VertexOnGeom &) const;
			void             WriteGeometry(BamgGeom *bamggeom, BamgOpts*bamgopts);
	};

}
#endif
