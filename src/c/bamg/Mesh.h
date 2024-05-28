#ifndef _MESH_H_
#define _MESH_H_

#include "./include.h"
#include "./BamgOpts.h"
#include "./BamgMesh.h"
#include "./BamgGeom.h"
#include "./Triangle.h"
#include "./VertexOnGeom.h"
#include "./VertexOnVertex.h"
#include "./VertexOnEdge.h"
#include "./ListofIntersectionTriangles.h"

namespace bamg {

	class Geometry;
	class CrackedEdge;
	class BamgQuadtree;
	class SubDomain;

	class Mesh {

		public:

			Geometry                    & Gh;                    // Geometry
			Mesh                        & BTh;                   // Background Mesh Bth== *this =>no background
			BamgVertex                   *vertices;
			Triangle                     *triangles;
			Edge                         *edges;
			BamgQuadtree                 *quadtree;
			BamgVertex                  **orderedvertices;
			SubDomain                    *subdomains;
			long                          NbRef;                 // counter of ref on the this class if 0 we can delete
			long                          maxnbv,maxnbt;         // nombre max de sommets , de triangles
			long                          nbv,nbt,nbe;           // nb of vertices, of triangles and edges
			long                          nbsubdomains;
			long                          nbtout;                // Nb of oudeside triangle

			R2                            pmin,pmax;             // extrema
			double                        coefIcoor;             // coef to integer
			ListofIntersectionTriangles   lIntTria;
			long									randomseed;            //used for random number generation

			long                          NbVerticesOnGeomVertex;
			VertexOnGeom                 *VerticesOnGeomVertex;
			long                          NbVerticesOnGeomEdge;
			VertexOnGeom                 *VerticesOnGeomEdge;
			long                          NbVertexOnBThVertex;
			VertexOnVertex               *VertexOnBThVertex;
			long                          NbVertexOnBThEdge;
			VertexOnEdge                 *VertexOnBThEdge;
			long                          NbCrackedVertices;
			long                         *CrackedVertices;
			long                          NbCrackedEdges;
			CrackedEdge                  *CrackedEdges;

			//Constructors/Destructors
			Mesh(BamgGeom* bamggeom,BamgMesh* bamgmesh,BamgOpts* bamgopts);
			Mesh(int* index,double* x,double* y,int nods,int nels,BamgOpts* bamgopts);/*MeshConvert*/
			Mesh(double* x,double* y,int nods,BamgOpts* bamgopts); /*BamgTriangulate*/
			Mesh(Mesh &,Geometry * pGh=0,Mesh* pBTh=0,long maxnbv_in=0 ); //copy operator
			Mesh(const Mesh &,const int *flag,const int *bb,BamgOpts* bamgopts); // truncature
			Mesh(long maxnbv,Mesh & BT,BamgOpts* bamgopts,int keepBackVertices=1);
			Mesh(long maxnbv,Geometry & G,BamgOpts* bamgopts);
			~Mesh(); 

			//Operators
			const BamgVertex &operator[](long i) const { return vertices[i];  };
			BamgVertex       &operator[](long i) { return vertices[i];        };
			const Triangle   &operator()(long i) const { return triangles[i]; };
			Triangle         &operator()(long  i) { return triangles[i];             };

			//Methods
			void SetIntCoor(const char * from =0);
			double MinimalHmin();
			double MaximalHmax();
			I2 R2ToI2(const R2 & P) const;
			R2 I2ToR2(const I2 & P) const;
			void AddVertex(BamgVertex & s,Triangle * t,long long *  =0) ;
			void Insert(BamgOpts* bamgopts);
			void Echo(void);
			void ForceBoundary(BamgOpts* bamgopts);
			void FindSubDomain(BamgOpts* bamgopts,int OutSide=0);
			long TriangleReferenceList(long*) const;
			void TriangleIntNumbering(long* renumbering);
			void CrackMesh(BamgOpts* bamgopts);
			void SmoothMetric(BamgOpts* bamgopts,double raisonmax) ;
			void BoundAnisotropy(BamgOpts* bamgopts,double anisomax,double hminaniso= 1e-100) ;
			Edge** MakeGeomEdgeToEdge();
			long SplitInternalEdgeWithBorderVertices();
			void MakeBamgQuadtree();
			void MaxSubDivision(BamgOpts* bamgopts,double maxsubdiv);
			void NewPoints(Mesh &,BamgOpts* bamgopts,int KeepVertices=1);
			long InsertNewPoints(long nbvold,long & NbTSwap,BamgOpts* bamgopts); 
			void TrianglesRenumberBySubDomain(bool justcompress=false);
			void SmoothingVertex(BamgOpts* bamgopts,int =3,double=0.3);
			Metric MetricAt (const R2 &);
			GeomEdge* ProjectOnCurve( Edge & AB, BamgVertex &  A, BamgVertex & B,double theta, BamgVertex & R,VertexOnEdge & BR,VertexOnGeom & GR);
			long GetId(const Triangle & t) const;
			long GetId(const Triangle * t) const;
			long GetId(const BamgVertex & t) const;
			long GetId(const BamgVertex * t) const;
			long GetId(const Edge & t) const;
			long GetId(const Edge * t) const;
			BamgVertex* NearestVertex(int i,int j) ;
			Triangle* TriangleFindFromCoord(const I2 & ,long long [3],Triangle *tstart=0);
			void ReadMesh(int* index,double* x,double* y,int nods,int nels,BamgOpts* bamgopts);
			void ReadMesh(BamgMesh* bamgmesh, BamgOpts* bamgopts);
			void WriteMesh(BamgMesh* bamgmesh,BamgOpts* bamgopts);
			void ReadMetric(const BamgOpts* bamgopts);
			void WriteMetric(BamgOpts* bamgopts);
			void WriteIndex(int** pindex,int* pnels);
			void AddMetric(BamgOpts* bamgopts);
			void BuildMetric0(BamgOpts* bamgopts);
			void BuildMetric1(BamgOpts* bamgopts);
			void BuildGeometryFromMesh(BamgOpts* bamgopts=NULL);
			long RandomNumber(long max);
			void ReconstructExistingMesh(BamgOpts* bamgopts);

			//Inline methods
			inline  void CreateSingleVertexToTriangleConnectivity(){
				for (int i=0;i<nbv;i++) vertices[i].IndexInTriangle=0, vertices[i].t=NULL;
				for (int i=0;i<nbt;i++) triangles[i].SetSingleVertexToTriangleConnectivity();
			}
			inline  void  UnMarkUnSwapTriangle(){
				for (int i=0;i<nbt;i++)
				 for(int j=0;j<3;j++)
				  triangles[i].SetUnMarkUnSwap(j);
			  }
			inline  void  SetVertexFieldOn(){
				for (int i=0;i<nbv;i++)                    vertices[i].GeomEdgeHook=NULL;
				for (int j=0;j<NbVerticesOnGeomVertex;j++) VerticesOnGeomVertex[j].SetOn();
				for (int k=0;k<NbVerticesOnGeomEdge;k++ )  VerticesOnGeomEdge[k].SetOn();
			}	       
			inline  void   SetVertexFieldOnBTh(){
				for (int i=0;i<nbv;i++)                 vertices[i].GeomEdgeHook=NULL;
				for (int j=0;j<NbVertexOnBThVertex;j++) VertexOnBThVertex[j].SetOnBTh();
				for (int k=0;k<NbVertexOnBThEdge;k++ )  VertexOnBThEdge[k].SetOnBTh();
			}

		private:
			void TriangulateFromGeom1(BamgOpts* bamgopts,int KeepVertices=1);// the real constructor mesh adaption
			void TriangulateFromGeom0(BamgOpts* bamgopts);// the real constructor mesh generator
			void Triangulate(double* x,double* y,int nods,BamgOpts* bamgopts);
			void Init(long);
			int ForceEdge(BamgVertex &a, BamgVertex & b,AdjacentTriangle & taret) ;
			int SwapForForcingEdge(BamgVertex   *  & pva ,BamgVertex  * &   pvb ,
						AdjacentTriangle & tt1,long long & dets1,
						long long & detsa,long long & detsb, int & nbswap);
	};

	/*Intermediary*/
	AdjacentTriangle CloseBoundaryEdge(I2 ,Triangle *, double &,double &) ;
	void  swap(Triangle *t1,short a1,
				Triangle *t2,short a2,
				BamgVertex *s1,BamgVertex *s2,long long det1,long long det2);

	inline AdjacentTriangle Previous(const AdjacentTriangle & ta){
		return AdjacentTriangle(ta.t,PreviousEdge[ta.a]);
	}
	inline AdjacentTriangle Next(const AdjacentTriangle & ta){
		return AdjacentTriangle(ta.t,NextEdge[ta.a]);
	}
	inline  AdjacentTriangle Adj(const AdjacentTriangle & a){
		return  a.Adj();
	}
	inline void Adj(GeomEdge * & on,int &i){
		int j=i;i=on->AdjVertexIndex[i];on=on->Adj[j];
	}
}
#endif
