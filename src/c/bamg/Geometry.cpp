#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "./bamgobjects.h"
#include "../shared/shared.h"

namespace bamg {

	/*Constructors/Destructors*/
	Geometry::Geometry(){/*{{{*/
		Init();
	}
	/*}}}*/
	Geometry::Geometry(BamgGeom* bamggeom, BamgOpts* bamgopts){/*{{{*/
		Init();
		ReadGeometry(bamggeom,bamgopts);
		PostRead();
	}
	/*}}}*/
	Geometry::Geometry(const Geometry & Gh) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/Geometry)*/

		long i;
		*this = Gh;
		NbRef =0;
		quadtree=0;
		vertices = nbv ? new GeomVertex[nbv] : NULL;
		edges = nbe ? new GeomEdge[nbe]:NULL;
		curves= nbcurves ? new Curve[nbcurves]:NULL;
		subdomains = nbsubdomains ? new GeomSubDomain[nbsubdomains]:NULL;
		for (i=0;i<nbe;i++)
		 edges[i].Set(Gh.edges[i],Gh,*this);
		for (i=0;i<nbcurves;i++)
		 curves[i].Set(Gh.curves[i],Gh,*this);
		for (i=0;i<nbsubdomains;i++)
		 subdomains[i].Set(Gh.subdomains[i],Gh,*this);
	}
	/*}}}*/
	Geometry::~Geometry() {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/~Geometry)*/
		if(NbRef>0){   _printf_("Trying to delete geometry and NbRef>0, probably due to an error. NbRef:" << NbRef); return;}
		if(vertices)   delete [] vertices;  vertices=0;
		if(edges)      delete [] edges;     edges=0;
		if(quadtree)   delete    quadtree;  quadtree=0;
		if(curves)     delete [] curves;    curves=0;nbcurves=0;
		if(subdomains) delete [] subdomains;subdomains=0;
		Init();
	}
	/*}}}*/

	/*IO*/
	void Geometry::ReadGeometry(BamgGeom* bamggeom,BamgOpts* bamgopts){/*{{{*/

		int verbose;
		nbcurves=0;

		double Hmin = HUGE_VAL;// the infinie value
		int i,j,n,i0,i1,i2,i3;

		/*initialize some variables*/
		verbose= bamgopts->verbose;
		nbv    = bamggeom->VerticesSize[0];
		nbe    = bamggeom->EdgesSize[0];

		//some checks
		if (bamggeom->Vertices==NULL) _error_("the domain provided does not contain any vertex");
		if (bamggeom->Edges==NULL) _error_("the domain provided does not contain any edge");

		//Vertices
		if (bamggeom->Vertices){
			if(verbose>5) _printf_("      processing Vertices\n");
			if (bamggeom->VerticesSize[1]!=3) _error_("Vertices should have 3 columns");
			vertices = new GeomVertex[nbv];
			for (i=0;i<nbv;i++) {
				vertices[i].r.x=(double)bamggeom->Vertices[i*3+0];
				vertices[i].r.y=(double)bamggeom->Vertices[i*3+1];
				vertices[i].ReferenceNumber=(long)bamggeom->Vertices[i*3+2];
				vertices[i].color =0;
				vertices[i].type=0;
			}
			/*find domain extrema (pmin,pmax) that will define the square box used for by the quadtree*/
			pmin =  vertices[0].r;
			pmax =  vertices[0].r;
			for (i=0;i<nbv;i++) {
				pmin.x = Min(pmin.x,vertices[i].r.x);
				pmin.y = Min(pmin.y,vertices[i].r.y);
				pmax.x = Max(pmax.x,vertices[i].r.x);
				pmax.y = Max(pmax.y,vertices[i].r.y);
			}
			/*Offset pmin and pmax to avoid round-off errors*/
			R2 offset = (pmax-pmin)*0.05;
			pmin -= offset;
			pmax += offset;
			/*coefIcoor is the coefficient used for integer coordinates:
			 *                       (x-pmin.x)
			 * Icoor x = (2^30 -1) ------------
			 *                          D
			 * where D is the longest side of the domain (direction x or y)
			 * so that (x-pmin.x)/D is in ]0 1[
			 *
			 * coefIcoor = (2^30 -1)/D
			 */
			int MaxICoord = 1073741823;
			coefIcoor=(MaxICoord)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
			if(coefIcoor<=0) _error_("coefIcoor should be positive");
		}
		else{
			_error_("No BamgVertex provided");
		}

		//Edges
		if (bamggeom->Edges){
			R2      zerovector(0,0);
			double* verticeslength=NULL;

			if(verbose>5) _printf_("      processing Edges\n");
			if (bamggeom->EdgesSize[1]!=3) _error_("Edges should have 3 columns");
			edges = new GeomEdge[nbe];

			//initialize verticeslength (sum of the lengths of the edges holding vertex)
			verticeslength = new double[nbv];
			for(i=0;i<nbv;i++) verticeslength[i]=0;

			/*Loop over the edges*/
			for (i=0;i<nbe;i++){

				i1=(int)bamggeom->Edges[i*3+0]-1; //-1 for C indexing
				i2=(int)bamggeom->Edges[i*3+1]-1; //-1 for C indexing
				edges[i].v[0]= vertices + i1;     //pointer toward vertex i1 (=&vertices[i1])
				edges[i].v[1]= vertices + i2;     //pointer toward vertex i2
				edges[i].ReferenceNumber=(long)bamggeom->Edges[i*3+2];

				//get length of edge
				R2     x12=vertices[i2].r-vertices[i1].r;
				double l12=sqrt((x12,x12));
				Hmin=Min(Hmin,l12);

				//initialize other fields
				edges[i].tg[0]=zerovector;
				edges[i].tg[1]=zerovector;
				edges[i].AdjVertexIndex[0] = edges[i].AdjVertexIndex[1] = -1;
				edges[i].Adj[0] = edges[i].Adj[1] = NULL;
				edges[i].type = 0;

				//Cracked?
				if (edges[i].ReferenceNumber!=1) edges[i].SetCracked();

				//prepare metric
				vertices[i1].color++;
				vertices[i2].color++;
				verticeslength[i1] += l12;
				verticeslength[i2] += l12;
			}

			// definition the default of the given mesh size
			for (i=0;i<nbv;i++) {
				if (vertices[i].color > 0)
				 vertices[i].m=Metric(verticeslength[i] /(double) vertices[i].color);
				else
				 vertices[i].m=Metric(Hmin);
			}
			delete [] verticeslength;

		}
		else{
			_error_("No edges provided");
		}

		//hVertices
		if(bamgopts->hVertices && bamgopts->hVerticesLength==nbv){
			if(verbose>5) _printf_("      processing hVertices\n");
			for (i=0;i< nbv;i++){
				if (!xIsNan<IssmPDouble>(bamgopts->hVertices[i])){
					vertices[i].m=Metric((double)bamgopts->hVertices[i]);
				}
			}
		}

		//MetricVertices
		if(bamgopts->metric && bamgopts->metric[0]==nbv){
			if(verbose>5) _printf_("      processing MetricVertices\n");
			for (i=0;i< nbv;i++) {
				vertices[i].m = Metric((double)bamgopts->metric[i*3+0],(double)bamgopts->metric[i*3+1],(double)bamgopts->metric[i*3+2]);
			}
		}

		//TangentAtEdges
		if (bamggeom->TangentAtEdges){
			if(verbose>5) _printf_("      processing TangentAtEdges");
			if (bamggeom->TangentAtEdgesSize[1]!=4) _error_("TangentAtEdges should have 4 columns");
			int n,i,j,k;
			R2 tg;

			n=bamggeom->TangentAtEdgesSize[0];
			for (k=0;k<n;k++) {
				i=(int)bamggeom->TangentAtEdges[k*4+0]-1; //for C indexing
				j=(int)bamggeom->TangentAtEdges[k*4+1]-1; //for C indexing
				tg.x=bamggeom->TangentAtEdges[k*4+2];
				tg.y=bamggeom->TangentAtEdges[k*4+3];
				if (i<0 || i>=nbe) _error_("TangentAtEdges first index exceeds matrix dimension");
				if (j!=0 && j!=1)  _error_("TangentAtEdges second index should be 1 or 2 only");
				edges[i].tg[j] = tg;
			}
		}

		//Corners
		if(bamggeom->Corners){
			if(verbose>5) _printf_("      processing Corners");
			if (bamggeom->CornersSize[1]!=1) _error_("Corners should have 1 column");
			n=bamggeom->CornersSize[0];
			for (i=0;i<n;i++) {
				j=(int)bamggeom->Corners[i]-1; //for C indexing
				if (j>nbv-1 || j<0) _error_("Bad corner definition: should in [0 " << nbv << "]");
				/*Required => at the same time SetRequired and SetCorner*/
				vertices[j].SetCorner();
				vertices[j].SetRequired();
			}
		}

		//RequiredVertices
		if(bamggeom->RequiredVertices){
			if(verbose>5) _printf_("      processing RequiredVertices\n");
			if (bamggeom->RequiredVerticesSize[1]!=1) _error_("RequiredVertices should have 1 column");
			n=bamggeom->RequiredVerticesSize[0];
			for (i=0;i<n;i++) {
				j=(int)bamggeom->RequiredVertices[i]-1; //for C indexing
				if (j>nbv-1 || j<0) _error_("Bad RequiredVerticess  definition: should in [0 " << nbv << "]");
				vertices[j].SetRequired();
			}
		}

		//RequiredEdges
		if(bamggeom->RequiredEdges){
			if(verbose>5) _printf_("      processing RequiredEdges\n");
			if (bamggeom->RequiredEdgesSize[1]!=1) _error_("RequiredEdges should have 1 column");
			n=bamggeom->RequiredEdgesSize[0];
			for (i=0;i<n;i++) {
				j=(int)bamggeom->RequiredEdges[i]-1; //for C indexing
				if (j>nbe-1 || j<0) _error_("Bad RequiredEdges definition: should in [0 " << nbe << "]");
				edges[j].SetRequired();
			}
		}

		//SubDomain
		if(bamggeom->SubDomains){
			if(verbose>5) _printf_("      processing SubDomains\n");
			if (bamggeom->SubDomainsSize[1]!=4) _error_("SubDomains should have 4 columns");
			nbsubdomains=bamggeom->SubDomainsSize[0];
			subdomains = new GeomSubDomain[nbsubdomains];
			for (i=0;i<nbsubdomains;i++){
				i0=(int)bamggeom->SubDomains[i*4+0];
				i1=(int)bamggeom->SubDomains[i*4+1];
				i2=(int)bamggeom->SubDomains[i*4+2];
				i3=(int)bamggeom->SubDomains[i*4+3];
				if (i0!=2) _error_("Bad Subdomain definition: first number should be 2 (for Edges)");
				if (i1>nbe || i1<=0) _error_("Bad Subdomain definition: second number should in [1 " << nbe << "] (edge number)");
				subdomains[i].edge=edges + (i1-1);
				subdomains[i].direction = (int) i2;
				subdomains[i].ReferenceNumber = i3;
			}
		}
	}
	/*}}}*/
	void Geometry::WriteGeometry(BamgGeom* bamggeom, BamgOpts* bamgopts){/*{{{*/

		int verbose;
		int nbreq=0;
		int nbreqv=0;
		int nbtan=0;
		int i,count;

		/*Get options*/
		verbose=bamgopts->verbose;

		/*Vertices*/
		if(verbose>5) _printf_("      writing Vertices\n");
		bamggeom->VerticesSize[0]=nbv;
		bamggeom->VerticesSize[1]=3;
		if (nbv){
			bamggeom->Vertices=xNew<double>(3*nbv);
			for (i=0;i<nbv;i++){
				bamggeom->Vertices[i*3+0]=vertices[i].r.x;
				bamggeom->Vertices[i*3+1]=vertices[i].r.y;
				bamggeom->Vertices[i*3+2]=vertices[i].GetReferenceNumber();

				//update counters
				if (vertices[i].Required()) nbreqv++;
			}
		}

		/*Edges*/
		if(verbose>5) _printf_("      writing Edges\n");
		bamggeom->EdgesSize[0]=nbe;
		bamggeom->EdgesSize[1]=3;
		if (nbe){
			bamggeom->Edges=xNew<double>(3*nbe);
			for (i=0;i<nbe;i++){
				bamggeom->Edges[i*3+0]=GetId(edges[i][0])+1; //back to Matlab indexing
				bamggeom->Edges[i*3+1]=GetId(edges[i][1])+1; //back to Matlab indexing
				bamggeom->Edges[i*3+2]=(double)edges[i].ReferenceNumber;

				//update counters
				if (edges[i].Required()) nbreq++;
				if (edges[i].TgA() && edges[i][0].Corner()) nbtan++;
				if (edges[i].TgB() && edges[i][1].Corner()) nbtan++;
			}
		}

		/*RequiredEdges*/
		if(verbose>5) _printf_("      writing " << nbreq << " RequiredEdges\n");
		bamggeom->RequiredEdgesSize[0]=nbreq;
		bamggeom->RequiredEdgesSize[1]=1;
		if (nbreq){
			bamggeom->RequiredEdges=xNew<double>(1*nbreq);
			count=0;
			for (i=0;i<nbe;i++){
				if (edges[i].Required()){
					bamggeom->RequiredEdges[count]=i+1; //back to Matlab indexing
					count=count+1;
				}
			}
		}

		//No corners

		/*RequiredVertices*/
		if(verbose>5) _printf_("      writing " << nbreqv << " RequiredVertices\n");
		bamggeom->RequiredVerticesSize[0]=nbreqv;
		bamggeom->RequiredVerticesSize[1]=1;
		if (nbreqv){
			bamggeom->RequiredVertices=xNew<double>(1*nbreqv);
			count=0;
			for (i=0;i<nbv;i++){
				if (vertices[i].Required()){
					bamggeom->RequiredVertices[count]=i+1; //back to Matlab indexing
					count=count+1;
				}
			}
		}

		/*SubDomains*/
		if(verbose>5) _printf_("      writing SubDomains\n");
		bamggeom->SubDomainsSize[0]=nbsubdomains;
		bamggeom->SubDomainsSize[1]=4;
		if (nbsubdomains){
			bamggeom->SubDomains=xNew<double>(4*nbsubdomains);
			for (i=0;i<nbsubdomains;i++){
				bamggeom->SubDomains[4*i+0]=2;
				bamggeom->SubDomains[4*i+1]=GetId(subdomains[i].edge)+1; //back to Matlab indexing
				bamggeom->SubDomains[4*i+2]=subdomains[i].direction;
				bamggeom->SubDomains[4*i+3]=subdomains[i].ReferenceNumber;
			}
		}

		/*TangentAtEdges*/
		if(verbose>5) _printf_("      writing TangentAtEdges\n");
		bamggeom->TangentAtEdgesSize[0]=nbtan;
		bamggeom->TangentAtEdgesSize[1]=4;
		if (nbtan){
			bamggeom->TangentAtEdges=xNew<double>(4*nbtan);
			for (i=0;i<nbe;i++){
				if (edges[i].TgA() && edges[i][0].Corner()){
					bamggeom->TangentAtEdges[4*i+0]=i+1; //back to Matlab indexing
					bamggeom->TangentAtEdges[4*i+1]=1;
					bamggeom->TangentAtEdges[4*i+2]=edges[i].tg[0].x;
					bamggeom->TangentAtEdges[4*i+3]=edges[i].tg[0].y;
				}
				if (edges[i].TgB() && edges[i][1].Corner()){
					bamggeom->TangentAtEdges[4*i+0]=i+1; //back to Matlab indexing
					bamggeom->TangentAtEdges[4*i+1]=2;
					bamggeom->TangentAtEdges[4*i+2]=edges[i].tg[1].x;
					bamggeom->TangentAtEdges[4*i+3]=edges[i].tg[1].y;
				}
			}
		}
	}
	/*}}}*/

	/*Methods*/
	void Geometry::Echo(void){/*{{{*/

		_printf_("Geometry:\n");
		_printf_("   nbv  (number of vertices) : " << nbv << "\n");
		_printf_("   nbe  (number of edges)    : " << nbe << "\n");
		_printf_("   nbsubdomains: " << nbsubdomains << "\n");
		_printf_("   nbcurves: " << nbcurves << "\n");
		_printf_("   vertices: " << vertices << "\n");
		_printf_("   edges: " << edges << "\n");
		_printf_("   quadtree: " << quadtree << "\n");
		_printf_("   subdomains: " << subdomains << "\n");
		_printf_("   curves: " << curves << "\n");
		_printf_("   pmin (x,y): (" << pmin.x << " " << pmin.y << ")\n");
		_printf_("   pmax (x,y): (" << pmax.x << " " << pmax.y << ")\n");
		_printf_("   coefIcoor: " << coefIcoor << "\n");

		return;
	}
	/*}}}*/
	void Geometry::Init(void){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/EmptyGeometry)*/

		NbRef=0;
		nbv=0;
		nbe=0;
		quadtree=NULL;
		curves=NULL;
		edges=NULL;
		vertices=NULL;
		nbsubdomains=0;
		nbcurves=0;
		subdomains=NULL;
	}
	/*}}}*/
	double Geometry::MinimalHmin() {/*{{{*/
		/* coeffIcoor = (2^30-1)/D
		 * We cannot go beyond hmin = D/2^30 because of the quadtree
		 * hmin is therefore approximately 2/coeffIcoor */
		return 2.0/coefIcoor;
	}/*}}}*/
	double Geometry::MaximalHmax() {/*{{{*/
		return Max(pmax.x-pmin.x,pmax.y-pmin.y);
	}/*}}}*/
	long Geometry::GetId(const GeomVertex & t) const  {/*{{{*/
		return &t - vertices;
	}/*}}}*/
	long Geometry::GetId(const GeomVertex * t) const  {/*{{{*/
		return t - vertices;
	}/*}}}*/
	long Geometry::GetId(const GeomEdge & t) const  {/*{{{*/
		return &t - edges;
	}/*}}}*/
	long Geometry::GetId(const GeomEdge * t) const  {/*{{{*/
		return t - edges;
	}/*}}}*/
	long Geometry::GetId(const Curve * c) const  {/*{{{*/
		return c - curves;
	}/*}}}*/
	void Geometry::PostRead(bool checkcurve){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/AfterRead)*/

		long          i          ,j,k;
		long         *head_v   = new long[nbv];
		long         *next_p   = new long[2 *nbe];
		float        *eangle   = new float[nbe];
		double        eps      = 1e-20;
		BamgQuadtree  quadtree;                            // build quadtree to find duplicates

		k=0;

		//build quadtree for this geometry
		for (i=0;i<nbv;i++){

			/*build integer coordinates (non unique)
			  these coordinates are used by the quadtree to group
			  the vertices by groups of 5:
			  All the coordinates are transformed to ]0,1[^2
			  then, the integer coordinates are computed using
			  the transformation ]0,1[^2 -> [0 2^30-1[^2 for a quadtree of depth 30*/
			vertices[i].i=R2ToI2(vertices[i].r);

			/*find nearest vertex already present in the quadtree (NULL if empty)*/
			BamgVertex* v=quadtree.NearestVertex(vertices[i].i.x,vertices[i].i.y);

			/*if there is a vertex found that is to close to vertices[i] -> error*/
			if( v && Norme1(v->r - vertices[i].r) < eps ){
				_printf_("reference numbers: " << v->ReferenceNumber << " " << vertices[i].ReferenceNumber << "\n");
				_printf_("Id: " << i+1 << "\n");
				_printf_("Coords: ["<<v->r.x<<" "<<v->r.y<<"] ["<<vertices[i].r.x<<" "<<vertices[i].r.y<<"]\n");

				delete [] next_p;
				delete [] head_v;
				delete [] eangle;
				_error_("two points of the geometry are very close to each other (see reference numbers above)");
			}

			/*Add vertices[i] to the quadtree*/
			quadtree.Add(vertices[i]);
		}

		/* Here we use a powerful chaining algorithm
		 *
		 * 1. What is a chaining algorithm?
		 *
		 * If F is a function that goes from i in [0 n] to j in [0 m]
		 * and we want to compute the reciprocal function F-1 of F
		 * (what are the antecedents of a given j in Im(F) )
		 * We use 2 lists:
		 *    head_F[j] that holds the head of lists
		 *    next_F[i] that holds the list of elements that have the same image
		 *
		 * Example:
		 *    i1, i2, ..., ip in [0,n] are all antecedents of a given j in [0 m]
		 *    head_F[j] = ip
		 *    next_F[ip]= ip-1
		 *    ....
		 *    next_F[i2]= i1
		 *    next_F[i1]= -1  //end of the list
		 *
		 * Algorithm:
		 *    for(j=0;j<m;j++)  head_F[j] = -1 //initialization
		 *    for(i=0;i<n;i++){
		 *       j=F[i];
		 *       next_F[i]= head_F[j];
		 *       head_F[j]=i;
		 *    }
		 *
		 *    Then, we can go through all the elements that have for image j:
		 *    for(i=head_F[j]; i!=-1; i=next_F[i])
		 *    initialization of i by i=head_F[j]
		 *    stop the loop when i=-1 (end of the chain)
		 *    iterate using i=next_F[i] (next element that have for image j)
		 *
		 * 2. How to use this algorithm here?
		 *
		 * Here F is a function that associates two vertices v0 and v1 for a given edge E
		 * We want to build the reciprocal function: what are the edges that contains
		 * a vertex v?
		 * To do so, we use the same chaining algorithm but there is a difficulty
		 * coming from the fact that for F we have a couple of vertices and not one
		 * vertex.
		 * To overcome this difficulty, we use a global indexing exactly like in
		 * C/C++ so that
		 * a member of a 2-column-table can be described by one index p=i*2+j
		 * i=(int)p/2 line number of p
		 * j=p%2      column number of p
		 *
		 * Algorithm:
		 *    for(i=0;i<nbv;i++)  head_v[i] = -1 //initialization
		 *    for(i=0;i<nbe;i++){
		 *       for(j=0;j<2;j++){
		 *          p=2*i+j;
		 *          v=edges(i,j);
		 *          next_p[p]= head_v[v];
		 *          head_v[v]=p;
		 *       }
		 *    }
		 */

		//initialize head_v as -1
		for (i=0;i<nbv;i++) head_v[i]=-1;
		k=0;
		for (i=0;i<nbe;i++) {
			//compute vector of edge i that goes from vertex 0 to vertex 1
			R2 v10=edges[i].v[1]->r - edges[i].v[0]->r;
			double lv10=Norme2(v10);
			//check that its length is not 0
			if(lv10==0){
				delete [] next_p;
				delete [] head_v;
				delete [] eangle;
                _printf_("Length of edge " << i << " is 0\n");
				_error_("Length of edge " << i << " is 0");
			}
			//compute angle in [-Pi Pi]
			eangle[i] = atan2(v10.y,v10.x);
			//build chains head_v and next_p
			for (j=0;j<2;j++){
				long v=GetId(edges[i].v[j]);
				next_p[k]=head_v[v];
				head_v[v]=k++; //post increment: head_v[v]=k; and then k=k+1;
			}
		}

		//sort head_v by order of increasing edges angle
		for (i=0;i<nbv;i++) {
			int exch=1,ord=0;

			//exchange vertices position in head_v and next_p till tey are sorted
			while (exch){
				long *p=head_v+i;
				long *po=p;
				long  n=*p;
				float angleold=-1000 ; // angle = - infinity
				ord=0; exch=0;

				// loop over the edges that contain the vertex i (till -1)
				while (n >=0){
					ord++;
					long  i1=n/2;       // i1 = floor (n/2)      -> Edge number
					long  j1=n%2;       // j1 = 1 if n is odd    -> Vertex index for this edge (0 or 1)
					long* pn=next_p+n;

					//Next vertex index
					n=*pn;

					//compute angle between horizontal axis and v0->v1
					float angle = j1 ? OppositeAngle(eangle[i1]):  eangle[i1];

					//exchange if the current edge angle is smaller than the previous one
					if (angleold > angle){
						exch=1;
						*pn=*po;  // next_p[n] = n + 1
						*po=*p;   //
						*p=n;     // next_p[n+1] = n
						po=pn;    // po now point toward pn (invert next and current)
					}

					//else, continue going to the next edge position
					else{                        //  to have : po -> p -> pn
						angleold=angle; // update maximum angle
						po=p;           // po now point toward p  (current position)
						p=pn;           // p  now point toward pn (next position)
					}
				}
			}

			/*Do we want to check for curve? Default is no, but if we are reconstructing, it's better to turn it on with a small threshold*/
			if(checkcurve){
				/*angular test on current vertex to guess whether it is a corner (ord = number of edges holding i) */
				if(ord==2){
					IssmDouble MaxCornerAngle = 1*Pi/180; /*default is 1 degree*/
					long  n1 = head_v[i];
					long  n2 = next_p[n1];
					long  i1 = n1/2, i2 = n2/2; // edge number
					long  j1 = n1%2, j2 = n2%2; // vertex in the edge
					float angle1=  j1 ? OppositeAngle(eangle[i1]) : eangle[i1];
					float angle2= !j2 ? OppositeAngle(eangle[i2]) : eangle[i2];
					float da12 = Abs(angle2-angle1);
					if (( da12 >= MaxCornerAngle ) && (da12 <= 2*Pi -MaxCornerAngle)) {
						vertices[i].SetCorner() ;
					}
					/* if the edge type/referencenumber a changing then is SetRequired();*/
					if(edges[i1].type != edges[i2].type || edges[i1].Required()){
						vertices[i].SetRequired();
					}
					if(edges[i1].ReferenceNumber != edges[i2].ReferenceNumber) {
						vertices[i].SetRequired();
					}
				}
				if(ord!=2) vertices[i].SetCorner();
			}
			else{
				/*all vertices provided in geometry are corners (ord = number of edges holding i)*/
				vertices[i].SetCorner() ;
				if(ord==2){
					long  n1 = head_v[i];
					long  n2 = next_p[n1];
					long  i1 = n1/2, i2 = n2/2; // edge number
					long  j1 = n1%2, j2 = n2%2; // vertex in the edge
					/* if the edge type/referencenumber a changing then is SetRequired();*/
					if (edges[i1].type != edges[i2].type || edges[i1].Required()){
						vertices[i].SetRequired();
					}
					if (edges[i1].ReferenceNumber != edges[i2].ReferenceNumber) {
						vertices[i].SetRequired();
					}
				}
			}

			/*close the list around the vertex to have a circular loop*/
			long no=-1, ne = head_v[i];
			while (ne >=0) ne = next_p[no=ne];
			if(no>=0) next_p[no] = head_v[i];
		}

		/*Check that the list makes sense (we have all the time the same vertex)
		 * and build adjacent edges*/
		k =0;
		for (i=0;i<nbe;i++){
			for (j=0;j<2;j++){

				long n1 = next_p[k++];
				long i1 = n1/2 ,j1=n1%2;

				if( edges[i1].v[j1] != edges[i].v[j]) _error_("Problem while processing edges: check the edge list");

				edges[i1].Adj[j1] = edges + i;
				edges[i1].AdjVertexIndex[j1] = j;
			}
		}

		/* generation of  all the tangents*/
		for (i=0;i<nbe;i++) {
			R2    AB =edges[i].v[1]->r -edges[i].v[0]->r;// AB = vertex0 -> vertex1
			double lAB=Norme2(AB);                       // Get length of AB
			double ltg2[2]={0.0};                        // initialize tangent

			//loop over the 2 vertices of the edge
			for (j=0;j<2;j++) {
				R2     tg =edges[i].tg[j];
				double ltg=Norme2(tg);

				//by default, tangent=[0 0]
				if(ltg==0){
					//if the current vertex of the edge is not a corner
					if(!edges[i].v[j]->Corner()){
						/*The tangent is set as the vector between the
						 * previous and next vertices connected to current vertex
						 * normed by the edge length*/
						tg = edges[i].v[1-j]->r - edges[i].Adj[j]->v[1-edges[i].AdjVertexIndex[j]]->r;
						ltg= Norme2(tg);
						tg = tg *(lAB/ltg);
						ltg= lAB;
					}
					//else:  a Corner no tangent => nothing to do
				}
				else{
					//tangent has already been computed
					tg = tg*(lAB/ltg),ltg=lAB;
				}

				ltg2[j] = ltg;

				if ((tg,AB)<0) tg = -tg;

				edges[i].tg[j]=tg;
			}
			if (ltg2[0]!=0) edges[i].SetTgA();
			if (ltg2[1]!=0) edges[i].SetTgB();
		}

		/* generation of  all curves (from corner to corner)*/
		/*We proceed in 2 steps: first allocate, second build*/
		for (int step=0;step<2;step++){

			//unmark all edges
			for (i=0;i<nbe;i++) edges[i].SetUnMark();
			long  nb_marked_edges=0;

			//initialize number of curves
			nbcurves = 0;

			for (int level=0;level<2 && nb_marked_edges!=nbe;level++){
				for (i=0;i<nbe;i++){

					GeomEdge & ei=edges[i];
					for(j=0;j<2;j++){
						/*If current edge ei is unmarked and (level=1 or vertex i is required (corner)):
						 * we do have the first edge of a new curve*/
						if (!ei.Mark() && (level || ei[j].Required())) {
							int k0=j,k1;
							GeomEdge   *e=&ei;
							GeomVertex *a=(*e)(k0); // begin
							if(curves){
								curves[nbcurves].FirstEdge=e;
								curves[nbcurves].FirstVertexIndex=k0;
							}
							int nee=0;
							for(;;){
								nee++;
								k1 = 1-k0; // next vertex of the edge
								e->SetMark();
								nb_marked_edges++;
								e->CurveNumber=nbcurves;
								GeomVertex *b=(*e)(k1);

								//break if we have reached the other end of the curve
								if (a==b || b->Required()){
									if(curves){
										curves[nbcurves].LastEdge=e;
										curves[nbcurves].LastVertexIndex=k1;
									}
									break;
								}
								//else: go to next edge (adjacent)
								else{
									k0 = e->AdjVertexIndex[k1];//  vertex in next edge
									e  = e->Adj[k1]; // next edge
								}
							}
							nbcurves++;
							if(level) a->SetRequired();
						}
					}
				}
			}
			_assert_(nb_marked_edges && nbe);
			//allocate if first step
			if(step==0) curves=new Curve[nbcurves];
		}

		/*clean up*/
		delete [] next_p;
		delete [] head_v;
		delete [] eangle;

	}
	/*}}}*/
	GeomEdge* Geometry::ProjectOnCurve(const Edge &e,double s,BamgVertex &V,VertexOnGeom &GV) const {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/ProjectOnCurve)*/
		/*Add a vertex on an existing geometrical edge according to the metrics of the two vertices constituting the edge*/
		/*FIXME: should go away*/

		double save_s=s;
		int NbTry=0;

retry:

		s=save_s;
		GeomEdge* on=e.GeomEdgeHook;
		if (!on){
			_error_("ProjectOnCurve error message: edge provided should be on geometry");
		}
		if (!e[0].GeomEdgeHook ||  !e[1].GeomEdgeHook){
			_error_("ProjectOnCurve error message: at least one of the vertex of the edge provided is not on geometry");
		}

		//Get the two vertices of the edge
		const BamgVertex &v0=e[0];
		const BamgVertex &v1=e[1];

		//Get position of V0, V1 and vector v0->v1
		R2 V0=v0,V1=v1,V01=V1-V0;

		//Get geometrical vertices corresponding to v0 and v1
		VertexOnGeom  vg0=*v0.GeomEdgeHook,  vg1=*v1.GeomEdgeHook;

		//build two pointers towrad current geometrical edge
		GeomEdge *eg0=on, *eg1=on;

		//Get edge direction and swap v0 and v1 if necessary
		R2 Ag=(R2)(*on)[0],Bg=(R2)(*on)[1],AB=Bg-Ag;
		int OppositeSens = (V01,AB)<0;
		int direction0=0,direction1=1;
		if (OppositeSens) s=1-s,Exchange(vg0,vg1),Exchange(V0,V1);

		//Compute metric of new vertex as a linear interpolation of the two others
		V.m=Metric(1.0-s,v0.m,s,v1.m);

		const int mxe=100;
		GeomEdge* ge[mxe+1];
		int     directionge[mxe+1];
		double  lge[mxe+1];
		int bge=mxe/2,tge=bge;
		ge[bge] = e.GeomEdgeHook;
		directionge[bge]=1;

		while(eg0!=(GeomEdge*)vg0 && (*eg0)(direction0)!=(GeomVertex*)vg0){
			if (bge<=0) {
				if(NbTry) {
					_printf_("Fatal Error: on the class Mesh before call Geometry::ProjectOnCurve\n");
					_printf_("That bug might come from:\n");
					_printf_(" 1)  a mesh edge  containing more than " << mxe/2 << " geometrical edges\n");
					_printf_(" 2)  code bug : be sure that we call   Mesh::SetVertexFieldOn() before\n");
					_printf_("To solve the problem do a coarsening of the geometrical mesh or change the constant value of mxe (dangerous)\n");
					_error_("see above");
				}
				NbTry++;
				goto retry;
			}
			GeomEdge* tmpge = eg0;
			ge[--bge] =eg0 = eg0->Adj[direction0];
			_assert_(bge>=0 && bge<=mxe);
			direction0 = 1-( directionge[bge] = tmpge->AdjVertexIndex[direction0]);
		}
		while (eg1 != (GeomEdge*) vg1  &&  (*eg1)(direction1) != (GeomVertex*) vg1) {
			if(tge>=mxe ) {
				_printf_("WARNING: on the class Mesh before call Geometry::ProjectOnCurve is having issues (isn't it Eric?)\n");
				NbTry++;
				if (NbTry<2) goto retry;
				_printf_("Fatal Error: on the class Mesh before call Geometry::ProjectOnCurve\n");
				_printf_("That bug might come from:\n");
				_printf_(" 1)  a mesh edge  contening more than " << mxe/2 << " geometrical edges\n");
				_printf_(" 2)  code bug : be sure that we call   Mesh::SetVertexFieldOn() before\n");
				_printf_("To solve the problem do a coarsening of the geometrical mesh or change the constant value of mxe (dangerous)\n");
				_error_("see above");
			}
			GeomEdge* tmpge = eg1;
			ge[++tge] =eg1 = eg1->Adj[direction1];
			directionge[tge]= direction1 = 1-tmpge->AdjVertexIndex[direction1];
			_assert_(tge>=0 && tge<=mxe);
		}

		if ((*eg0)(direction0)==(GeomVertex*)vg0)
		 vg0=VertexOnGeom(*(BamgVertex*) vg0,*eg0,direction0); //vg0 = absisce

		if ((*eg1)(direction1)==(GeomVertex*)vg1)
		 vg1=VertexOnGeom(*(BamgVertex*) vg1,*eg1,direction1);

		double sg;
		if (eg0 == eg1) {
			double s0=vg0,s1=vg1;
			sg =  s0*(1.0-s) +  s*s1;
			on=eg0;
		}
		else {
			R2 AA=V0,BB;
			double s0,s1;
			int i=bge;
			double ll=0;
			for(i=bge;i<tge;i++){
				_assert_(i>=0 && i<=mxe);
				BB =  (*ge[i])[directionge[i]];
				lge[i]=ll += Norme2(AA-BB);
				AA=BB ;}
				lge[tge]=ll+=Norme2(AA-V1);
				// search the geometrical edge
				_assert_(s<=1.0);
				double ls= s*ll;
				on =0;
				s0 = vg0;
				s1= directionge[bge];
				double l0=0,l1;
				i=bge;
				while (  (l1=lge[i]) < ls ) {
					_assert_(i>=0 && i<=mxe);
					i++,s0=1-(s1=directionge[i]),l0=l1;
				}
				on=ge[i];
				if (i==tge)
				 s1=vg1;

				s  =(ls-l0)/(l1-l0);
				sg =s0*(1.0-s)+s*s1;
		}
		_assert_(on);
		V.r= on->F(sg);
		GV=VertexOnGeom(V,*on,sg);
		return on;
	}
	/*}}}*/
	I2 Geometry::R2ToI2(const R2 & P) const {/*{{{*/
		/*coefIcoor is the coefficient used for integer coordinates:
		 *                       (x-pmin.x)
		 * Icoor x = (2^30 -1) ------------
		 *                          D
		 * where D is the longest side of the domain (direction x or y)
		 * so that (x-pmin.x)/D is in ]0 1[
		 *
		 * coefIcoor = (2^30 -1)/D
		 */
		return  I2( (int) (coefIcoor*(P.x-pmin.x)) ,(int) (coefIcoor*(P.y-pmin.y)) );
	}/*}}}*/
	void Geometry::UnMarkEdges() {/*{{{*/
		for (int i=0;i<nbe;i++) edges[i].SetUnMark();
	}/*}}}*/
}
