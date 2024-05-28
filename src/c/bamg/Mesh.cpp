#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "../shared/shared.h"
#include "./bamgobjects.h"
#include "./det.h"

namespace bamg {

	/*Constructors/Destructors*/
	Mesh::Mesh(BamgGeom* bamggeom,BamgMesh* bamgmesh, BamgOpts* bamgopts):Gh(*(new Geometry())),BTh(*this){ /*{{{*/

		/*Initialize fields*/
		Init(0);

		/*Read Geometry if provided*/
		if(bamggeom->Edges) {
			Gh.ReadGeometry(bamggeom,bamgopts);
			Gh.PostRead();
		}

		/*Read background mesh*/
		ReadMesh(bamgmesh,bamgopts);

		/*Build Geometry if not provided*/
		if(bamggeom->Edges==NULL) {
			/*Recreate geometry if needed*/
			_printf_("WARNING: mesh present but no geometry found. Reconstructing...\n");
			BuildGeometryFromMesh(bamgopts);
			Gh.PostRead(true);
		}

		/*Set integer coordinates*/
		SetIntCoor();

		/*Fill holes and generate mesh properties*/
		ReconstructExistingMesh(bamgopts);
	}
	/*}}}*/
	Mesh::Mesh(int* index,double* x,double* y,int nods,int nels,BamgOpts* bamgopts):Gh(*(new Geometry())),BTh(*this){/*{{{*/

		Init(0);
		ReadMesh(index,x,y,nods,nels,bamgopts);
		SetIntCoor();
		ReconstructExistingMesh(bamgopts);
	}
	/*}}}*/
	Mesh::Mesh(double* x,double* y,int nods,BamgOpts* bamgopts):Gh(*(new Geometry())),BTh(*this){/*{{{*/
		Triangulate(x,y,nods,bamgopts);
	}
	/*}}}*/
	Mesh::Mesh(const Mesh & Tho,const int *flag ,const int *bb,BamgOpts* bamgopts) : Gh(*(new Geometry())), BTh(*this) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Triangles)*/

		  int i,k,itadj;
		  int kt=0;
		  int * kk    = new int [Tho.nbv];
		  long * reft = new long[Tho.nbt];
		  long nbInT =    Tho.TriangleReferenceList(reft);
		  long * refv = new long[Tho.nbv];

		  for (i=0;i<Tho.nbv;i++)
			kk[i]=-1;
		  for (i=0;i<Tho.nbv;i++)
			refv[i]=0;
		  int nbNewBedge =0;
		  //  int nbOldBedge =0;  
		  for (i=0;i<Tho.nbt;i++)
			if(  reft[i] >=0 && flag[i]) 
			  {
				const Triangle & t = Tho.triangles[i];
				kt++;
				kk[Tho.GetId(t[0])]=1;
				kk[Tho.GetId(t[1])]=1;
				kk[Tho.GetId(t[2])]=1;
				itadj=Tho.GetId(t.TriangleAdj(0));
				if (  reft[itadj] >=0 && !flag[itadj])
				  { nbNewBedge++;
					refv[Tho.GetId(t[VerticesOfTriangularEdge[0][0]])]=bb[i];
					refv[Tho.GetId(t[VerticesOfTriangularEdge[0][1]])]=bb[i];
				  }
				itadj=Tho.GetId(t.TriangleAdj(1));
				if (  reft[itadj] >=0 && !flag[itadj])
				  { nbNewBedge++;
					refv[Tho.GetId(t[VerticesOfTriangularEdge[1][0]])]=bb[i];
					refv[Tho.GetId(t[VerticesOfTriangularEdge[1][1]])]=bb[i];}
					itadj=Tho.GetId(t.TriangleAdj(2));
					if (  reft[itadj] >=0 && !flag[itadj])
					  { nbNewBedge++;
						refv[Tho.GetId(t[VerticesOfTriangularEdge[2][0]])]=bb[i];
						refv[Tho.GetId(t[VerticesOfTriangularEdge[2][1]])]=bb[i];}
			  }
		  k=0;
		  for (i=0;i<Tho.nbv;i++){
			  if (kk[i]>=0) kk[i]=k++;
			}
		  _printf_("   number of vertices " << k << ", remove = " << Tho.nbv - k << "\n");
		  _printf_("   number of triangles " << kt << ", remove = " << nbInT-kt << "\n");
		  _printf_("   number of New boundary edge " << nbNewBedge << "\n");
		  long imaxnbv =k;
		  Init(imaxnbv);
		  for (i=0;i<Tho.nbv;i++)
			if (kk[i]>=0) 
			  {
				vertices[nbv] = Tho.vertices[i];
				if (!vertices[nbv].GetReferenceNumber())
				 vertices[nbv].ReferenceNumber = refv[i];
				nbv++;
			  }
		  if (imaxnbv != nbv){
			  delete [] kk;
			  delete [] refv;
			  _error_("imaxnbv != nbv");
		  }
		  for (i=0;i<Tho.nbt;i++)
			if(reft[i] >=0 && flag[i]){
				const Triangle & t = Tho.triangles[i];
				int i0 = Tho.GetId(t[0]);
				int i1 = Tho.GetId(t[1]);
				int i2 = Tho.GetId(t[2]);
				if(i0<0 || i1<0 || i2<0){
					delete [] refv;
					_error_("i0<0 || i1<0 || i2< 0");
				}
				if(i0>=Tho.nbv || i1>=Tho.nbv || i2>=Tho.nbv){
					delete [] refv;
					_error_("i0>=Tho.nbv || i1>=Tho.nbv || i2>=Tho.nbv");
				}
				triangles[nbt] = Triangle(this,kk[i0],kk[i1],kk[i2]);
				triangles[nbt].color = Tho.subdomains[reft[i]].ReferenceNumber; 
				nbt++;           
			  }
		  if (nbt==0 && nbv==0){
			  delete [] refv;
			  _error_("All triangles have been removed");
		  }
		  delete [] kk;
		  delete [] reft;
		  delete [] refv;
		  BuildGeometryFromMesh(bamgopts);
		  Gh.PostRead(); 
		  SetIntCoor();
		  ReconstructExistingMesh(bamgopts);

		  /*Final checks*/
		  _assert_(kt==nbt);
		  _assert_(nbsubdomains);
		  _assert_(subdomains[0].head && subdomains[0].head->link); 
	  }
	/*}}}*/
	Mesh::Mesh(Mesh & Th,Geometry * pGh,Mesh * pBth,long maxnbv_in)/*{{{*/
	  : Gh(*(pGh?pGh:&Th.Gh)), BTh(*(pBth?pBth:this)) {
		  /*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Triangles)*/
		  Gh.NbRef++;
		  maxnbv_in = Max(maxnbv_in,Th.nbv); 
		  long i;
		  // do all the allocation to be sure all the pointer existe

		  Init(maxnbv_in);// to make the allocation 
		  // copy of triangles
		  nbv = Th.nbv;
		  nbt = Th.nbt;
		  nbe = Th.nbe;
		  nbsubdomains = Th.nbsubdomains;
		  nbtout = Th.nbtout;
		  NbVerticesOnGeomVertex = Th.NbVerticesOnGeomVertex;
		  if(NbVerticesOnGeomVertex)
			VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
		  NbVerticesOnGeomEdge = Th.NbVerticesOnGeomEdge;
		  if (NbVerticesOnGeomEdge)
			VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge] ;
		  if (& BTh == & Th.BTh){ // same background 
			  BTh.NbRef++;
			  NbVertexOnBThVertex = Th.NbVertexOnBThVertex;
			  if(NbVertexOnBThVertex)
				VertexOnBThVertex = new VertexOnVertex[NbVertexOnBThVertex];
			  NbVertexOnBThEdge = Th.NbVertexOnBThEdge;
			  if(NbVertexOnBThEdge)
				VertexOnBThEdge = new VertexOnEdge[NbVertexOnBThEdge];
			 }
		  else { // no add on background mesh 
			  BTh.NbRef++;
			  NbVertexOnBThVertex=0;
			  VertexOnBThVertex=0;
			  NbVertexOnBThEdge=0;
			  VertexOnBThEdge=0;
			 }

		  if(nbe)
			edges = new Edge[nbe];
		  if(nbsubdomains)
			subdomains = new SubDomain[nbsubdomains];
		  pmin = Th.pmin;
		  pmax = Th.pmax;
		  coefIcoor = Th.coefIcoor;
		  for(i=0;i<nbt;i++)
			triangles[i].Set(Th.triangles[i],Th,*this);
		  for(i=0;i<nbe;i++)
			edges[i].Set(Th,i,*this);
		  for(i=0;i<nbv;i++)
			vertices[i].Set(Th.vertices[i],Th,*this);
		  for(i=0;i<nbsubdomains;i++)  
			subdomains[i].Set(Th,i,*this);
		  for (i=0;i<NbVerticesOnGeomVertex;i++)
			VerticesOnGeomVertex[i].Set(Th.VerticesOnGeomVertex[i],Th,*this);
		  for (i=0;i<NbVerticesOnGeomEdge;i++)
			VerticesOnGeomEdge[i].Set(Th.VerticesOnGeomEdge[i],Th,*this);
		  quadtree=0;

	  }
	/*}}}*/
	Mesh::Mesh(long imaxnbv,Mesh & BT,BamgOpts* bamgopts,int keepBackVertices) :Gh(BT.Gh),BTh(BT) {/*{{{*/
		this->Init(imaxnbv);
		TriangulateFromGeom1(bamgopts,keepBackVertices);
	}
	/*}}}*/
	Mesh::Mesh(long imaxnbv,Geometry & G,BamgOpts* bamgopts):Gh(G),BTh(*this){/*{{{*/
		Init(imaxnbv);
		TriangulateFromGeom0(bamgopts);
	}
	/*}}}*/
	Mesh::~Mesh() {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Triangles)*/

		if(vertices)             delete [] vertices;
		if(edges)                delete [] edges;
		if(triangles)            delete [] triangles;
		if(quadtree)             delete    quadtree;
		if(orderedvertices)      delete [] orderedvertices;
		if(subdomains)           delete []  subdomains;
		if(VerticesOnGeomEdge)   delete [] VerticesOnGeomEdge;
		if(VerticesOnGeomVertex) delete [] VerticesOnGeomVertex;
		if(VertexOnBThVertex)    delete [] VertexOnBThVertex;
		if(VertexOnBThEdge)      delete [] VertexOnBThEdge;

		if (Gh.NbRef>0) Gh.NbRef--;
		else if (Gh.NbRef==0) delete &Gh;
		if(&BTh != this){
			if (BTh.NbRef>0) BTh.NbRef--;
			else if (BTh.NbRef==0) delete &BTh;
		}
		Init(0); // set all to zero 
	}
	/*}}}*/

	/*IO*/
	void Mesh::ReadMesh(int* index,double* x,double* y,int nods,int nels,BamgOpts* bamgopts){/*{{{*/

		long i1,i2,i3;
		long i;
		Metric M1(1);
		int verbose=0;
		bool* nodeflags=NULL;

		nbv=nods;
		maxnbv=nbv;
		nbt=nels;
		if(bamgopts) verbose=bamgopts->verbose;

		//Vertices
		if (verbose) _printf_("Reading vertices (" << nbv << ")\n");
		vertices=xNew<BamgVertex>(nbv);
		orderedvertices=xNew<BamgVertex*>(nbv);
		for (i=0;i<nbv;i++){
			vertices[i].r.x=x[i];
			vertices[i].r.y=y[i];
			vertices[i].ReferenceNumber=1;
			vertices[i].m=M1;
			vertices[i].color=0;
		}
		maxnbt=2*maxnbv-2; // for filling The Holes and quadrilaterals 

		//Triangles
		if (verbose) _printf_("Reading triangles (" << nbt << ")\n");
		triangles =new Triangle[maxnbt]; //we cannot allocate only nbt triangles since 
		nodeflags=xNew<bool>(nbv);
		for(i=0;i<nbv;i++) nodeflags[i]=false;
		//other triangles will be added for each edge
		for (i=0;i<nbt;i++){
			Triangle & t = triangles[i];
			i1=(long)index[i*3+0]-1; //for C indexing
			i2=(long)index[i*3+1]-1; //for C indexing
			i3=(long)index[i*3+2]-1; //for C indexing
			t=Triangle(this,i1,i2,i3);
			t.color=1;
			nodeflags[i1]=nodeflags[i2]=nodeflags[i3]=true;
		}

		/*Recreate geometry: */
		if (verbose) _printf_("Building Geometry\n");
		BuildGeometryFromMesh();
		if (verbose) _printf_("Completing geometry\n");
		Gh.PostRead();

		/*Check that there is no orphan*/
		bool isorphan=false;
		for(i=0;i<nbv;i++){
			if(!nodeflags[i]){
				_printf_("Vertex " << i+1 << " does not belong to any element\n");
				isorphan=true;
			}
		}
		if(isorphan) _error_("Orphan found in mesh, see ids above");

		/*Clean up*/
		xDelete<bool>(nodeflags);
	}
	/*}}}*/
	void Mesh::ReadMesh(BamgMesh* bamgmesh, BamgOpts* bamgopts){/*{{{*/

		int    verbose=0;
		double Hmin = HUGE_VAL;    // the infinie value
		long   i1,i2,i3;
		long   i,j;
		Metric M1(1);

		/*Check needed pointer*/
		_assert_(bamgmesh);

		if(bamgopts) verbose=bamgopts->verbose;

		nbv=bamgmesh->VerticesSize[0];
		maxnbv=nbv;
		nbt=bamgmesh->TrianglesSize[0];

		//Vertices
		if(bamgmesh->Vertices){
			if(verbose>5) _printf_("      processing Vertices\n");

			vertices=xNew<BamgVertex>(nbv);
			orderedvertices=xNew<BamgVertex*>(nbv);

			for (i=0;i<nbv;i++){
				vertices[i].r.x=bamgmesh->Vertices[i*3+0];
				vertices[i].r.y=bamgmesh->Vertices[i*3+1];
				vertices[i].ReferenceNumber=(long)bamgmesh->Vertices[i*3+2];
				vertices[i].m=M1;
				vertices[i].color=0;
			}
			maxnbt=2*maxnbv-2; // for filling The Holes and quadrilaterals 
		}
		else{
			if(verbose>5) _error_("no Vertices found in the initial mesh");
		}

		//Triangles
		if(bamgmesh->Triangles){
			if(verbose>5) _printf_("      processing Triangles\n");
			triangles =new Triangle[maxnbt]; //we cannot allocate only nbt triangles since 
			//other triangles will be added for each edge
			for (i=0;i<nbt;i++){
				Triangle &t=triangles[i];
				i1=(long)bamgmesh->Triangles[i*4+0]-1; //for C indexing
				i2=(long)bamgmesh->Triangles[i*4+1]-1; //for C indexing
				i3=(long)bamgmesh->Triangles[i*4+2]-1; //for C indexing
				t=Triangle(this,i1,i2,i3);
				t.color=(long)bamgmesh->Triangles[i*4+3];
			}
		}
		else{
			if(verbose>5) _error_("no Triangles found in the initial mesh");
		}

		//VerticesOnGeomEdge
		if(bamgmesh->VerticesOnGeomEdge){
			if(verbose>5) _printf_("      processing VerticesOnGeomEdge\n");
			NbVerticesOnGeomEdge=bamgmesh->VerticesOnGeomEdgeSize[0];
			VerticesOnGeomEdge= new  VertexOnGeom[NbVerticesOnGeomEdge] ;
			for (i=0;i<NbVerticesOnGeomEdge;i++){
				long  i1,i2;
				double s;
				i1=(long)  bamgmesh->VerticesOnGeomEdge[i*3+0]-1; //for C indexing
				i2=(long)  bamgmesh->VerticesOnGeomEdge[i*3+1]-1; //for C indexing
				s =(double)bamgmesh->VerticesOnGeomEdge[i*3+2];
				VerticesOnGeomEdge[i]=VertexOnGeom(vertices[i1],Gh.edges[i2],s);
			}
		}

		//VerticesOnGeomVertex
		if(bamgmesh->VerticesOnGeomVertexSize[0]){
			if(verbose>5) _printf_("      processing VerticesOnGeomVertex\n");
			NbVerticesOnGeomVertex=bamgmesh->VerticesOnGeomVertexSize[0];
			VerticesOnGeomVertex  = new  VertexOnGeom[NbVerticesOnGeomVertex] ;
			for (i=0;i<NbVerticesOnGeomVertex;i++){
				long  i1,i2;
				i1=(long)bamgmesh->VerticesOnGeomVertex[i*2+0]-1; //for C indexing
				i2=(long)bamgmesh->VerticesOnGeomVertex[i*2+1]-1; //for C indexing
				VerticesOnGeomVertex[i]=VertexOnGeom(vertices[i1],Gh.vertices[i2]);
			}
		}

		//Edges
		if (bamgmesh->Edges){
			int i1,i2;
			double* len=NULL;

			if(verbose>5) _printf_("      processing Edges\n");
			nbe=bamgmesh->EdgesSize[0];
			edges= new Edge[nbe];
			//initialize length of each edge (used to provided metric)
			len= new double[nbv];
			for(i=0;i<nbv;i++) len[i]=0;

			for (i=0;i<nbe;i++){
				i1=(int)bamgmesh->Edges[i*3+0]-1; //-1 for C indexing
				i2=(int)bamgmesh->Edges[i*3+1]-1; //-1 for C indexing
				edges[i].ReferenceNumber=(long)bamgmesh->Edges[i*3+2];
				edges[i].v[0]= vertices +i1;
				edges[i].v[1]= vertices +i2;
				edges[i].adj[0]=NULL;
				edges[i].adj[1]=NULL;
				R2 x12=vertices[i2].r-vertices[i1].r;
				double l12=sqrt((x12,x12));

				//prepare metric
				vertices[i1].color++;
				vertices[i2].color++;
				len[i1]+=l12;
				len[i2]+=l12;
				Hmin = Min(Hmin,l12);
			}

			// definition  the default of the given mesh size 
			for (i=0;i<nbv;i++){
				if (vertices[i].color>0) 
				 vertices[i].m=Metric(len[i]/(double)vertices[i].color);
				else 
				 vertices[i].m=Metric(Hmin);
			}
			delete [] len;

			// construction of edges[].adj 
			for (i=0;i<nbv;i++){ 
				vertices[i].color=(vertices[i].color ==2) ?-1:-2;
			}
			for (i=0;i<nbe;i++){
				for (j=0;j<2;j++) { 
					BamgVertex *v=edges[i].v[j];
					long i0=v->color,j0;
					if(i0==-1){
						v->color=i*2+j;
					}
					else if (i0>=0) {// i and i0 edge are adjacent by the vertex v
						j0 = i0%2;
						i0 = i0/2;
						_assert_(v==edges[i0 ].v[j0]);
						edges[i ].adj[j ] =edges +i0;
						edges[i0].adj[j0] =edges +i ;
						v->color = -3;
					}
				}
			}
		}

		//EdgeOnGeomEdge
		if(bamgmesh->EdgesOnGeomEdge){
			if(verbose>5) _printf_("      processing EdgesOnGeomEdge\n");
			int i1,i2,i,j;
			i2=bamgmesh->EdgesOnGeomEdgeSize[0];
			for (i1=0;i1<i2;i1++) {
				i=(int)bamgmesh->EdgesOnGeomEdge[i1*2+0]-1; //C indexing
				j=(int)bamgmesh->EdgesOnGeomEdge[i1*2+1]-1; //C indexing
				//Check value
				if(!(i>=0 && j>=0 && i<nbe && j<Gh.nbe)) {
					_error_("ReadMesh error: EdgesOnGeomEdge edge provided (line " << i1+1 << ": [" << i+1 << " " << j+1 << "]) is incorrect (must be positive, [0<i<nbe=" << nbe << " 0<j<Gh.nbe=" << Gh.nbe << "]");
				}
				edges[i].GeomEdgeHook=Gh.edges+j;
			}
		}

		//SubDomain
		if(bamgmesh->SubDomains){
			long i3,head,direction;
			if(verbose>5) _printf_("      processing SubDomains\n");
			nbsubdomains=bamgmesh->SubDomainsSize[0];
			subdomains = new SubDomain [ nbsubdomains ];
			for (i=0;i<nbsubdomains;i++) {
				i3  =(int)bamgmesh->SubDomains[i*4+0];
				head=(int)bamgmesh->SubDomains[i*4+1]-1;//C indexing
				direction=(int)bamgmesh->SubDomains[i*4+2];
				if (i3!=3) _error_("Bad Subdomain definition: first number should be 3");
				if (head<0 || head>=nbt) _error_("Bad Subdomain definition: head should in [1 " << nbt << "] (triangle number)");
				subdomains[i].head = triangles+head;
				subdomains[i].direction = direction;
				subdomains[i].ReferenceNumber = i3;
			}
		}

	}
	/*}}}*/
	void Mesh::WriteMesh(BamgMesh* bamgmesh,BamgOpts* bamgopts){/*{{{*/

		/*Intermediary*/
		int i,j,k,num,i1,i2;
		long n;
		int* head_1=NULL;
		int* next_1=NULL;
		int* connectivitysize_1=NULL;
		int  connectivitymax_1=0;
		int verbose=0;

		/*Check needed pointer*/
		_assert_(bamgmesh);

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		/*Build reft that holds the number the subdomain number of each triangle, and the real numbering of the elements*/
		long* reft = new long[nbt];
		long* numt = new long[nbt];
		long nbInT = TriangleReferenceList(reft);
		TriangleIntNumbering(numt);

		/*Chaining algorithm used to generate connectivity tables and other outputs*/

		//Memory Allocation
		head_1=xNew<int>(nbv);
		next_1=xNew<int>(3*nbt);
		connectivitysize_1=xNew<int>(nbv);

		//Initialization
		for (i=0;i<nbv;i++) head_1[i]=-1;
		for (i=0;i<nbv;i++) connectivitysize_1[i]=0;
		k=0;
		//Chains generation
		for (i=0;i<nbt;i++) {
			//Do not take into account outside triangles (reft<0)
			if (reft[i]>=0){
				for (j=0;j<3;j++){
					int v=GetId(triangles[i][j]); //jth vertex of the ith triangle
					if (k>3*nbt-1 || k<0) _error_("k = " << k << ", nbt = " << nbt);
					next_1[k]=head_1[v];
					if (v>nbv-1 || v<0)   _error_("v = " << v << ", nbv = " << nbv);
					head_1[v]=k++;
					connectivitysize_1[v]+=1;
				}
			}
		}
		//Get maximum connectivity
		connectivitymax_1=0;
		for (i=0;i<nbv;i++){
			if (connectivitysize_1[i]>connectivitymax_1) connectivitymax_1=connectivitysize_1[i];
		}

		/*OK, now build outputs*/

		/*Vertices*/
		if(verbose>5) _printf_("      writing Vertices\n");
		bamgmesh->VerticesSize[0]=nbv;
		bamgmesh->VerticesSize[1]=3;
		if(nbv){
			bamgmesh->Vertices=xNew<double>(3*nbv);
			bamgmesh->PreviousNumbering=xNew<double>(nbv);
			for (i=0;i<nbv;i++){
				bamgmesh->Vertices[i*3+0]=vertices[i].r.x;
				bamgmesh->Vertices[i*3+1]=vertices[i].r.y;
				bamgmesh->Vertices[i*3+2]=vertices[i].GetReferenceNumber();
				bamgmesh->PreviousNumbering[i]=vertices[i].PreviousNumber;
			}
		}

		/*Edges*/
		if(verbose>5) _printf_("      writing Edges\n");
		bamgmesh->EdgesSize[0]=nbe;
		bamgmesh->EdgesSize[1]=3;
		int NumIssmSegments=0;
		if (nbe){
			bamgmesh->Edges=xNew<double>(3*nbe);
			for (i=0;i<nbe;i++){
				bamgmesh->Edges[i*3+0]=GetId(edges[i][0])+1; //back to M indexing
				bamgmesh->Edges[i*3+1]=GetId(edges[i][1])+1; //back to M indexing
				bamgmesh->Edges[i*3+2]=edges[i].ReferenceNumber;
				if(edges[i].GeomEdgeHook){
					NumIssmSegments++;
				}
			}
		}

		/*Element edges*/
		if(verbose>5) _printf_("      writing element edges\n");
		SetOfEdges4* edge4=new SetOfEdges4(nbt*3,nbv);
		double* elemedge=NULL;
		elemedge=xNew<double>(3*nbt);
		for (i=0;i<3*nbt;i++) elemedge[i]=-2.;//will become -1
		k=0;
		for (i=0;i<nbt;i++){
			//Do not take into account outside triangles (reft<0)
			if (reft[i]>=0){
				for  (j=0;j<3;j++) {
					i1=GetId(triangles[i][VerticesOfTriangularEdge[j][0]]);
					i2=GetId(triangles[i][VerticesOfTriangularEdge[j][1]]);
					n =edge4->SortAndFind(i1,i2);
					if (n==-1){
						//first time
						n=edge4->SortAndAdd(i1,i2);
						elemedge[n*2+0]=double(k);
					}
					else{
						//second time
						elemedge[n*2+1]=double(k);
					}
				}
				k++;
			}
		}
		bamgmesh->IssmEdgesSize[0]=edge4->nb();
		bamgmesh->IssmEdgesSize[1]=4;
		bamgmesh->IssmEdges=xNew<double>(4*edge4->nb());
		for (i=0;i<edge4->nb();i++){
			/*Invert first two vertices if necessary*/
			bool found=false;
			for (j=0;j<3;j++){
				if (triangles[(int)elemedge[2*i+0]](j)==vertices+edge4->i(i)){
					if (triangles[(int)elemedge[2*i+0]]((j+1)%3)==vertices+edge4->j(i)){
						//trigonometric direction
						bamgmesh->IssmEdges[i*4+0]=edge4->i(i)+1;// back to M indexing
						bamgmesh->IssmEdges[i*4+1]=edge4->j(i)+1;// back to M indexing
					}
					else{
						bamgmesh->IssmEdges[i*4+0]=edge4->j(i)+1;// back to M indexing
						bamgmesh->IssmEdges[i*4+1]=edge4->i(i)+1;// back to M indexing
					}
					found=true;
					break;
				}
			}
			_assert_(found);
			bamgmesh->IssmEdges[i*4+2]=elemedge[2*i+0]+1; // back to M indexing
			bamgmesh->IssmEdges[i*4+3]=elemedge[2*i+1]+1; // back to M indexing
		}
		//clean up
		delete edge4;
		xDelete<double>(elemedge);

		/*IssmSegments*/
		if(verbose>5) _printf_("      writing IssmSegments\n");
		bamgmesh->IssmSegmentsSize[0]=NumIssmSegments;
		bamgmesh->IssmSegmentsSize[1]=4;
		bamgmesh->IssmSegments=xNew<double>(4*NumIssmSegments);
		num=0;
		for (i=0;i<nbe;i++){
			if(edges[i].GeomEdgeHook){
				//build segment
				int i1=GetId(edges[i][0]);
				int i2=GetId(edges[i][1]);
				bool stop=false;
				for(j=head_1[i1];j!=-1;j=next_1[j]){
					for(k=0;k<3;k++){
						if (GetId(triangles[(int)j/3][k])==i1){
							if (GetId(triangles[(int)j/3][(int)((k+1)%3)])==i2){
								bamgmesh->IssmSegments[num*4+0]=GetId(edges[i][0])+1; //back to M indexing
								bamgmesh->IssmSegments[num*4+1]=GetId(edges[i][1])+1; //back to M indexing
								bamgmesh->IssmSegments[num*4+2]=(int)j/3+1;            //back to M indexing
								bamgmesh->IssmSegments[num*4+3]=edges[i].ReferenceNumber;
								num+=1;
								stop=true;
								break;
							}
							if (GetId(triangles[(int)j/3][(int)((k+2)%3)])==i2){
								bamgmesh->IssmSegments[num*4+0]=GetId(edges[i][1])+1; //back to M indexing
								bamgmesh->IssmSegments[num*4+1]=GetId(edges[i][0])+1; //back to M indexing
								bamgmesh->IssmSegments[num*4+2]=(int)j/3+1;            //back to M indexing
								bamgmesh->IssmSegments[num*4+3]=edges[i].ReferenceNumber;
								num+=1;
								stop=true;
								break;
							}
						}
					}
					if(stop) break;
				}
				if (!stop){
					_error_("Element holding segment [" << i1+1 << " " << i2+1 << "] not found...");
				}
			}
		}

		/*Triangles*/
		if(verbose>5) _printf_("      writing Triangles\n");
		k=nbInT;
		num=0;
		bamgmesh->TrianglesSize[0]=k;
		bamgmesh->TrianglesSize[1]=4;
		if (k){
			bamgmesh->Triangles=xNew<double>(4*k);
			for (i=0;i<nbt;i++){
				Triangle &t=triangles[i];
				//reft[i]=-1 for outside triangle
				if (reft[i]>=0 && !( t.Hidden(0) || t.Hidden(1) || t.Hidden(2) )){
					bamgmesh->Triangles[num*4+0]=GetId(t[0])+1; //back to M indexing
					bamgmesh->Triangles[num*4+1]=GetId(t[1])+1; //back to M indexing
					bamgmesh->Triangles[num*4+2]=GetId(t[2])+1; //back to M indexing
					bamgmesh->Triangles[num*4+3]=subdomains[reft[i]].ReferenceNumber;
					num=num+1;
				}
			}
		}

		/*SubDomains*/
		if(verbose>5) _printf_("      writing SubDomains\n");
		bamgmesh->SubDomainsSize[0]=nbsubdomains;
		bamgmesh->SubDomainsSize[1]=4;
		if (nbsubdomains){
			bamgmesh->SubDomains=xNew<double>(4*nbsubdomains);
			for (i=0;i<nbsubdomains;i++){
				bamgmesh->SubDomains[i*4+0]=3;
				bamgmesh->SubDomains[i*4+1]=reft[GetId(subdomains[i].head)]+1;//MATLAB indexing
				bamgmesh->SubDomains[i*4+2]=1;
				bamgmesh->SubDomains[i*4+3]=subdomains[i].ReferenceNumber;
			}
		}

		/*SubDomainsFromGeom*/
		if(verbose>5) _printf_("      writing SubDomainsFromGeom\n");
		bamgmesh->SubDomainsFromGeomSize[0]=Gh.nbsubdomains;
		bamgmesh->SubDomainsFromGeomSize[1]=4;
		if (Gh.nbsubdomains){
			bamgmesh->SubDomainsFromGeom=xNew<double>(4*Gh.nbsubdomains);
			for (i=0;i<Gh.nbsubdomains;i++){
				bamgmesh->SubDomainsFromGeom[i*4+0]=2;
				bamgmesh->SubDomainsFromGeom[i*4+1]=GetId(subdomains[i].edge)+1; //back to Matlab indexing
				bamgmesh->SubDomainsFromGeom[i*4+2]=subdomains[i].direction;
				bamgmesh->SubDomainsFromGeom[i*4+3]=Gh.subdomains[i].ReferenceNumber;
			}
		}

		/*VerticesOnGeomVertex*/
		if(verbose>5) _printf_("      writing VerticesOnGeomVertex\n");
		bamgmesh->VerticesOnGeomVertexSize[0]=NbVerticesOnGeomVertex;
		bamgmesh->VerticesOnGeomVertexSize[1]=2;
		if (NbVerticesOnGeomVertex){
			bamgmesh->VerticesOnGeomVertex=xNew<double>(2*NbVerticesOnGeomVertex);
			for (i=0;i<NbVerticesOnGeomVertex;i++){
				VertexOnGeom &v=VerticesOnGeomVertex[i];
				_assert_(v.OnGeomVertex());
				bamgmesh->VerticesOnGeomVertex[i*2+0]=GetId((BamgVertex*)v)+1; //back to Matlab indexing
				bamgmesh->VerticesOnGeomVertex[i*2+1]=Gh.GetId((GeomVertex*)v)+1; //back to Matlab indexing
			}
		}

		/*VertexOnGeomEdge*/
		if(verbose>5) _printf_("      writing VerticesOnGeomEdge\n");
		bamgmesh->VerticesOnGeomEdgeSize[0]=NbVerticesOnGeomEdge;
		bamgmesh->VerticesOnGeomEdgeSize[1]=3;
		if (NbVerticesOnGeomEdge){
			bamgmesh->VerticesOnGeomEdge=xNew<double>(3*NbVerticesOnGeomEdge);
			for (i=0;i<NbVerticesOnGeomEdge;i++){
				const VertexOnGeom &v=VerticesOnGeomEdge[i];
				if (!v.OnGeomEdge()){
					_error_("A vertices supposed to be OnGeomEdge is actually not");
				}
				bamgmesh->VerticesOnGeomEdge[i*3+0]=GetId((BamgVertex*)v)+1; //back to Matlab indexing
				bamgmesh->VerticesOnGeomEdge[i*3+1]=Gh.GetId((const GeomEdge*)v)+1; //back to Matlab indexing
				bamgmesh->VerticesOnGeomEdge[i*3+2]=(double)v; //absisce
			}
		}

		/*EdgesOnGeomEdge*/
		if(verbose>5) _printf_("      writing EdgesOnGeomEdge\n");
		k=0;
		for (i=0;i<nbe;i++){
			if (edges[i].GeomEdgeHook) k=k+1;
		}
		bamgmesh->EdgesOnGeomEdgeSize[0]=k;
		bamgmesh->EdgesOnGeomEdgeSize[1]=2;
		if (k){
			bamgmesh->EdgesOnGeomEdge=xNew<double>(2*(int)k);
			int count=0;
			for (i=0;i<nbe;i++){
				if (edges[i].GeomEdgeHook){
					bamgmesh->EdgesOnGeomEdge[count*2+0]=(double)i+1; //back to Matlab indexing
					bamgmesh->EdgesOnGeomEdge[count*2+1]=(double)Gh.GetId(edges[i].GeomEdgeHook)+1; //back to Matlab indexing
					count=count+1;
				}
			}
		}

		/*Element Connectivity*/
		if(verbose>5) _printf_("      writing Element connectivity\n");
		bamgmesh->ElementConnectivitySize[0]=nbt-nbtout;
		bamgmesh->ElementConnectivitySize[1]=3;
		bamgmesh->ElementConnectivity=xNew<double>(3*(nbt-nbtout));
		for (i=0;i<3*(nbt-nbtout);i++) bamgmesh->ElementConnectivity[i]=NAN;
		num=0;
		for (i=0;i<nbt;i++){
			if (reft[i]>=0){
				for (j=0;j<3;j++){
					k=GetId(triangles[i].TriangleAdj(j));
					if (reft[k]>=0){
						_assert_(3*num+j<3*(nbt-nbtout));
						bamgmesh->ElementConnectivity[3*num+j]=k+1; // back to Matlab indexing
					}
				}
				num+=1;
			}
		}

		/*ElementNodal Connectivity*/
		if(verbose>5) _printf_("      writing Nodal element connectivity\n");
		bamgmesh->NodalElementConnectivitySize[0]=nbv;
		bamgmesh->NodalElementConnectivitySize[1]=connectivitymax_1;
		bamgmesh->NodalElementConnectivity=xNew<double>(connectivitymax_1*nbv);
		for (i=0;i<connectivitymax_1*nbv;i++) bamgmesh->NodalElementConnectivity[i]=NAN;
		for (i=0;i<nbv;i++){
			k=0;
			for(j=head_1[i];j!=-1;j=next_1[j]){
				_assert_(connectivitymax_1*i+k < connectivitymax_1*nbv);
				bamgmesh->NodalElementConnectivity[connectivitymax_1*i+k]=floor((double)j/3)+1;
				k++;
			}
		}

		/*Nodal Connectivity*/
		if(verbose>5) _printf_("      writing Nodal connectivity\n");
		//chaining algorithm (again...)
		int* head_2=NULL;
		int* next_2=NULL;
		int* connectivitysize_2=NULL;
		int  connectivitymax_2=0;
		i1=bamgmesh->IssmEdgesSize[0];
		i2=bamgmesh->IssmEdgesSize[1];
		head_2=xNew<int>(nbv);
		next_2=xNew<int>(2*i1);
		connectivitysize_2=xNew<int>(nbv);
		//Initialization
		for (i=0;i<nbv;i++) head_2[i]=-1;
		for (i=0;i<nbv;i++) connectivitysize_2[i]=0;
		k=0;
		//Chains generation
		for (i=0;i<i1;i++) {
			for (j=0;j<2;j++){
				int v=(int)bamgmesh->IssmEdges[i*i2+j]-1; //back to C indexing
				if (k>2*i1-1 || k<0) _error_("Index exceed matrix dimensions (k=" << k << " not in [0 " << 2*i1-1 << "]");
				next_2[k]=head_2[v];
				if (v>nbv-1 || v<0)   _error_("Index exceed matrix dimensions (v=" << v << " not in [0 " << nbv-1 << "])");
				head_2[v]=k++;
				connectivitysize_2[v]+=1;
			}
		}
		//Get maximum connectivity
		for (i=0;i<nbv;i++){
			if (connectivitysize_2[i]>connectivitymax_2) connectivitymax_2=connectivitysize_2[i];
		}
		//Build output
		connectivitymax_2++;//add last column for size
		bamgmesh->NodalConnectivitySize[0]=nbv;
		bamgmesh->NodalConnectivitySize[1]=connectivitymax_2;
		bamgmesh->NodalConnectivity=xNew<double>(connectivitymax_2*nbv);
		for (i=0;i<connectivitymax_2*nbv;i++) bamgmesh->NodalConnectivity[i]=0;
		for (i=0;i<nbv;i++){
			k=0;
			for(j=head_2[i];j!=-1;j=next_2[j]){
				_assert_(connectivitymax_2*i+k < connectivitymax_2*nbv);
				num=(int)bamgmesh->IssmEdges[int(j/2)*i2+0];
				if (i+1==num){ //carefull, ElementEdge is in M indexing
					//i is the first vertex of the edge, it is therefore connected to the second vertex
					bamgmesh->NodalConnectivity[connectivitymax_2*i+k]=bamgmesh->IssmEdges[int(j/2)*i2+1];
				}
				else{
					bamgmesh->NodalConnectivity[connectivitymax_2*i+k]=num;
				}
				k++;
			}
			bamgmesh->NodalConnectivity[connectivitymax_2*(i+1)-1]=k;
		}

		/*Cracked vertices*/
		if(verbose>5) _printf_("      writing Cracked vertices\n");
		bamgmesh->CrackedVerticesSize[0]=NbCrackedVertices;
		bamgmesh->CrackedVerticesSize[1]=2;
		if (NbCrackedVertices){
			bamgmesh->CrackedVertices=xNew<double>(2*NbCrackedVertices);
			for (i=0;i<NbCrackedVertices;i++){
				bamgmesh->CrackedVertices[i*2+0]=CrackedVertices[i*2+0]+1; //M indexing
				bamgmesh->CrackedVertices[i*2+1]=CrackedVertices[i*2+1]+1; //M indexing
			}
		}

		/*Cracked vertices*/
		if(verbose>5) _printf_("      writing Cracked vertices\n");
		bamgmesh->CrackedEdgesSize[0]=NbCrackedEdges;
		bamgmesh->CrackedEdgesSize[1]=4;
		if (NbCrackedEdges){
			bamgmesh->CrackedEdges=xNew<double>(2*NbCrackedEdges);
			for (i=0;i<NbCrackedEdges;i++){
				bamgmesh->CrackedEdges[i*2+0]=0;//CrackedEdges[i]->+1; //M indexing
				bamgmesh->CrackedEdges[i*2+1]=0;//CrackedEdges[i]-]->+1; //M indexing
			}
		}

		//clean up
		xDelete<int>(connectivitysize_1);
		xDelete<int>(head_1);
		xDelete<int>(next_1);
		xDelete<int>(connectivitysize_2);
		xDelete<int>(head_2);
		xDelete<int>(next_2);
		delete [] reft;
		delete [] numt;
	}
	/*}}}*/
	void Mesh::ReadMetric(const BamgOpts* bamgopts) {/*{{{*/

		/*Intermediary*/
		int  i,j;
		int verbose=0;

		/*Check pointer*/
		_assert_(bamgopts);

		/*Get options*/
		verbose=bamgopts->verbose;

		if(verbose>3) _printf_("      processing metric\n");
		double hmin = Max(bamgopts->hmin,MinimalHmin());
		double hmax = Min(bamgopts->hmax,MaximalHmax());
		double coef = bamgopts->coeff;

		//for now we only use j==3
		j=3;

		for (i=0;i<nbv;i++){
			double h;
			if (j == 1){
				h=bamgopts->metric[i];
				vertices[i].m=Metric(Max(hmin,Min(hmax, h*coef)));
			}
			else if (j==3){
				//do not erase metric computed by hVertices
				if (vertices[i].m.a11==1 && vertices[i].m.a21==0 && vertices[i].m.a22==1){
					double a,b,c;	     
					a=bamgopts->metric[i*3+0];
					b=bamgopts->metric[i*3+1];
					c=bamgopts->metric[i*3+2];
					Metric M(a,b,c);
					EigenMetric Vp(M/coef);

					Vp.Maxh(hmax);
					Vp.Minh(hmin);
					vertices[i].m = Vp;
				}
			}
		}
	}
	/*}}}*/
	void Mesh::WriteMetric(BamgOpts* bamgopts) {/*{{{*/
		int i;
		_assert_(bamgopts);
		xDelete<double>(bamgopts->metric);
		bamgopts->metric=xNew<double>(3*nbv);
		for (i=0;i<nbv;i++){
			bamgopts->metric[i*3+0]=vertices[i].m.a11;
			bamgopts->metric[i*3+1]=vertices[i].m.a21;
			bamgopts->metric[i*3+2]=vertices[i].m.a22;
		}
	}
	/*}}}*/
	void Mesh::WriteIndex(int** pindex,int* pnels){/*{{{*/

		/*Intermediary*/
		int i,k;

		/*output*/
		int* index=NULL;
		int  num=0;

		/*Get number of triangles*/
		k=0;
		for (i=0;i<nbt;i++){
			Triangle &t=triangles[i];
			if(t.det>0) k++;
		}

		if (k){
			index=xNew<int>(3*k);
			for (i=0;i<nbt;i++){
				Triangle &t=triangles[i];
				if (t.det>0 && !(t.Hidden(0)||t.Hidden(1) || t.Hidden(2) )){
					/*Remove triangles that have a bad aspect ratio*/
					//if(t.Anisotropy()<2 & t.Length()<1.e+5){
						index[num*3+0]=GetId(t[0])+1; //back to M indexing
						index[num*3+1]=GetId(t[1])+1; //back to M indexing
						index[num*3+2]=GetId(t[2])+1; //back to M indexing
						num=num+1;
					//}
				}
			}
		}

		/*Assign output pointers*/
		*pindex=index;
		*pnels=num;
	}
	/*}}}*/

	/*Methods*/
	void Mesh::AddMetric(BamgOpts* bamgopts){/*{{{*/
		//  Hessiantype = 0 =>  H is computed using double L2 projection
		//  Hessiantype = 1 =>  H is computed with green formula

		/*Check pointer*/
		_assert_(bamgopts);

		/*Options*/
		int Hessiantype=bamgopts->Hessiantype;

		if (Hessiantype==0){
			BuildMetric0(bamgopts);
		}
		else if (Hessiantype==1){
			BuildMetric1(bamgopts);
		}
		else{
			_error_("Hessiantype " << Hessiantype << " not supported yet (1->use Green formula, 0-> double L2 projection)");
		}
	}
	/*}}}*/
	void Mesh::AddVertex( BamgVertex &s,Triangle* t, long long* det3){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Add)*/
		// -------------------------------
		//             s2
		//                               !
		//             /|\               !
		//            / | \              !
		//           /  |  \             !
		//    tt1   /   |   \ tt0        !
		//         /    |s   \           !
		//        /     .     \          !
		//       /  .      `   \         !
		//      / .           ` \        !
		//      ----------------         !
		//   s0       tt2       s1
		//-------------------------------

		/*Intermediaries*/
		Triangle* tt[3];       //the three triangles
		long long det3local[3];   //three determinants (integer)
		int nbzerodet =0;      //number of zeros in det3
		int izerodet=-1;       //egde containing the vertex s
		int iedge; 

		/*three vertices of t*/
		BamgVertex* s0=t->GetVertex(0);
		BamgVertex* s1=t->GetVertex(1);
		BamgVertex* s2=t->GetVertex(2);

		//determinant of t
		long long detOld=t->det;

		/* infvertexindex = index of the infinite vertex (NULL)
			if no infinite vertex (NULL) infvertexindex=-1
			else if v_i is infinite, infvertexindex=i*/
		int infvertexindex = s0 ? ((s1? (s2?-1:2):1)):0;

		//some checks
		if(((infvertexindex <0 ) && (detOld <0)) ||  ((infvertexindex >=0) && (detOld >0)) ){
			_error_("inconsistent configuration (Contact ISSM developers)");
		}

		// if det3 does not exist, build it 
		if (!det3){ 
			//allocate
			det3 = det3local;
			//if no infinite vertex
			if (infvertexindex<0 ) {
				det3[0]=bamg::det(s .GetIntegerCoordinates(),s1->GetIntegerCoordinates(),s2->GetIntegerCoordinates());
				det3[1]=bamg::det(s0->GetIntegerCoordinates(),s .GetIntegerCoordinates(),s2->GetIntegerCoordinates());
				det3[2]=bamg::det(s0->GetIntegerCoordinates(),s1->GetIntegerCoordinates(),s.GetIntegerCoordinates());}
			else { 
				// one of &s1  &s2  &s0 is NULL
				det3[0]= s0 ? -1 : bamg::det(s.GetIntegerCoordinates(),s1->GetIntegerCoordinates(),s2->GetIntegerCoordinates()) ;
				det3[1]= s1 ? -1 : bamg::det(s0->GetIntegerCoordinates(),s.GetIntegerCoordinates(),s2->GetIntegerCoordinates()) ;
				det3[2]= s2 ? -1 : bamg::det(s0->GetIntegerCoordinates(),s1->GetIntegerCoordinates(),s.GetIntegerCoordinates()) ;
			}
		}

		if (!det3[0]) izerodet=0,nbzerodet++;
		if (!det3[1]) izerodet=1,nbzerodet++;
		if (!det3[2]) izerodet=2,nbzerodet++;

		//if nbzerodet>0, point s is on an egde or on a vertex 
		if  (nbzerodet>0){ 
			/*s is on an edge*/
			if (nbzerodet==1) {
				iedge = OppositeEdge[izerodet];
				AdjacentTriangle ta = t->Adj(iedge);

				/*if the point is one the boundary 
				  add the point in outside part */
				if (t->det>=0){ // inside triangle
					if (((Triangle*)ta)->det<0 ) {
						// add in outside triangle 
						AddVertex(s,( Triangle *)ta);
						return;
					}
				}
			}
			else{
				_error_("Cannot add a vertex more than once. Check duplicates");
			}
		}

		// remove de MarkUnSwap edge
		t->SetUnMarkUnSwap(0);
		t->SetUnMarkUnSwap(1);
		t->SetUnMarkUnSwap(2);

		tt[0]= t;
		tt[1]= &triangles[nbt++];
		tt[2]= &triangles[nbt++];

		if (nbt>maxnbt) _error_("Not enough triangles");

		*tt[1]=*tt[2]=*t;
		tt[0]->link=tt[1];
		tt[1]->link=tt[2]; 

		(*tt[0])(OppositeVertex[0])=&s;
		(*tt[1])(OppositeVertex[1])=&s;
		(*tt[2])(OppositeVertex[2])=&s;

		tt[0]->det=det3[0];
		tt[1]->det=det3[1];
		tt[2]->det=det3[2];         

		//  update adj des triangles externe 
		tt[0]->SetAdjAdj(0);
		tt[1]->SetAdjAdj(1);
		tt[2]->SetAdjAdj(2);
		//  update des adj des 3 triangle interne
		const int i0 = 0;
		const int i1= NextEdge[i0];
		const int i2 = PreviousEdge[i0];

		tt[i0]->SetAdj2(i2,tt[i2],i0);
		tt[i1]->SetAdj2(i0,tt[i0],i1);
		tt[i2]->SetAdj2(i1,tt[i1],i2);

		tt[0]->SetSingleVertexToTriangleConnectivity();
		tt[1]->SetSingleVertexToTriangleConnectivity();
		tt[2]->SetSingleVertexToTriangleConnectivity();

		// swap if the point s is on a edge
		if(izerodet>=0) {
			int rswap=tt[izerodet]->swap(iedge);

			if (!rswap) {
				_error_("swap the point s is on a edge");
			}
		}

	}/*}}}*/
	void Mesh::BoundAnisotropy(BamgOpts* bamgopts,double anisomax,double hminaniso) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/BoundAnisotropy)*/

		int verbose=0;
		double lminaniso = 1./ (Max(hminaniso*hminaniso,1e-100));

		//Get options
		if(bamgopts) verbose=bamgopts->verbose;

		//display info
		if (verbose > 1)  _printf_("   BoundAnisotropy by " << anisomax << "\n");

		double h1=1.e30,h2=1e-30;
		double coef = 1./(anisomax*anisomax);
		double hn1=1.e30,hn2=1e-30,rnx =1.e-30,rx=0;  

		//loop over all vertices
		for (int i=0;i<nbv;i++){
			EigenMetric Vp(vertices[i]);
			double lmax=Vp.lmax();
			Vp*=Min(lminaniso,lmax)/lmax;
			Vp.BoundAniso2(coef);
			vertices[i].m = Vp;

			//info to be displayed
			if (verbose>2){
				h1 =Min(h1,Vp.lmin());
				h2 =Max(h2,Vp.lmax());
				hn1=Min(hn1,Vp.lmin());
				hn2=Max(hn2,Vp.lmax());
				rx =Max(rx,Vp.Aniso2());
				rnx= Max(rnx,Vp.Aniso2());
			}
		}

		//display info
		if (verbose>2){
			_printf_("      input:  Hmin = " << pow(h2,-0.5)  << ", Hmax = " << pow(h1,-0.5) << ", factor of anisotropy max  = " << pow(rx,0.5) << "\n");
			_printf_("      output: Hmin = " << pow(hn2,-0.5) << ", Hmax = " << pow(hn1,-0.5)<< ", factor of anisotropy max  = " <<pow(rnx,0.5) << "\n");
		}
	}
	/*}}}*/
	void Mesh::BuildGeometryFromMesh(BamgOpts* bamgopts){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/ConsGeometry)*/

		/*Reconstruct Geometry from Mesh*/

		/*Intermediary*/
		int i,j,k,kk,it,jt;
		int    verbose=0;

		/*Recover options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//display info
		if (verbose>1) _printf_("   construction of the geometry from the 2d mesh\n");

		//check that the mesh is not empty
		if(nbt<=0 || nbv <=0 ) _error_("nbt or nbv is negative (Mesh empty?)");

		/*Construction of the edges*/

		/*initialize st and edge4*/
		SetOfEdges4* edge4= new SetOfEdges4(nbt*3,nbv);
		long*        st   = new long[nbt*3];

		/*initialize st as -1 (chaining algorithm)*/
		for (i=0;i<nbt*3;i++) st[i]=-1;

		/*build edge4 (chain)*/
		for (i=0;i<nbe;i++){
			edge4->SortAndAdd(GetId(edges[i][0]),GetId(edges[i][1]));
		}
		/*check that there is no double edge*/
		if (nbe !=  edge4->nb()){ 
			delete [] st;
			_error_("Some Double edge in the mesh, the number is " << nbe << ", nbe4=" << edge4->nb()); 
		}
		//keep nbe in nbeold
		long nbeold = nbe;

		//Go through the triangles and add the edges in edge4 if they are not there yet
		for (i=0;i<nbt;i++){
			//3 edges per triangle
			for  (j=0;j<3;j++) {
				//Add Edge to edge4 (k=numberofedges in edge4)
				long k =edge4->SortAndAdd(GetId(triangles[i][VerticesOfTriangularEdge[j][0]]), GetId(triangles[i][VerticesOfTriangularEdge[j][1]]));
				long invisible = triangles[i].Hidden(j);

				//if st[k] has not been changed yet, add 3*i+j (= vertex position in the index)
				if(st[k]==-1) st[k]=3*i+j;

				//else st[k]>=0 -> the edge already exist, check
				else if(st[k]>=0) {
					//check that it is not an edge on boundary (should not already exist)
					if (triangles[i].TriangleAdj(j) || triangles[st[k]/3].TriangleAdj((int) (st[k]%3))){
						_error_("problem in Geometry reconstruction: an edge on boundary is duplicated (double element?)");
					}
					/*OK, the element is not on boundary, is belongs to 2 triangles -> build Adjacent triangles list*/
					triangles[i].SetAdj2(j,triangles + st[k] / 3,(int) (st[k]%3));
					if (invisible)  triangles[i].SetHidden(j);
					/* if k < nbe mark the edge as on Boundary (Locked)*/
					if (k<nbe) {
						triangles[i].SetLocked(j);
					}
					/*set st[k] as negative so that the edge should not be called again*/
					st[k]=-2-st[k]; 
				}
				else {
					/*else (see 3 lines above), the edge has been called more than twice: return error*/
					_printf_("The edge (" << GetId(triangles[i][VerticesOfTriangularEdge[j][0]]) << "," << GetId(triangles[i][VerticesOfTriangularEdge[j][1]]) << ") belongs to more than 2 triangles (" << k << ")\n");
					_printf_("Edge " << j << " of triangle " << i << "\n");
					_printf_("Edge " << (-st[k]+2)%3 << " of triangle " << (-st[k]+2)/3 << "\n");
					_printf_("Edge " << triangles[(-st[k]+2)/3].NuEdgeTriangleAdj((int)((-st[k]+2)%3)) << " of triangle " << GetId(triangles[(-st[k]+2)/3].TriangleAdj((int)((-st[k]+2)%3))) << "\n");
					_error_("An edge belongs to more than 2 triangles");
				}	
			}
		}

		//delete edge4
		long nbedges = edge4->nb(); // the total number of edges 
		delete edge4; edge4=NULL;

		//display info
		if(verbose>5) {
			_printf_("         info on Mesh:\n");
			_printf_("            - number of vertices    = " << nbv << "\n"); 
			_printf_("            - number of triangles   = " << nbt << "\n"); 
			_printf_("            - number of given edges = " << nbe << "\n"); 
			_printf_("            - number of all edges   = " << nbedges << "\n"); 
			_printf_("            - Euler number 1 - nb of holes = " << nbt-nbedges+nbv << "\n"); 
		}

		/*check consistency of edge[].adj and geometrical required  vertices*/
		k=0; kk=0;
		for (i=0;i<nbedges;i++){
			//internal edge
			if (st[i] <-1) {
				//get triangle number back
				it =  (-2-st[i])/3;
				//get edge position back
				j  =  (int) ((-2-st[i])%3);
				Triangle &tt=*triangles[it].TriangleAdj(j);
				if (triangles[it].color != tt.color|| i < nbeold) k++;
			}
			//boundary edge (alone)
			else if (st[i] >=0) 
			 kk++;
		}

		/*Constructions of edges*/
		k += kk;
		kk=0;
		if (k) {
			nbe = k;
			Edge* edgessave=edges;
			edges = new Edge[nbe];
			k =0;

			//display info
			if(verbose>4) _printf_("   Construction of the edges " << nbe << "\n");

			for (i=0;i<nbedges;i++){ 
				long  add= -1;

				//internal edge (belongs to two triangles)
				if (st[i] <-1){ 
					it =  (-2-st[i])/3;
					j  =  (int) ((-2-st[i])%3);
					Triangle & tt = * triangles[it].TriangleAdj(j);
					if (triangles[it].color !=  tt.color || i < nbeold) add=k++;
				}
				//boundary edge
				else if (st[i] >=0){
					it = st[i]/3;
					j  = (int) (st[i]%3);
					add=k++;
				}
				if (add>=0 && add < nbe){
					edges[add].v[0] = &triangles[it][VerticesOfTriangularEdge[j][0]];
					edges[add].v[1] = &triangles[it][VerticesOfTriangularEdge[j][1]];
					edges[add].GeomEdgeHook=NULL; 
					//if already existed
					if (i<nbeold){
						edges[add].ReferenceNumber=edgessave[i].ReferenceNumber; 		      
						edges[add].GeomEdgeHook=edgessave[i].GeomEdgeHook; //  HACK to get required edges
						_printf_("oh no...\n");
					}
					else
					 edges[add].ReferenceNumber=Min(edges[add].v[0]->GetReferenceNumber(),edges[add].v[1]->GetReferenceNumber());
				  }
			}

			//check that we have been through all edges
			if (k!=nbe){
				_error_("problem in edge construction process: k!=nbe (should not happen)");
			}
			//delete edgessave
			if (edgessave) delete [] edgessave;
		}

		/*Color the vertices*/

		//initialize color of all vertices as 0
		for (i=0;i<nbv;i++) vertices[i].color =0;

		//go through the edges and add a color to corresponding vertices
		//(A vertex in 4 edges will have a color 4)
		for (i=0;i<nbe;i++){
		 for (j=0;j<2;j++) edges[i].v[j]->color++;
		}

		//change the color: if a vertex belongs to 2 edges -1, else -2
		for (i=0;i<nbv;i++) {
			vertices[i].color=(vertices[i].color ==2)? -1 : -2;
		}

		/*Build edges[i].adj: adjacency of each edge (if on the same curve)*/
		for (i=0;i<nbe;i++){
			for (j=0;j<2;j++){ 
				//get current vertex
				BamgVertex* v=edges[i].v[j];
				//get vertex color (i0)
				long i0=v->color;
				long j0;

				//if color<0 (first time), no adjacent edge
				if(i0<0) edges[i].adj[j]=NULL;

				//if color=-1 (corner),change the vertex color as 2*i+j (position of the vertex in edges)
				if(i0==-1) v->color=i*2+j;

				//if color>=0 (i and i0 edge are adjacent by the vertex v)
				else if (i0>=0) {
					//get position of v in edges back
					j0 =  i0%2; //column in edges
					i0 =  i0/2; //line in edges

					//check that we have the correct vertex
					if (v!=edges[i0 ].v[j0]){
						_error_("v!=edges[i0 ].v[j0]: this should not happen as the vertex belongs to this edge");
					}

					//Add adjacence
					edges[i ].adj[j ]=edges +i0;
					edges[i0].adj[j0]=edges +i ;

					//change color to -3
					v->color = -3;
				}
			}
		}

		/*Reconstruct subdomains info*/

		//check that nbsubdomains is empty
		if(nbsubdomains) _error_("nbsubdomains should be 0");
		nbsubdomains=0;

		//color the subdomains
		long* colorT= new long[nbt];
		Triangle *tt,*t;

		//initialize the color of each triangle as -1
		for (it=0;it<nbt;it++) colorT[it]=-1;

		//loop over the triangles
		for (it=0;it<nbt;it++){

			//if the triangle has not been colored yet:
			if (colorT[it]<0){

				//color = number of subdomains
				colorT[it]=nbsubdomains;

				//color all the adjacent triangles of T that share a non marked edge
				int level =1;
				int kolor=triangles[it].color;
				st[0]=it; // stack 
				st[1]=0;
				k=1;
				while (level>0){
					if( (j=st[level]++)<3 ){ 
						t = &triangles[st[level-1]];
						tt=t->TriangleAdj((int)j);

						//color the adjacent triangle
						if ( ! t->Locked(j) && tt && (colorT[jt = GetId(tt)] == -1) && ( tt->color==kolor)) {
							colorT[jt]=nbsubdomains;
							st[++level]=jt;
							st[++level]=0;
							k++;
						}
					}
					else level-=2;
				}
				nbsubdomains++;
			}
		}
		if (verbose> 3) _printf_("      The Number of sub domain = " << nbsubdomains << "\n"); 

		//build subdomains
		long isd;
		subdomains = new SubDomain[nbsubdomains];

		//initialize subdomains[isd].head as 0
		for (isd=0;isd<nbsubdomains;isd++) subdomains[isd].head =0;

		k=0;
		for (it=0;it<nbt;it++){
			for (int j=0;j<3;j++){
				tt=triangles[it].TriangleAdj(j);
				if ((!tt || tt->color != triangles[it].color) && !subdomains[isd=colorT[it]].head){
					subdomains[isd].head = triangles+it;
					subdomains[isd].ReferenceNumber =  triangles[it].color;
					subdomains[isd].direction = j; // hack
					subdomains[isd].edge = 0;
					k++;
				}
			}
		}
		//check that we have been through all subdomains
		if (k!= nbsubdomains){
			delete [] colorT;
			_error_("k!= nbsubdomains");
		}
		//delete colorT and st
		delete [] colorT;
		delete [] st;

		/*Reconstruct Geometry Gh*/

		//build colorV -1 for all vertex and 0 for the vertices belonging to edges
		long* colorV = new long[nbv];
		for (i=0;i<nbv;i++) colorV[i]=-1;
		for (i=0;i<nbe;i++){
		 for ( j=0;j<2;j++) colorV[GetId(edges[i][j])]=0;
		}
		//number the vertices belonging to edges
		k=0;
		for (i=0;i<nbv;i++){
		 if(!colorV[i]) colorV[i]=k++;
		}

		//Build Gh
		Gh.nbv=k;
		Gh.nbe = nbe;
		Gh.vertices = new GeomVertex[k];
		Gh.edges = new GeomEdge[nbe];
		Gh.nbsubdomains = nbsubdomains;
		Gh.subdomains = new GeomSubDomain[nbsubdomains];
		if (verbose>3) _printf_("   number of vertices = " << Gh.nbv << "\n   number of edges = " << Gh.nbe << "\n");
		NbVerticesOnGeomVertex = Gh.nbv;
		VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];
		NbVerticesOnGeomEdge =0;
		VerticesOnGeomEdge =0;

		//Build VertexOnGeom
		for (i=0;i<nbv;i++){
			if((j=colorV[i])>=0){
				BamgVertex & v = Gh.vertices[j];
				v = vertices[i];
				v.color =0;
				VerticesOnGeomVertex[j] = VertexOnGeom(vertices[i], Gh.vertices[j]);
			}
		}

		//Buid pmin and pmax of Gh (extrema coordinates)
		Gh.pmin =  Gh.vertices[0].r;
		Gh.pmax =  Gh.vertices[0].r;
		// recherche des extrema des vertices pmin,pmax
		for (i=0;i<Gh.nbv;i++) {
			Gh.pmin.x = Min(Gh.pmin.x,Gh.vertices[i].r.x);
			Gh.pmin.y = Min(Gh.pmin.y,Gh.vertices[i].r.y);
			Gh.pmax.x = Max(Gh.pmax.x,Gh.vertices[i].r.x);
			Gh.pmax.y = Max(Gh.pmax.y,Gh.vertices[i].r.y);
		}
		R2 DD05 = (Gh.pmax-Gh.pmin)*0.05;
		Gh.pmin -=  DD05;
		Gh.pmax +=  DD05;

		//Build Gh.coefIcoor
		long MaxICoord = 1073741823; //2^30 - 1 = =111...111 (29 times one)
		Gh.coefIcoor= (MaxICoord)/(Max(Gh.pmax.x-Gh.pmin.x,Gh.pmax.y-Gh.pmin.y));
		if (Gh.coefIcoor<=0){
			delete [] colorV;
			_error_("Gh.coefIcoor<=0 in infered Geometry (this should not happen)");
		}

		/*Build Gh.edges*/

		//initialize len as 0
		double * len = new double[Gh.nbv];
		for(i=0;i<Gh.nbv;i++) len[i]=0;

		//initialize edge4 again
		edge4= new SetOfEdges4(nbe,nbv);  
		double hmin = HUGE_VAL;
		int kreq=0;
		for (i=0;i<nbe;i++){

			long i0 = GetId(edges[i][0]);
			long i1 = GetId(edges[i][1]);
			long j0 = colorV[i0];
			long j1 = colorV[i1];

			Gh.edges[i].v[0] = Gh.vertices +  j0;
			Gh.edges[i].v[1] = Gh.vertices +  j1;

			Gh.edges[i].type = 0;

			Gh.edges[i].tg[0]=R2();
			Gh.edges[i].tg[1]=R2();

			bool required= edges[i].GeomEdgeHook; 
			if(required) kreq++;
			edges[i].GeomEdgeHook =  Gh.edges + i;
			if(required){
				Gh.edges[i].v[0]->SetRequired();
				Gh.edges[i].v[1]->SetRequired();
				Gh.edges[i].SetRequired();
			}

			R2 x12 = Gh.vertices[j0].r-Gh.vertices[j1].r;
			double l12=Norme2(x12);        
			hmin = Min(hmin,l12);

			Gh.vertices[j1].color++;
			Gh.vertices[j0].color++;

			len[j0]+= l12;
			len[j1] += l12;
			hmin = Min(hmin,l12);
			Gh.edges[i].ReferenceNumber  = edges[i].ReferenceNumber;

			k = edge4->SortAndAdd(i0,i1);
			if (k != i){
				delete [] len;
				delete [] colorV;
				_error_("problem in Edge4 construction: k != i");
			}
		}

		//Build metric for all vertices of Gh
		for (i=0;i<Gh.nbv;i++){
		 if (Gh.vertices[i].color > 0) 
		  Gh.vertices[i].m=  Metric(len[i] /(double) Gh.vertices[i].color);
		 else 
		  Gh.vertices[i].m=  Metric(hmin);
		}
		//delete len
		delete [] len;

		//Build Gh.subdomains
		for (i=0;i<nbsubdomains;i++){
			it = GetId(subdomains[i].head);
			j = subdomains[i].direction;
			long i0 = GetId(triangles[it][VerticesOfTriangularEdge[j][0]]);
			long i1 = GetId(triangles[it][VerticesOfTriangularEdge[j][1]]);
			k = edge4->SortAndFind(i0,i1);
			if(k<0){
				delete [] colorV;
				_error_("This should not happen");
			}
			subdomains[i].direction = (vertices + i0 == edges[k].v[0]) ? 1 : -1;
			subdomains[i].edge = edges+k;
			Gh.subdomains[i].edge = Gh.edges + k;
			Gh.subdomains[i].direction  =  subdomains[i].direction;
			Gh.subdomains[i].ReferenceNumber =  subdomains[i].ReferenceNumber;
		}

		delete edge4;
		delete [] colorV;

		//unset adj
		for (i=0;i<nbt;i++){
			for (j=0;j<3;j++){
				triangles[i].SetAdj2(j,0,triangles[i].GetAllflag(j));
			}
		}

	}
	/*}}}*/
	void Mesh::BuildMetric0(BamgOpts* bamgopts){/*{{{*/

		/*Options*/
		double* s=NULL;
		long    nbsol;
		int     verbose=0;

		int   i,j,k,iA,iB,iC;
		int   iv;

		/*Check pointer*/
		_assert_(bamgopts);

		/*Recover options*/
		verbose=bamgopts->verbose;

		/*Get and process fields*/
		s=bamgopts->field;
		nbsol=bamgopts->fieldSize[1];

		/*Check size*/
		if (bamgopts->fieldSize[0] != nbv) _error_("'field' should have " << nbv << " rows");

		//initialization of some variables
		double* ss=(double*)s;
		double  sA,sB,sC;
		double*  detT = new double[nbt];
		double*  sumareas = new double[nbv];
		double*  alpha= new double[nbt*3];
		double*  beta = new double[nbt*3];
		double*  dx_elem    = new double[nbt];
		double*  dy_elem    = new double[nbt];
		double*  dx_vertex  = new double[nbv];
		double*  dy_vertex  = new double[nbv];
		double*  dxdx_elem  = new double[nbt];
		double*  dxdy_elem  = new double[nbt];
		double*  dydy_elem  = new double[nbt];
		double*  dxdx_vertex= new double[nbv];
		double*  dxdy_vertex= new double[nbv];
		double*  dydy_vertex= new double[nbv];

		//display infos
		if(verbose>1) {
			_printf_("   Construction of Metric: number of field: " << nbsol << " (nbt=" << nbt << ", nbv=" << nbv << ")\n");
		}

		//first, build the chains that will be used for the Hessian computation, as weel as the area of each element
		int* head_s=NULL;
		head_s=xNew<int>(nbv);
		int* next_p=NULL;
		next_p=xNew<int>(3*nbt);
		int  p=0;
		//initialization
		for(i=0;i<nbv;i++){
			sumareas[i]=0;
			head_s[i]=-1;
		}
		for(i=0;i<nbt;i++){

			//lopp over the real triangles (no boundary elements)
			if(triangles[i].link){ 

				//get current triangle t
				const Triangle &t=triangles[i];

				// coor of 3 vertices 
				R2 A=t[0];
				R2 B=t[1];
				R2 C=t[2];

				//compute triangle determinant (2*Area)
				double dett = bamg::Area2(A,B,C);
				detT[i]=dett;

				/*The nodal functions are such that for a vertex A:
				 *    N_A(x,y)=alphaA x + beta_A y +gamma_A
				 *    N_A(A) = 1,   N_A(B) = 0,   N_A(C) = 0
				 * solving this system of equation (determinant = 2Area(T) != 0 if A,B and C are not inlined)
				 * leads to:
				 *    N_A = (xB yC - xC yB + x(yB-yC) +y(xC-xB))/(2*Area(T))
				 * and this gives:
				 *    alpha_A = (yB-yC)/(2*Area(T))*/
				alpha[i*3+0]=(B.y-C.y)/dett;
				alpha[i*3+1]=(C.y-A.y)/dett;
				alpha[i*3+2]=(A.y-B.y)/dett;
				beta[ i*3+0]=(C.x-B.x)/dett;
				beta[ i*3+1]=(A.x-C.x)/dett;
				beta[ i*3+2]=(B.x-A.x)/dett;

				//compute chains
				for(j=0;j<3;j++){
					k=GetId(triangles[i][j]);
					next_p[p]=head_s[k];
					head_s[k]=p++;

					//add area to sumareas
					sumareas[k]+=dett;
				}

			}
		}

		/*Check err, we need to make sure it has the right size!*/
		bool    deleteerr = false;
		double* err       = NULL;
		if(bamgopts->errSize[0]==1){
			/*Let's copy this value for all vertices*/
			err =  new double[nbv*bamgopts->errSize[1]];
			for(int i=0;i<nbv;i++){
				for(int j=0;j<bamgopts->errSize[1];j++){
					err[i*bamgopts->errSize[1]+j] = bamgopts->err[j];
				}
			}
			deleteerr = true;
		}
		else if(bamgopts->errSize[0]==this->nbv){
			/*Nothing to do, already right size*/
			err = bamgopts->err;
		}
		else{
			_error_("number of rows in 'err' not supported: size "<<bamgopts->errSize[0]<<"x"<<bamgopts->errSize[1]<<" (nbv is "<<this->nbv<<")");
		}

		/*for all Solutions*/
		for (int nusol=0;nusol<nbsol;nusol++) {
			double smin=ss[nusol],smax=ss[nusol];

			//get min(s), max(s) and initialize Hessian (dxdx,dxdy,dydy)
			for ( iv=0,k=0; iv<nbv; iv++){
				smin=Min(smin,ss[iv*nbsol+nusol]);
				smax=Max(smax,ss[iv*nbsol+nusol]);
			}
			double sdelta=smax-smin;
			double absmax=Max(Abs(smin),Abs(smax));

			//display info
			if(verbose>2) _printf_("      Solution " << nusol << ", Min = " << smin << ", Max = " << smax << ", Delta = " << sdelta << "\n");

			//skip constant field
			if (sdelta < 1.0e-10*Max(absmax,1e-20)){
				_printf_("      Solution " << nusol << " is constant, skipping...\n");
				continue;
			}

			//initialize the hessian and gradient matrices
			for ( iv=0,k=0; iv<nbv; iv++) dxdx_vertex[iv]=dxdy_vertex[iv]=dydy_vertex[iv]=dx_vertex[iv]=dy_vertex[iv]=0;

			//1: Compute gradient for each element (exact)
			for (i=0;i<nbt;i++){
				if(triangles[i].link){
					// number of the 3 vertices
					iA = GetId(triangles[i][0]);
					iB = GetId(triangles[i][1]);
					iC = GetId(triangles[i][2]);

					// value of the P1 fonction on 3 vertices 
					sA = ss[iA*nbsol+nusol];
					sB = ss[iB*nbsol+nusol];
					sC = ss[iC*nbsol+nusol];

					//gradient = (sum alpha_i s_i, sum_i beta_i s_i)
					dx_elem[i]=sA*alpha[3*i+0]+sB*alpha[3*i+1]+sC*alpha[3*i+2];
					dy_elem[i]=sA*beta[ 3*i+0]+sB*beta[ 3*i+1]+sC*beta[ 3*i+2];
				}
			}

			//2: then compute a gradient for each vertex using a P2 projection
			for(i=0;i<nbv;i++){
				for(p=head_s[i];p!=-1;p=next_p[p]){
					//Get triangle number
					k=(long)(p/3);
					dx_vertex[i]+=dx_elem[k]*detT[k]/sumareas[i];
					dy_vertex[i]+=dy_elem[k]*detT[k]/sumareas[i];
				}
			}

			//3: compute Hessian matrix on each element
			for (i=0;i<nbt;i++){
				if(triangles[i].link){
					// number of the 3 vertices
					iA = GetId(triangles[i][0]);
					iB = GetId(triangles[i][1]);
					iC = GetId(triangles[i][2]);

					//Hessian
					dxdx_elem[i]=dx_vertex[iA]*alpha[3*i+0]+dx_vertex[iB]*alpha[3*i+1]+dx_vertex[iC]*alpha[3*i+2];
					dxdy_elem[i]=dy_vertex[iA]*alpha[3*i+0]+dy_vertex[iB]*alpha[3*i+1]+dy_vertex[iC]*alpha[3*i+2];
					dydy_elem[i]=dy_vertex[iA]*beta[3*i+0]+dy_vertex[iB]*beta[3*i+1]+dy_vertex[iC]*beta[3*i+2];
				}
			}

			//4: finaly compute Hessian on each vertex using the second P2 projection
			for(i=0;i<nbv;i++){
				for(p=head_s[i];p!=-1;p=next_p[p]){
					//Get triangle number
					k=(long)(p/3);
					dxdx_vertex[i]+=dxdx_elem[k]*detT[k]/sumareas[i];
					dxdy_vertex[i]+=dxdy_elem[k]*detT[k]/sumareas[i];
					dydy_vertex[i]+=dydy_elem[k]*detT[k]/sumareas[i];
				}
			}

			/*Compute Metric from Hessian*/
			for(iv=0;iv<nbv;iv++){
				vertices[iv].MetricFromHessian(dxdx_vertex[iv],dxdy_vertex[iv],dydy_vertex[iv],smin,smax,ss[iv*nbsol+nusol],err[iv*nbsol+nusol],bamgopts);
			}

		}//for all solutions

		//clean up
		xDelete<int>(head_s);
		xDelete<int>(next_p);
		delete [] detT;
		delete [] alpha;
		delete [] beta;
		delete [] sumareas;
		delete [] dx_elem;
		delete [] dy_elem;
		delete [] dx_vertex;
		delete [] dy_vertex;
		delete [] dxdx_elem;
		delete [] dxdy_elem;
		delete [] dydy_elem;
		delete [] dxdx_vertex;
		delete [] dxdy_vertex;
		delete [] dydy_vertex;
		if(deleteerr) delete [] err;
	}
	/*}}}*/
	void Mesh::BuildMetric1(BamgOpts* bamgopts){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/IntersectConsMetric)*/

		/*Options*/
		double* s=NULL;
		long nbsol;
		int NbJacobi;
		int verbose;

		/*Check pointer*/
		_assert_(bamgopts);

		/*Recover options*/
		verbose=bamgopts->verbose;
		NbJacobi=bamgopts->nbjacobi;

		/*Get and process fields*/
		s=bamgopts->field;
		nbsol=bamgopts->fieldSize[1];

		/*Check size*/
		if (bamgopts->fieldSize[0] != nbv) _error_("'field' should have " << nbv << " rows");

		//initialization of some variables
		long    i,k,iA,iB,iC,iv;
		R2      O(0,0);
		double* ss=(double*)s;
		double  sA,sB,sC;
		double*  detT = new double[nbt];
		double*  Mmass= new double[nbv];
		double*  Mmassxx= new double[nbv];
		double*  dxdx= new double[nbv];
		double*  dxdy= new double[nbv];
		double*  dydy= new double[nbv];
		double*  workT= new double[nbt];
		double*  workV= new double[nbv];
		int*    OnBoundary = new int[nbv];

		//display infos
		if(verbose>1) {
			_printf_("   Construction of Metric: number of field: " << nbsol << " (nbt=" << nbt << ", nbv=" << nbv << ")\n");
		}

		//initialize Mmass, OnBoundary and Massxx by zero
		for (iv=0;iv<nbv;iv++){
			Mmass[iv]=0;
			OnBoundary[iv]=0;
			Mmassxx[iv]=0;
		}

		//Build detT Mmas Mmassxx workT and OnBoundary
		for (i=0;i<nbt;i++){ 

			//lopp over the real triangles (no boundary elements)
			if(triangles[i].link){ 

				//get current triangle t
				const Triangle &t=triangles[i];

				// coor of 3 vertices 
				R2 A=t[0];
				R2 B=t[1];
				R2 C=t[2];

				// number of the 3 vertices
				iA = GetId(t[0]);
				iB = GetId(t[1]);
				iC = GetId(t[2]);

				//compute triangle determinant (2*Area)
				double dett = bamg::Area2(A,B,C);
				detT[i]=dett;
				dett /= 6;

				// construction of OnBoundary (flag=1 if on boundary, else 0)
				int nbb=0;
				for(int j=0;j<3;j++){
					//get adjacent triangle
					Triangle *ta=t.Adj(j);
					//if there is no adjacent triangle, the edge of the triangle t is on boundary
					if ( !ta || !ta->link){
						//mark the two vertices of the edge as OnBoundary
						OnBoundary[GetId(t[VerticesOfTriangularEdge[j][0]])]=1;
						OnBoundary[GetId(t[VerticesOfTriangularEdge[j][1]])]=1;
						nbb++;
					}
				}

				//number of vertices on boundary for current triangle t
				workT[i] = nbb;

				//Build Mmass Mmass[i] = Mmass[i] + Area/3
				Mmass[iA] += dett;
				Mmass[iB] += dett;
				Mmass[iC] += dett;

				//Build Massxx = Mmass
				Mmassxx[iA] += dett;
				Mmassxx[iB] += dett;
				Mmassxx[iC] += dett;
			}

			//else: the triangle is a boundary triangle -> workT=-1
			else workT[i]=-1;
		}

		/*Check err, we need to make sure it has the right size!*/
		bool    deleteerr = false;
		double* err       = NULL;
		if(bamgopts->errSize[0]==1){
			/*Let's copy this value for all vertices*/
			err =  new double[nbv*bamgopts->errSize[1]];
			for(int i=0;i<nbv;i++){
				for(int j=0;j<bamgopts->errSize[1];j++){
					err[i*bamgopts->errSize[1]+j] = bamgopts->err[j];
				}
			}
			deleteerr = true;
		}
		else if(bamgopts->errSize[0]==this->nbv){
			/*Nothing to do, already right size*/
			err = bamgopts->err;
		}
		else{
			_error_("number of rows in 'err' not supported: size "<<bamgopts->errSize[0]<<"x"<<bamgopts->errSize[1]<<" (nbv is "<<this->nbv<<")");
		}

		//for all Solution  
		for (int nusol=0;nusol<nbsol;nusol++) {

			double smin=ss[nusol],smax=ss[nusol];
			double h1=1.e30,h2=1e-30,rx=0;
			double hn1=1.e30,hn2=1e-30,rnx =1.e-30;  

			//get min(s), max(s) and initialize Hessian (dxdx,dxdy,dydy)
			for ( iv=0,k=0; iv<nbv; iv++ ){
				dxdx[iv]=dxdy[iv]=dydy[iv]=0;
				smin=Min(smin,ss[iv*nbsol+nusol]);
				smax=Max(smax,ss[iv*nbsol+nusol]);
			}
			double sdelta=smax-smin;
			double absmax=Max(Abs(smin),Abs(smax));

			//display info
			if(verbose>2) _printf_("      Solution " << nusol << ", Min = " << smin << ", Max = " << smax << ", Delta = " << sdelta << ", number of fields = " << nbsol << "\n");

			//skip constant field
			if (sdelta < 1.0e-10*Max(absmax,1e-20) ){
				if (verbose>2) _printf_("      Solution " << nusol << " is constant, skipping...\n");
				continue;
			}

			//pointer toward ss that is also a pointer toward s (solutions)
			double* sf=ss; 

				//initialize the hessian matrix
				for ( iv=0,k=0; iv<nbv; iv++) dxdx[iv]=dxdy[iv]=dydy[iv]=0;

				//loop over the triangles
				for (i=0;i<nbt;i++){

					//for real all triangles 
					if(triangles[i].link){

						// coor of 3 vertices 
						R2 A=triangles[i][0];
						R2 B=triangles[i][1];
						R2 C=triangles[i][2];

						//warning: the normal is internal and the size is the length of the edge
						R2 nAB = Orthogonal(B-A);
						R2 nBC = Orthogonal(C-B);
						R2 nCA = Orthogonal(A-C);
						//note that :  nAB + nBC + nCA == 0 

						// number of the 3 vertices
						iA = GetId(triangles[i][0]);
						iB = GetId(triangles[i][1]);
						iC = GetId(triangles[i][2]);

						// for the test of  boundary edge
						// the 3 adj triangles 
						Triangle *tBC = triangles[i].TriangleAdj(OppositeEdge[0]);
						Triangle *tCA = triangles[i].TriangleAdj(OppositeEdge[1]);
						Triangle *tAB = triangles[i].TriangleAdj(OppositeEdge[2]);

						// value of the P1 fonction on 3 vertices 
						sA = ss[iA*nbsol+nusol];
						sB = ss[iB*nbsol+nusol];
						sC = ss[iC*nbsol+nusol];

						/*The nodal functions are such that for a vertex A:
						  N_A(x,y)=alphaA x + beta_A y +gamma_A
						  N_A(A) = 1,   N_A(B) = 0,   N_A(C) = 0
						  solving this system of equation (determinant = 2Area(T) != 0 if A,B and C are not inlined)
						  leads to:
						  N_A = (xB yC - xC yB + x(yB-yC) +y(xC-xB))/(2*Area(T))
						  and this gives:
						  alpha_A = (yB-yC)/(2*Area(T))
						  beta_A = (xC-xB)/(2*Area(T))
						  and therefore:
						  grad N_A = nA / detT
						  for an interpolation of a solution s:
						  grad(s) = s * sum_{i=A,B,C} grad(N_i) */

						R2 Grads=(nAB*sC+nBC*sA+nCA*sB)/detT[i];

						//Use Green to compute Hessian Matrix

						// if edge on boundary no contribution  => normal = 0
						if ( !tBC || !tBC->link ) nBC=O;
						if ( !tCA || !tCA->link ) nCA=O;
						if ( !tAB || !tAB->link ) nAB=O;

						// remark we forgot a 1/2 because
						//       int_{edge} w_i = 1/2 if i is in edge 
						//                         0  if not
						// if we don't take the  boundary 
						dxdx[iA] += ( nCA.x + nAB.x ) *Grads.x;
						dxdx[iB] += ( nAB.x + nBC.x ) *Grads.x;
						dxdx[iC] += ( nBC.x + nCA.x ) *Grads.x;

						//warning optimization (1) the division by 2 is done on the metric construction
						dxdy[iA] += (( nCA.y + nAB.y ) *Grads.x + ( nCA.x + nAB.x ) *Grads.y) ;
						dxdy[iB] += (( nAB.y + nBC.y ) *Grads.x + ( nAB.x + nBC.x ) *Grads.y) ;
						dxdy[iC] += (( nBC.y + nCA.y ) *Grads.x + ( nBC.x + nCA.x ) *Grads.y) ; 

						dydy[iA] += ( nCA.y + nAB.y ) *Grads.y;
						dydy[iB] += ( nAB.y + nBC.y ) *Grads.y;
						dydy[iC] += ( nBC.y + nCA.y ) *Grads.y;

					} // for real all triangles 
				}

				long kk=0;
				for ( iv=0,k=0 ; iv<nbv; iv++){
					if(Mmassxx[iv]>0){
						dxdx[iv] /= 2*Mmassxx[iv];
						// warning optimization (1) on term dxdy[iv]*ci/2 
						dxdy[iv] /= 4*Mmassxx[iv];
						dydy[iv] /= 2*Mmassxx[iv];
						// Compute the matrix with abs(eigen value)
						Metric M(dxdx[iv], dxdy[iv], dydy[iv]);
						EigenMetric Vp(M);
						Vp.Abs();
						M = Vp;
						dxdx[iv] = M.a11;
						dxdy[iv] = M.a21;
						dydy[iv] = M.a22;
					}
					else kk++;
				}

				// correction of second derivative
				// by a laplacien
				double* dd;
				for (int xy = 0;xy<3;xy++) {
					if      (xy==0) dd=dxdx;
					else if (xy==1) dd=dxdy;
					else if (xy==2) dd=dydy;
					else{
						delete [] detT;
						delete [] Mmass;
						delete [] workT;
						delete [] workV;
						delete [] Mmassxx;
						delete [] OnBoundary;
						_error_("not supported yet");
					}
					// do leat 2 iteration for boundary problem
					for (int ijacobi=0;ijacobi<Max(NbJacobi,2);ijacobi++){
						for (i=0;i<nbt;i++) 
						 if(triangles[i].link){// the real triangles 
							 // number of the 3 vertices
							 iA = GetId(triangles[i][0]);
							 iB = GetId(triangles[i][1]);
							 iC = GetId(triangles[i][2]);
							 double cc=3;
							 if(ijacobi==0)
							  cc = Max((double) ((Mmassxx[iA]>0)+(Mmassxx[iB]>0)+(Mmassxx[iC]>0)),1.);
							 workT[i] = (dd[iA]+dd[iB]+dd[iC])/cc;
						 }
						for (iv=0;iv<nbv;iv++) workV[iv]=0;

						for (i=0;i<nbt;i++){ 
							if(triangles[i].link){ // the real triangles 
								// number of the 3 vertices
								iA = GetId(triangles[i][0]);
								iB = GetId(triangles[i][1]);
								iC = GetId(triangles[i][2]);
								double cc =  workT[i]*detT[i];
								workV[iA] += cc;
								workV[iB] += cc;
								workV[iC] += cc;
							}
						}

						for (iv=0;iv<nbv;iv++){
							if( ijacobi<NbJacobi || OnBoundary[iv]){
								dd[iv] = workV[iv]/(Mmass[iv]*6);
							}
						}
					}
				}

				/*Compute Metric from Hessian*/
				for ( iv=0;iv<nbv;iv++){
					vertices[iv].MetricFromHessian(dxdx[iv],dxdy[iv],dydy[iv],smin,smax,ss[iv*nbsol+nusol],err[iv*nbsol+nusol],bamgopts);
				}

		}// end for all solution 

		delete [] detT;
		delete [] Mmass;
		delete [] dxdx;
		delete [] dxdy;
		delete [] dydy;
		delete [] workT;
		delete [] workV;
		delete [] Mmassxx;
		delete [] OnBoundary;
		if(deleteerr) delete [] err;

	}
	/*}}}*/
	void Mesh::CrackMesh(BamgOpts* bamgopts) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/CrackMesh)*/

		/*Intermediary*/
		int i,j,k,num,count;
		int i1,i2;
		int j1,j2;
		int verbose=0;

		/*Options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//  computed the number of cracked edge
		for (k=i=0;i<nbe;i++){
			if(edges[i].GeomEdgeHook->Cracked()) k++;
		}

		//Return if no edge is cracked
		if(k==0) return;
		if (verbose>4) _printf_("      number of Cracked Edges = " << k << "\n");

		//Initialize Cracked edge
		NbCrackedEdges=k;
		CrackedEdges=new CrackedEdge[k];

		//Compute number of Cracked Vertices
		k=0;
		NbCrackedVertices=0;

		int* splitvertex=new int[nbv];
		for (i=0;i<nbv;i++) splitvertex[i]=0;

		for (i=0;i<nbe;i++){
			if(edges[i].GeomEdgeHook->Cracked()){

				//Fill edges fields of CrackedEdges
				CrackedEdges[k  ].E =edges[i].GeomEdgeHook;
				CrackedEdges[k++].e1=&edges[i];

				//Get number of the two vertices on the edge
				i1=GetId(edges[i][0]);
				i2=GetId(edges[i][1]);
				_assert_(i1>=0 && i1<nbv && i2>=0 && i2<nbv);
				splitvertex[i1]++;
				splitvertex[i2]++;

				//If the vertex has already been flagged once, it is a cracked vertex (tip otherwise)
				if (splitvertex[i1]==2) NbCrackedVertices++;
				if (splitvertex[i2]==2) NbCrackedVertices++;

				//The vertex cannot be marked more than twice
				if (splitvertex[i1]==3 || splitvertex[i2]==3){
					delete [] splitvertex;
					_error_("Crossing rifts not supported yet");
				}
			}
		}
		_assert_(k==NbCrackedEdges);

		//Add new vertices
		if (verbose>4) _printf_("      number of Cracked Vertices = " << NbCrackedVertices << "\n");
		if (NbCrackedVertices){
			CrackedVertices=xNew<long>(2*NbCrackedVertices);
			num=0;
			for (i=0;i<nbv;i++){
				if (splitvertex[i]==2){
					CrackedVertices[num*2+0]=i;      //index of first vertex
					CrackedVertices[num*2+1]=nbv+num;//index of new vertex
					num++;
				}
			}
			_assert_(num==NbCrackedVertices);
		}
		delete [] splitvertex;

		//Now, find the triangles that hold a cracked edge
		CreateSingleVertexToTriangleConnectivity();

		long* Edgeflags=new long[NbCrackedEdges];
		for(i=0;i<NbCrackedEdges;i++) Edgeflags[i]=0;

		for(i=0;i<NbCrackedEdges;i++){
			//Get the numbers of the 2 vertices of the crren cracked edge
			i1=GetId((*CrackedEdges[i].e1)[0]);
			i2=GetId((*CrackedEdges[i].e1)[1]);

			//find a triangle holding the vertex i1 (first vertex of the ith cracked edge)
			Triangle* tbegin=vertices[i1].t;
			k=vertices[i1].IndexInTriangle;//local number of i in triangle tbegin
			_assert_(GetId((*tbegin)[k])==GetId(vertices[i1]));

			//Now, we are going to go through the adjacent triangle that hold i1 till
			//we find one that has the cracked edge
			AdjacentTriangle ta(tbegin,EdgesVertexTriangle[k][0]);
			count=0;
			do {
				for(j=0;j<3;j++){
					//Find the position of i1 in the triangle index
					if (GetId((*ta.t)[j])==i1){
						j1=j;
						break;
					}
				}
				for(j=0;j<3;j++){
					//Check wether i2 is also in the triangle index
					if (GetId((*ta.t)[j])==i2){
						j2=j;
						//Invert j1 and j2 if necessary
						if ((j1+1)%3==j2){
							int j3=j1;
							j1=j2;
							j2=j3;
						}
						if (Edgeflags[i]==0){
							//first element
							CrackedEdges[i].a=ta.t;
							CrackedEdges[i].length=Norme2((*ta.t)[j1].r-(*ta.t)[j2].r);
							CrackedEdges[i].normal=Orthogonal((*ta.t)[j1].r-(*ta.t)[j2].r);
						}
						else{
							//Second element -> to renumber
							CrackedEdges[i].b=ta.t;
							CrackedEdges[i].length=Norme2((*ta.t)[j1].r-(*ta.t)[j2].r);
							CrackedEdges[i].normal=Orthogonal((*ta.t)[j1].r-(*ta.t)[j2].r);
						}
						Edgeflags[i]++;
						break;
					}
				}
				//_printf_(element_renu[GetId(ta.t)] << " -> " << GetId((*ta.t)[0])+1 << " " << GetId((*ta.t)[1])+1 << " " << GetId((*ta.t)[2])+1 << ", edge [" << i1 << "->" << j1 << " " << i2 << "->" << j2 << "]\n");
				ta = Next(ta).Adj(); 
				if (count++>50) _error_("Maximum number of iteration exceeded");
			}while ((tbegin != ta)); 
		}

		//Check EdgeFlag
		for(i=0;i<NbCrackedEdges;i++){
			if (Edgeflags[i]!=2){
				_error_("A problem occured: at least one crack edge (number " << i+1 << ") does not belong to 2 elements");
			}
		}
		delete [] Edgeflags;

		//Reset BamgVertex to On
		SetVertexFieldOn();

	}
	/*}}}*/
	void Mesh::Echo(void) {/*{{{*/

		int i;

		_printf_("Mesh Echo:\n");
		_printf_("   nbv = " << nbv << "\n");
		_printf_("   nbt = " << nbt << "\n");
		_printf_("   nbe = " << nbe << "\n");
		_printf_("   index:\n");
		for (i=0;i<nbt;i++){
			_printf_("   " << setw(4) << i+1 << ": [" 
						<< setw(4) << (((BamgVertex*)triangles[i](0))?GetId(triangles[i][0])+1:0) << " " 
						<< setw(4) << (((BamgVertex*)triangles[i](1))?GetId(triangles[i][1])+1:0) << " " 
						<< setw(4) << (((BamgVertex*)triangles[i](2))?GetId(triangles[i][2])+1:0) << "]\n");
		}
		_printf_("   coordinates:\n");
		for (i=0;i<nbv;i++){
			_printf_("   " << setw(4) << i+1 << ": [" << vertices[i].r.x << " " << vertices[i].r.y << "] and [" << vertices[i].i.x << " " << vertices[i].i.y << "]\n");
		}

	}
	/*}}}*/
	void Mesh::ForceBoundary(BamgOpts* bamgopts) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/ForceBoundary)*/

		int verbose=0;
		int k=0;
		int nbfe=0,nbswp=0,Nbswap=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//display
		if (verbose > 2) _printf_("   ForceBoundary  nb of edge: " << nbe << "\n");

		//check that there is no triangle with 0 determinant
		for (int t = 0; t < nbt; t++){
			if (!triangles[t].det) k++;
		}
		if (k!=0) {
			_error_("there is " << k << " triangles of mes = 0");
		}

		//Force Edges
		AdjacentTriangle ta(0,0);
		for (int i = 0; i < nbe; i++){

			//Force edge i
			nbswp =  ForceEdge(edges[i][0],edges[i][1],ta);
			if (nbswp<0) k++;
			else Nbswap += nbswp;

			if (nbswp) nbfe++;
			if ( nbswp < 0 && k < 5){
				_error_("Missing Edge " << i << ", v0=" << GetId(edges[i][0]) << ",v1=" << GetId(edges[i][1]));
			}
		}

		if (k!=0) {
			_error_("There are " << k << " lost edges, the boundary might be crossing");
		}
		for (int j=0;j<nbv;j++){
			Nbswap +=  vertices[j].Optim(1,0);
		}
		if (verbose > 3) _printf_("      number of inforced edge = " << nbfe << ", number of swap= " << Nbswap << "\n"); 
	}
	/*}}}*/
	void Mesh::FindSubDomain(BamgOpts* bamgopts,int OutSide) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/FindSubDomain)*/

		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		if (verbose >2){
			if (OutSide) _printf_("   Find all external sub-domain\n"); 
			else _printf_("   Find all internal sub-domain\n");
		}
		short * HeapArete = new short[nbt];
		Triangle  **  HeapTriangle = new Triangle*  [nbt];
		Triangle *t,*t1;
		long k,it;

		/*No color by default*/
		for(int itt=0;itt<nbt;itt++) triangles[itt].link=0;

		long  NbSubDomTot =0;
		for(it=0;it<nbt;it++)  { 
			if( !triangles[it].link ){
				t = triangles + it;
				NbSubDomTot++;; // new composante connexe
				long i = 0; // niveau de la pile 
				t->link = t ; // sd forme d'un triangle cicular link

				HeapTriangle[i] =t ; 
				HeapArete[i] = 3;

				while(i >= 0) // boucle sur la pile
				  { while (HeapArete[i]--) // boucle sur les 3 aretes 
					  { 
						int na =  HeapArete[i];
						Triangle * tc =  HeapTriangle[i]; // triangle courant
						if( ! tc->Locked(na)) // arete non frontiere
						  {
							Triangle * ta = tc->TriangleAdj(na) ; // næ triangle adjacent
							if (ta->link == 0 ) // non deja chainer => on enpile
							  { 
								i++;
								ta->link = t->link ;  // on chaine les triangles
								t->link = ta ;  // d'un meme sous domaine          
								HeapArete[i] = 3; // pour les 3 triangles adjacents
								HeapTriangle[i] = ta;
							  }}
					  } // deplie fin de boucle sur les 3 adjacences
					i--;
				  }          
			}      
		}

		// supression de tous les sous domaine infini <=>  contient le sommet NULL
		it =0;
		nbtout = 0;
		while (it<nbt) {
			if (triangles[it].link) 
			  { 
				if (!( triangles[it](0) &&  triangles[it](1) &&  triangles[it](2) )) 
				  {
					// infini triangle 
					NbSubDomTot --;
					t=&triangles[it];
					nbtout--;  // on fait un coup de trop. 
					while  (t){
						nbtout++;
						t1=t;
						t=t->link;
						t1->link=0;
					}
				  }
			  }   
			it++;} // end while (it<nbt)
			if (nbt == nbtout ||  !NbSubDomTot) {
				delete [] HeapArete;
				delete [] HeapTriangle;
				_error_("The boundary is not close: all triangles are outside");
			}

			delete [] HeapArete;
			delete [] HeapTriangle;

			if (OutSide|| !Gh.subdomains || !Gh.nbsubdomains ) 
			  { // No geom sub domain
				long i;
				if (subdomains) delete [] subdomains;
				subdomains = new SubDomain[ NbSubDomTot];
				nbsubdomains=  NbSubDomTot;
				for ( i=0;i<nbsubdomains;i++) {
					subdomains[i].head=NULL;
					subdomains[i].ReferenceNumber=i+1;
				}
				long * mark = new long[nbt];
				for (it=0;it<nbt;it++)
				 mark[it]=triangles[it].link ? -1 : -2;

				it =0;
				k = 0;
				while (it<nbt) {
					if (mark[it] == -1) {
						t1 = & triangles[it];
						t = t1->link;
						mark[it]=k;
						subdomains[k].head = t1;
						do {
							mark[GetId(t)]=k;
							t=t->link;
						} while (t!=t1);
						mark[it]=k++;}
						//    else if(mark[it] == -2 ) triangles[it].Draw(999);
						it++;} // end white (it<nbt)
						if (k!=nbsubdomains){
							delete [] mark;
							_error_("k!=nbsubdomains");
						}
						if(OutSide){
							//  to remove all the sub domain by parity adjacents
							//  because in this case we have only the true boundary edge
							//  so the boundary is manifold
							long nbk = nbsubdomains;
							while (nbk)
							 for (it=0;it<nbt && nbk ;it++)
							  for (int na=0;na<3 && nbk ;na++)
								 {
								  Triangle *ta = triangles[it].TriangleAdj(na);
								  long kl = ta ? mark[GetId(ta)] : -2;
								  long kr = mark[it];
								  if(kr !=kl) {
									  if (kl >=0 && subdomains[kl].ReferenceNumber <0 && kr >=0 && subdomains[kr].ReferenceNumber>=0)
										nbk--,subdomains[kr].ReferenceNumber=subdomains[kl].ReferenceNumber-1;
									  if (kr >=0 && subdomains[kr].ReferenceNumber <0 && kl >=0 && subdomains[kl].ReferenceNumber>=0)
										nbk--,subdomains[kl].ReferenceNumber=subdomains[kr].ReferenceNumber-1;
									  if(kr<0 && kl >=0 && subdomains[kl].ReferenceNumber>=0)
										nbk--,subdomains[kl].ReferenceNumber=-1;
									  if(kl<0 && kr >=0 && subdomains[kr].ReferenceNumber>=0)
										nbk--,subdomains[kr].ReferenceNumber=-1;
								  }
								 }
							long  j=0;
							for ( i=0;i<nbsubdomains;i++)
							 if((-subdomains[i].ReferenceNumber) %2) { // good 
								 if(i != j) 
								  Exchange(subdomains[i],subdomains[j]);
								 j++;}
							 else{ 
								 t= subdomains[i].head;
								 while (t){
									 nbtout++;
									 t1=t;
									 t=t->link;
									 t1->link=0;
								 }//while (t)
								}
							if(verbose>4) _printf_("      Number of removes subdomains (OutSideMesh) = " << nbsubdomains-j << "\n");
							nbsubdomains=j;
						  }

						delete []  mark; 

			  }
			else{ // find the head for all subdomains
				if (Gh.nbsubdomains != nbsubdomains && subdomains)
				 delete [] subdomains, subdomains=0;
				if (! subdomains  ) 
				 subdomains = new SubDomain[ Gh.nbsubdomains];
				nbsubdomains =Gh.nbsubdomains;
				CreateSingleVertexToTriangleConnectivity();
				long * mark = new long[nbt];
				Edge **GeomEdgetoEdge = MakeGeomEdgeToEdge();

				for (it=0;it<nbt;it++)
				 mark[it]=triangles[it].link ? -1 : -2;
				long inew =0;
				for (int i=0;i<nbsubdomains;i++) {
					GeomEdge &eg = *Gh.subdomains[i].edge;
					subdomains[i].ReferenceNumber = Gh.subdomains[i].ReferenceNumber;
					// by carefull is not easy to find a edge create from a GeomEdge 
					// see routine MakeGeomEdgeToEdge
					Edge &e = *GeomEdgetoEdge[Gh.GetId(eg)];
					_assert_(&e);
					BamgVertex * v0 =  e(0),*v1 = e(1);
					Triangle *t  = v0->t;
					int direction = Gh.subdomains[i].direction;
					// test if ge and e is in the same direction 
					if (((eg[0].r-eg[1].r),(e[0].r-e[1].r))<0) direction = -direction ;
					subdomains[i].direction = direction;
					subdomains[i].edge = &e;
					_assert_(t && direction);

					AdjacentTriangle  ta(t,EdgesVertexTriangle[v0->IndexInTriangle][0]);// previous edges

					while (1) {
						_assert_(v0==ta.EdgeVertex(1));
						if (ta.EdgeVertex(0) == v1) { // ok we find the edge
							if (direction>0)  
							 subdomains[i].head=t=Adj(ta);
							else 
							 subdomains[i].head=t=ta;
							if(t<triangles || t >= triangles+nbt || t->det < 0 || t->link == 0) {
								_error_("bad definition of SubSomain " << i);
							}
							long it = GetId(t);
							if (mark[it] >=0) {
								break;
							}
							if(i != inew) 
							 Exchange(subdomains[i],subdomains[inew]);
							inew++;
							Triangle *tt=t;
							long kkk=0;
							do 
							  {
								kkk++;
								if (mark[GetId(tt)]>=0){
									_error_("mark[GetId(tt)]>=0");
								}
								mark[GetId(tt)]=i;
								tt=tt->link;
							  } while (tt!=t);
							break;
						}
						ta = Previous(Adj(ta));         
						if(t == (Triangle *) ta) {
							_error_("bad definition of SubSomain " << i);
						}
					}
				}

				if (inew < nbsubdomains) {
					if (verbose>5) _printf_("WARNING: " << nbsubdomains-inew << " SubDomains are being removed\n");
					nbsubdomains=inew;}

					for (it=0;it<nbt;it++)
					 if ( mark[it] ==-1 ) 
					  nbtout++,triangles[it].link =0;
					delete [] GeomEdgetoEdge;
					delete [] mark;

			  }
			nbtout=0;
			for (it=0;it<nbt;it++) 
			 if(!triangles[it].link)  nbtout++;
	}
	/*}}}*/
	long Mesh::GetId(const Triangle & t) const  { /*{{{*/
		return &t - triangles;
	}
	/*}}}*/
	long Mesh::GetId(const Triangle * t) const  { /*{{{*/
		return t - triangles;
	}
	/*}}}*/
	long Mesh::GetId(const BamgVertex & t) const  { /*{{{*/
		return &t - vertices;
	}
	/*}}}*/
	long Mesh::GetId(const BamgVertex * t) const  { /*{{{*/
		return t - vertices;
	}
	/*}}}*/
	long Mesh::GetId(const Edge & t) const  { /*{{{*/
		return &t - edges;
	}
	/*}}}*/
	long Mesh::GetId(const Edge * t) const  { /*{{{*/
		return t - edges;
	}
	/*}}}*/
	void Mesh::Init(long maxnbv_in) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/PreInit)*/

		/*Initialize fields*/
		this->NbRef                  = 0;
		this->quadtree               = NULL;
		this->nbv                    = 0;
		this->nbt                    = 0;
		this->nbe                    = 0;
		this->edges                  = NULL;
		this->nbsubdomains           = 0;
		this->subdomains             = NULL;
		this->maxnbv                 = maxnbv_in;
		this->maxnbt                 = 2 *maxnbv_in-2;
		this->NbVertexOnBThVertex    = 0;
		this->VertexOnBThVertex      = NULL;
		this->NbVertexOnBThEdge      = 0;
		this->VertexOnBThEdge        = NULL;
		this->NbCrackedVertices      = 0;
		this->CrackedVertices        = NULL;
		this->NbCrackedEdges         = 0;
		this->CrackedEdges           = NULL;
		this->NbVerticesOnGeomVertex = 0;
		this->VerticesOnGeomVertex   = NULL;
		this->NbVerticesOnGeomEdge   = 0;
		this->VerticesOnGeomEdge     = NULL;

		/*Initialize random seed*/
		this->randomseed = 1;

		/*Allocate if maxnbv_in>0*/
		if(maxnbv_in){
			this->vertices=new BamgVertex[this->maxnbv];
			this->orderedvertices=new BamgVertex* [this->maxnbv];
			this->triangles=new Triangle[this->maxnbt];
         _assert_(this->vertices);
         _assert_(this->orderedvertices);
			_assert_(this->triangles);
		}
		else{
			this->vertices        = NULL;
			this->orderedvertices = NULL;
			this->triangles       = NULL;
			this->maxnbt          = 0;
		} 
	}
	/*}}}*/
	void Mesh::Insert(BamgOpts* bamgopts){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Insert)*/

		/*Insert points in the existing Geometry*/

		//Intermediary
		long i;
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//Display info
		if (verbose>2) _printf_("   Insert initial " << nbv << " vertices\n");

		//Compute integer coordinates for the existing vertices
		SetIntCoor();

		/*Now we want to build a list (orderedvertices) of the vertices in a random
		 * order. To do so, we use the following method:
		 *
		 * From an initial k0 in [0 nbv[ random (vertex number)
		 * the next k (vertex number) is computed using a big
		 * prime number (PN>>nbv) following:
		 *
		 * k_{i+1} = k_i + PN  [nbv]
		 *
		 * let's show that:
		 *
		 *   for all j in [0 nbv[, there is a unique i in [0 nbv[ such that k_i=j
		 *
		 * Let's assume that there are 2 distinct j1 and j2 such that
		 * k_j1 = k_j2
		 *
		 * This means that
		 *  
		 *  k0+j1*PN = k0+j2*PN [nbv]
		 *  (j1-j2)*PN =0       [nbv]
		 * since PN is a prime number larger than nbv, and nbv!=1
		 *  j1-j2=0             [nbv]
		 * BUT
		 *  j1 and j2 are in [0 nbv[ which is impossible.
		 *
		 *  We hence have built a random list of nbv elements of
		 *  [0 nbv[ all distincts*/

		//Get Prime number
		const long PrimeNumber= BigPrimeNumber(nbv);
		long k0=this->RandomNumber(nbv);
		if (verbose>4) _printf_("      Prime Number = "<<PrimeNumber<<"\n");
		if (verbose>4) _printf_("      k0 = "<<k0<<"\n");

		//Build orderedvertices
		for (i=0; i<nbv; i++){
			orderedvertices[i]=&vertices[k0=(k0+PrimeNumber)%nbv];
		}

		/*Modify orderedvertices such that the first 3 vertices form a triangle*/

		//get first vertex i such that [0,1,i] are not aligned
		for (i=2; det(orderedvertices[0]->i,orderedvertices[1]->i,orderedvertices[i]->i)==0;){
			//if i is higher than nbv, it means that all the determinants are 0,
			//all vertices are aligned!
			if  (++i>=nbv) _error_("all the vertices are aligned");
		}
		if (verbose>4) _printf_("      i = "<<i<<"\n");
		// exchange i et 2 in "orderedvertices" so that
		// the first 3 vertices are not aligned (real triangle)
		Exchange(orderedvertices[2], orderedvertices[i]);

		/*Take the first edge formed by the first two vertices and build
		 * the initial simple mesh from this edge and 2 boundary triangles*/

		BamgVertex *v0=orderedvertices[0], *v1=orderedvertices[1];

		nbt = 2;
		triangles[0](0) = NULL;//infinite vertex
		triangles[0](1) = v0;
		triangles[0](2) = v1;
		triangles[1](0) = NULL;//infinite vertex
		triangles[1](2) = v0;
		triangles[1](1) = v1;

		//Build adjacence
		const int e0 = OppositeEdge[0];
		const int e1 = NextEdge[e0];
		const int e2 = PreviousEdge[e0];
		triangles[0].SetAdj2(e0, &triangles[1] ,e0);
		triangles[0].SetAdj2(e1, &triangles[1] ,e2);
		triangles[0].SetAdj2(e2, &triangles[1] ,e1);

		triangles[0].det = -1;  //boundary triangle: det = -1
		triangles[1].det = -1;  //boundary triangle: det = -1

		triangles[0].SetSingleVertexToTriangleConnectivity();
		triangles[1].SetSingleVertexToTriangleConnectivity();

		triangles[0].link=&triangles[1];
		triangles[1].link=&triangles[0];

		//build quadtree
		if (!quadtree)  quadtree = new BamgQuadtree(this,0);
		quadtree->Add(*v0);
		quadtree->Add(*v1);

		/*Now, add the vertices One by One*/
		long NbSwap=0;
		if (verbose>3) _printf_("   Begining of insertion process...\n");
		if (verbose>4) _printf_("      nbv = "<<nbv<<"\n");

		for (long icount=2; icount<nbv; icount++) {

			//Get new vertex
			BamgVertex *newvertex=orderedvertices[icount];

			//Find the triangle in which newvertex is located
			long long det3[3];
			Triangle* tcvi = TriangleFindFromCoord(newvertex->i,det3); //(newvertex->i = integer coordinates)

			//Add newvertex to the quadtree
			quadtree->Add(*newvertex); 

			//Add newvertex to the existing mesh
			AddVertex(*newvertex,tcvi,det3);

			//Make the mesh Delaunay around newvertex by swaping the edges
			NbSwap += newvertex->Optim(1,0);
		}

		//Display info
		if (verbose>3) {
			_printf_("      NbSwap of insertion: " << NbSwap << "\n");
			_printf_("      NbSwap/nbv:          " << NbSwap/nbv << "\n");
		}
	}
	/*}}}*/
	long Mesh::InsertNewPoints(long nbvold,long & NbTSwap,BamgOpts* bamgopts) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/InsertNewPoints)*/

		int verbose=0;
		double seuil= 1.414/2.;// for two close point 
		long i;
		long NbSwap=0;
		long long det3[3];

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//number of new points
		const long nbvnew=nbv-nbvold;

		//display info if required
		if (verbose>5) _printf_("      Try to Insert " << nbvnew << " new points\n");

		//return if no new points
		if (!nbvnew) return 0; 

		/*construction of a random order*/
		const long PrimeNumber= BigPrimeNumber(nbv)  ;
		long k3 = this->RandomNumber(nbvnew);
		//loop over the new points
		for (int is3=0; is3<nbvnew; is3++){
			long j=nbvold +(k3 = (k3+PrimeNumber)%nbvnew);
			long i=nbvold+is3; 
			orderedvertices[i]= vertices + j;
			orderedvertices[i]->ReferenceNumber=i;
		}

		// for all the new point
		long iv=nbvold;
		for(i=nbvold;i<nbv;i++){
			BamgVertex &vi=*orderedvertices[i];
			vi.i=R2ToI2(vi.r);
			vi.r=I2ToR2(vi.i);
			double hx,hy;
			vi.m.Box(hx,hy);
			int hi=(int) (hx*coefIcoor),hj=(int) (hy*coefIcoor);
			if(!quadtree->TooClose(&vi,seuil,hi,hj)){
				// a good new point 
				BamgVertex &vj = vertices[iv];
				long  j=vj.ReferenceNumber; 
				if (&vj!=orderedvertices[j]){
					_error_("&vj!= orderedvertices[j]");
				}
				if(i!=j){ 
					Exchange(vi,vj);
					Exchange(orderedvertices[j],orderedvertices[i]);
				}
				vj.ReferenceNumber=0; 
				Triangle *tcvj=TriangleFindFromCoord(vj.i,det3);
				if (tcvj && !tcvj->link){
					_printf_("While trying to add the following point:\n");
					vj.Echo();
					_printf_("BAMG determined that it was inside the following triangle, which probably lies outside of the geometric domain\n");
					tcvj->Echo();
					_error_("problem inserting point in InsertNewPoints (tcvj=" << tcvj << " and tcvj->link=" << tcvj->link << ")");
				}
				quadtree->Add(vj);
				AddVertex(vj,tcvj,det3);
				NbSwap += vj.Optim(1);          
				iv++;
			}
			else{
				vi.PreviousNumber = 0;
			}
		} 
		if (verbose>3) {
			_printf_("         number of new points: " << iv << "\n");
			_printf_("         number of to close (?) points: " << nbv-iv << "\n");
			_printf_("         number of swap after: " << NbSwap << "\n");
		}
		nbv = iv;

		for (i=nbvold;i<nbv;i++) NbSwap += vertices[i].Optim(1);  
		if (verbose>3) _printf_("   NbSwap=" << NbSwap << "\n");

		NbTSwap +=  NbSwap ;
		return nbv-nbvold;
	}
	/*}}}*/
	Edge** Mesh::MakeGeomEdgeToEdge() {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/MakeGeomEdgeToEdge)*/

		if (!Gh.nbe){
			_error_("!Gh.nbe");
		}
		Edge **e= new Edge* [Gh.nbe];

		long i;
		for ( i=0;i<Gh.nbe ; i++)
		 e[i]=NULL;
		for ( i=0;i<nbe ; i++) 
		  { 
			Edge * ei = edges+i;
			GeomEdge *GeomEdgeHook = ei->GeomEdgeHook; 
			e[Gh.GetId(GeomEdgeHook)] = ei;    
		  }
		for ( i=0;i<nbe ; i++) 
		 for (int ii=0;ii<2;ii++) { 
			 Edge * ei = edges+i;
			 GeomEdge *GeomEdgeHook = ei->GeomEdgeHook;
			 int j= ii;
			 while (!(*GeomEdgeHook)[j].Required()) { 
				 Adj(GeomEdgeHook,j); // next geom edge
				 j=1-j;
				 if (e[Gh.GetId(GeomEdgeHook)])  break; // optimisation
				 e[Gh.GetId(GeomEdgeHook)] = ei; 
			 }
		 }

		int kk=0;
		for ( i=0;i<Gh.nbe ; i++){
			if (!e[i]){
				kk++;
				if(kk<10) _printf_("BUG: the geometrical edge " << i << " is on no edge curve\n");
			}
		}
		if(kk){
			delete [] e;
			_error_("See above");
		}

		return e;
	}
	/*}}}*/
	void Mesh::MakeBamgQuadtree() {  /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/MakeBamgQuadtree)*/
		if(!quadtree) quadtree = new BamgQuadtree(this);
	}
	/*}}}*/
	double Mesh::MaximalHmax() {/*{{{*/
		return Max(pmax.x-pmin.x,pmax.y-pmin.y);
	}
	/*}}}*/
	void  Mesh::MaxSubDivision(BamgOpts* bamgopts,double maxsubdiv) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/MaxSubDivision)*/

		/*Intermediaries*/
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		const  double maxsubdiv2 = maxsubdiv*maxsubdiv;
		if(verbose>1) _printf_("   Limit the subdivision of a edges in the new mesh by " << maxsubdiv << "\n");
		// for all the edges 
		// if the len of the edge is to long 
		long it,nbchange=0;    
		double lmax=0;
		for (it=0;it<nbt;it++){
			Triangle &t=triangles[it];
			for (int j=0;j<3;j++){
				Triangle *ptt=t.TriangleAdj(j);
				Triangle &tt = *ptt;
				if ( (!ptt ||  it < GetId(tt)) && ( tt.link || t.link)){
					BamgVertex &v0 = t[VerticesOfTriangularEdge[j][0]];
					BamgVertex &v1 = t[VerticesOfTriangularEdge[j][1]];
					R2 AB= (R2) v1-(R2) v0;
					Metric M = v0;
					double l = M(AB,AB);
					lmax = Max(lmax,l);
					if(l> maxsubdiv2){
						R2 AC = M.Orthogonal(AB);// the ortogonal vector of AB in M
						double lc = M(AC,AC);
						D2xD2 Rt(AB,AC);// Rt.x = AB , Rt.y = AC;
						D2xD2 Rt1(Rt.inv());
						D2xD2 D(maxsubdiv2,0,0,lc);
						D2xD2 MM = Rt1*D*Rt1.t();
						v0.m =  M = Metric(MM.x.x,MM.y.x,MM.y.y);
						nbchange++;
					}
					M = v1;
					l = M(AB,AB);
					lmax = Max(lmax,l);
					if(l> maxsubdiv2){
						R2 AC = M.Orthogonal(AB);// the ortogonal vector of AB in M
						double lc = M(AC,AC);
						D2xD2 Rt(AB,AC);// Rt.x = AB , Rt.y = AC;
						D2xD2 Rt1(Rt.inv());
						D2xD2 D(maxsubdiv2,0,0,lc);
						D2xD2  MM = Rt1*D*Rt1.t();
						v1.m =  M = Metric(MM.x.x,MM.y.x,MM.y.y);
						nbchange++;
					}
				}
			}
		}
		if(verbose>3){
			_printf_("      number of metric changes = " << nbchange << ", maximum number of subdivision of a edges before change = " << pow(lmax,0.5) << "\n");
		}
	}
	/*}}}*/
	Metric Mesh::MetricAt(const R2 & A){ /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/MetricAt)*/

		I2 a = R2ToI2(A);
		long long deta[3];
		Triangle * t =TriangleFindFromCoord(a,deta);
		if (t->det <0) { // outside
			double ba,bb;
			AdjacentTriangle edge= CloseBoundaryEdge(a,t,ba,bb) ;
			return Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));}
		else { // inside
			double   aa[3];
			double s = deta[0]+deta[1]+deta[2];
			aa[0]=deta[0]/s;
			aa[1]=deta[1]/s;
			aa[2]=deta[2]/s;
			return Metric(aa,(*t)[0],(*t)[1],(*t)[2]);
		}
	}
	/*}}}*/
	double Mesh::MinimalHmin() {/*{{{*/
		return 2.0/coefIcoor;
	}
	/*}}}*/
	BamgVertex* Mesh::NearestVertex(int i,int j) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/NearestVertex)*/
		return  quadtree->NearestVertex(i,j); 
	} 
	/*}}}*/
	void  Mesh::NewPoints(Mesh & Bh,BamgOpts* bamgopts,int KeepVertices){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/NewPoints)*/

		int i,j,k;
		int verbose=0;
		long NbTSwap=0;
		long nbtold=nbt;
		long nbvold=nbv;
		long Headt=0;
		long next_t;
		long* first_np_or_next_t=new long[maxnbt];
		Triangle* t=NULL;

		/*Recover options*/
		if(bamgopts) verbose=bamgopts->verbose;

		/*First, insert old points if requested*/
		if(KeepVertices && (&Bh != this) && (nbv+Bh.nbv< maxnbv)){
			if (verbose>5) _printf_("         Inserting initial mesh points\n");
			bool pointsoutside = false;
			for(i=0;i<Bh.nbv;i++){ 
				BamgVertex &bv=Bh[i];
				/*Do not insert if the point is outside*/
				long long det3[3];
				Triangle* tcvj=TriangleFindFromCoord(bv.i,det3);
				if(tcvj->det<0 || !tcvj->link){
					pointsoutside = true;
					continue;
				}
				IssmDouble area_1=((bv.r.x -(*tcvj)(2)->r.x)*((*tcvj)(1)->r.y-(*tcvj)(2)->r.y) 
						- (bv.r.y -(*tcvj)(2)->r.y)*((*tcvj)(1)->r.x-(*tcvj)(2)->r.x));
				IssmDouble area_2=(((*tcvj)(0)->r.x -(*tcvj)(2)->r.x)*(bv.r.y -(*tcvj)(2)->r.y) 
						- ((*tcvj)(0)->r.y -(*tcvj)(2)->r.y)*(bv.r.x -(*tcvj)(2)->r.x));
				IssmDouble area_3 =((bv.r.x -(*tcvj)(1)->r.x)*((*tcvj)(0)->r.y-(*tcvj)(1)->r.y)
						- (bv.r.y -(*tcvj)(1)->r.y)*((*tcvj)(0)->r.x-(*tcvj)(1)->r.x));
				if(area_1<0 || area_2<0 || area_3<0){
					pointsoutside = true;
					continue;
				}
				if(!bv.GeomEdgeHook){
					vertices[nbv].r              = bv.r;
					vertices[nbv].PreviousNumber = i+1;
					vertices[nbv++].m = bv.m;
				}
			}
			//if(pointsoutside) _printf_("WARNING: One or more points of the initial mesh fall outside of the geometric boundary\n");
			Bh.CreateSingleVertexToTriangleConnectivity();     
			InsertNewPoints(nbvold,NbTSwap,bamgopts);
		}
		else Bh.CreateSingleVertexToTriangleConnectivity();     

		// generation of the list of next Triangle 
		for(i=0;i<nbt;i++) first_np_or_next_t[i]=-(i+1);
		// the next traingle of i is -first_np_or_next_t[i]

		// Big loop (most time consuming)
		int iter=0;
		if (verbose>5) _printf_("         Big loop\n");
		do {
			/*Update variables*/
			iter++;
			nbtold=nbt;
			nbvold=nbv;

			/*We test all triangles*/
			i=Headt;
			next_t=-first_np_or_next_t[i];
			for(t=&triangles[i];i<nbt;t=&triangles[i=next_t],next_t=-first_np_or_next_t[i]){

				//check i
				if (i<0 || i>=nbt ){
					_error_("Index problem in NewPoints (i=" << i << " not in [0 " << nbt-1 << "])");
				}
				//change first_np_or_next_t[i]
				first_np_or_next_t[i] = iter; 

				//Loop over the edges of t
				for(j=0;j<3;j++){
					AdjacentTriangle tj(t,j);
					BamgVertex &vA = *tj.EdgeVertex(0);
					BamgVertex &vB = *tj.EdgeVertex(1);

					//if t is a boundary triangle, or tj locked, continue
					if (!t->link)     continue;
					if (t->det <0)    continue;
					if (t->Locked(j)) continue;

					AdjacentTriangle tadjj = t->Adj(j);	  
					Triangle* ta=tadjj;

					//if the adjacent triangle is a boundary triangle, continue
					if (ta->det<0) continue;	  

					R2 A=vA;
					R2 B=vB;
					k=GetId(ta);

					//if this edge has already been done, go to next edge of triangle
					if(first_np_or_next_t[k]==iter) continue;

					lIntTria.SplitEdge(Bh,A,B);
					lIntTria.NewPoints(vertices,nbv,maxnbv);
				} // end loop for each edge 
			}// for triangle   

			if (!InsertNewPoints(nbvold,NbTSwap,bamgopts)) break;
			for (i=nbtold;i<nbt;i++) first_np_or_next_t[i]=iter;
			Headt = nbt; // empty list 

			// for all the triangle containing the vertex i
			for (i=nbvold;i<nbv;i++){ 
				BamgVertex*          s  = vertices + i;
				AdjacentTriangle ta(s->t, EdgesVertexTriangle[s->IndexInTriangle][1]);
				Triangle*        tbegin= (Triangle*) ta;
				long kt;
				do { 
					kt = GetId((Triangle*) ta);
					if (first_np_or_next_t[kt]>0){
						first_np_or_next_t[kt]=-Headt;
						Headt=kt;
					}
					if (ta.EdgeVertex(0)!=s){
						_error_("ta.EdgeVertex(0)!=s");
					}
					ta = Next(Adj(ta));
				} while ( (tbegin != (Triangle*) ta)); 
			}

		}while(nbv!=nbvold);
		delete [] first_np_or_next_t;

		long NbSwapf =0;
		for(i=0;i<nbv;i++) NbSwapf += vertices[i].Optim(0);
	}/*}}}*/
	GeomEdge*   Mesh::ProjectOnCurve( Edge & BhAB, BamgVertex &  vA, BamgVertex & vB,/*{{{*/
				double theta,BamgVertex & R,VertexOnEdge &  BR,VertexOnGeom & GR) {
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshQuad.cpp/ProjectOnCurve)*/

		void *pA=0,*pB=0;
		double tA=0,tB=0;
		R2 A=vA,B=vB;
		BamgVertex * pvA=&vA, * pvB=&vB;
		if (vA.IndexInTriangle == IsVertexOnVertex){
			pA=vA.BackgroundVertexHook;
		}
		else if (vA.IndexInTriangle == IsVertexOnEdge){
			pA=vA.BackgroundEdgeHook->be;
			tA=vA.BackgroundEdgeHook->abcisse;
		}
		else {
			_error_("ProjectOnCurve On BamgVertex " << BTh.GetId(vA) << " forget call to SetVertexFieldOnBTh");
		} 

		if (vB.IndexInTriangle == IsVertexOnVertex){
			pB=vB.BackgroundVertexHook;
		}
		else if(vB.IndexInTriangle == IsVertexOnEdge){
			pB=vB.BackgroundEdgeHook->be;
			tB=vB.BackgroundEdgeHook->abcisse;
		}
		else {
			_error_("ProjectOnCurve On BamgVertex " << BTh.GetId(vB) << " forget call to SetVertexFieldOnBTh");
		} 
		Edge * e = &BhAB;
		if (!pA || !pB || !e){
			_error_("!pA || !pB || !e");
		}
		// be carefull the back ground edge e is on same geom edge 
		// of the initiale edge def by the 2 vertex A B;
		//check Is a background Mesh;   
		if (e<BTh.edges || e>=BTh.edges+BTh.nbe){
			_error_("e<BTh.edges || e>=BTh.edges+BTh.nbe");
		}
		// walk on BTh edge 
		//not finish ProjectOnCurve with BackGround Mesh);
		// 1 first find a back ground edge contening the vertex A
		// 2 walk n back gound boundary to find the final vertex B

		if( vA.IndexInTriangle == IsVertexOnEdge) 
		  { // find the start edge 
			e = vA.BackgroundEdgeHook->be;	 

		  } 
		else if (vB.IndexInTriangle == IsVertexOnEdge) 
		  {
			theta = 1-theta;
			Exchange(tA,tB);
			Exchange(pA,pB);
			Exchange(pvA,pvB);
			Exchange(A,B);
			e =  vB.BackgroundEdgeHook->be;

		  } 
		else{ // do the search by walking 
			_error_("case not supported yet");
		}

		// find the direction of walking with direction of edge and pA,PB;
		R2 AB=B-A;

		double cosE01AB = (( (R2) (*e)[1] - (R2) (*e)[0] ) , AB);
		int kkk=0;
		int direction = (cosE01AB>0) ? 1 : 0;

		//   double l=0; // length of the edge AB
		double abscisse = -1;

		for (int step=0;step<2;step++){
			// 2 times algo:
			//    1 for computing the length l
			//    2 for find the vertex 
			int  iii;
			BamgVertex  *v0=pvA,*v1; 
			Edge *neee,*eee;
			double lg =0; // length of the curve 
			double te0;
			// we suppose take the curve's abcisse 
			for ( eee=e,iii=direction,te0=tA;
						eee && ((( void*) eee) != pB) && (( void*) (v1=&((*eee)[iii]))) != pB ;
						neee = eee->adj[iii],iii = 1-neee->Intersection(*eee),eee = neee,v0=v1,te0=1-iii ) { 

				kkk=kkk+1;
				_assert_(kkk<100);
				_assert_(eee);
				double lg0 = lg;
				double dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
				lg += dp;
				if (step && abscisse <= lg) { // ok we find the geom edge 
					double sss  =   (abscisse-lg0)/dp;
					double thetab = te0*(1-sss)+ sss*iii;
					_assert_(thetab>=0 && thetab<=1);
					BR = VertexOnEdge(&R,eee,thetab);
					return  Gh.ProjectOnCurve(*eee,thetab,R,GR);
				}
			}
			// we find the end 
			if (v1 != pvB){
				if (( void*) v1 == pB)
				 tB = iii;

				double lg0 = lg;
				_assert_(eee);
				v1 = pvB;
				double dp = LengthInterpole(v0->m,v1->m,(R2) *v1 - (R2) *v0);
				lg += dp;	
				abscisse = lg*theta;
				if (abscisse <= lg && abscisse >= lg0 ) // small optimisation we know the lenght because end
				  { // ok we find the geom edge 
					double sss  =   (abscisse-lg0)/dp;
					double thetab = te0*(1-sss)+ sss*tB;
					_assert_(thetab>=0 && thetab<=1);
					BR = VertexOnEdge(&R,eee,thetab);
					return  Gh.ProjectOnCurve(*eee,thetab,R,GR);
				  }
			}
			abscisse = lg*theta;

		}
		_error_("Big bug...");
		return 0; // just for the compiler 
	}                  
	/*}}}*/
	void Mesh::ReconstructExistingMesh(BamgOpts* bamgopts){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/FillHoleInMesh)*/

		/*This routine reconstruct an existing mesh to make it CONVEX:
		 * -all the holes are filled
		 * -concave boundaries are filled
		 * A convex mesh is required for a lot of operations. This is why every mesh
		 * goes through this process.
		 * This routine also generates mesh properties such as adjencies,...
		 */

		/*Intermediary*/
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		// generation of the integer coordinate

		// find extrema coordinates of vertices pmin,pmax
		long i;
		if(verbose>2) _printf_("      Reconstruct mesh of " << nbv << " vertices\n"); 

		//initialize orderedvertices
		_assert_(orderedvertices);
		for (i=0;i<nbv;i++) orderedvertices[i]=0;

		//Initialize nbsubdomains
		nbsubdomains =0;

		/* generation of triangles adjacency*/

		//First add existing edges
		long kk =0;
		SetOfEdges4* edge4= new SetOfEdges4(nbt*3,nbv);
		for (i=0;i<nbe;i++){
			kk=kk+(i==edge4->SortAndAdd(GetId(edges[i][0]),GetId(edges[i][1])));
		}
		if (kk != nbe){ 
			_error_("There are " << kk-nbe << " double edges in the mesh");
		}

		//Add edges of all triangles in existing mesh
		long* st = new long[nbt*3];
		for (i=0;i<nbt*3;i++) st[i]=-1;
		for (i=0;i<nbt;i++){
			for (int j=0;j<3;j++){

				//Add current triangle edge to edge4
				long k =edge4->SortAndAdd(GetId(triangles[i][VerticesOfTriangularEdge[j][0]]),GetId(triangles[i][VerticesOfTriangularEdge[j][1]]));

				long invisible=triangles[i].Hidden(j);

				//If the edge has not been added to st, add it
				if(st[k]==-1) st[k]=3*i+j;

				//If the edge already exists, add adjacency
				else if(st[k]>=0) {
					_assert_(!triangles[i].TriangleAdj(j));
					_assert_(!triangles[st[k]/3].TriangleAdj((int) (st[k]%3)));

					triangles[i].SetAdj2(j,triangles+st[k]/3,(int)(st[k]%3));
					if (invisible) triangles[i].SetHidden(j);
					if (k<nbe)     triangles[i].SetLocked(j);

					//Make st[k] negative so that it will throw an error message if it is found again
					st[k]=-2-st[k]; 
				}

				//An edge belongs to 2 triangles
				else {
					_error_("The edge (" << GetId(triangles[i][VerticesOfTriangularEdge[j][0]]) << " , " << GetId(triangles[i][VerticesOfTriangularEdge[j][1]]) << ") belongs to more than 2 triangles");
				}
			}
		}

		//Display info if required
		if(verbose>5) {
			_printf_("         info of Mesh:\n");
			_printf_("            - number of vertices    = " << nbv << " \n"); 
			_printf_("            - number of triangles   = " << nbt << " \n"); 
			_printf_("            - number of given edges = " << nbe << " \n"); 
			_printf_("            - number of all edges   = " << edge4->nb() << "\n"); 
			_printf_("            - Euler number 1 - nb of holes = " << nbt-edge4->nb()+nbv << "\n"); 
		}

		//check the consistency of edge[].adj and the geometrical required vertex
		long k=0;
		for (i=0;i<edge4->nb();i++){
			if (st[i]>=0){ // edge alone 
				if (i<nbe){
					long i0=edge4->i(i);
					orderedvertices[i0] = vertices+i0;
					long i1=edge4->j(i);
					orderedvertices[i1] = vertices+i1;
				}
				else {
					k=k+1;
					if (k<10) {
						//print only 10 edges
						_printf_("Lost boundary edges " << i << " : " << edge4->i(i) << " " << edge4->j(i) << "\n");
					}
					else if (k==10){
						_printf_("Other lost boundary edges not shown...\n");
					}
				}
			}
		}
		if(k) {
			_error_(k << " boundary edges (from the geometry) are not defined as mesh edges");
		}

		/* mesh generation with boundary points*/
		long nbvb=0;
		for (i=0;i<nbv;i++){ 
			vertices[i].t=0;
			vertices[i].IndexInTriangle=0;
			if (orderedvertices[i]) orderedvertices[nbvb++]=orderedvertices[i];
		}

		Triangle* savetriangles=triangles;
		long savenbt=nbt;
		long savemaxnbt=maxnbt;
		SubDomain* savesubdomains=subdomains;
		subdomains=0;

		long  Nbtriafillhole=2*nbvb;
		Triangle* triafillhole=new Triangle[Nbtriafillhole];
		triangles = triafillhole;

		nbt=2;
		maxnbt= Nbtriafillhole;

		//Find a vertex that is not aligned with vertices 0 and 1
		for (i=2;det(orderedvertices[0]->i,orderedvertices[1]->i,orderedvertices[i]->i)==0;) 
		 if  (++i>=nbvb) {
			 _error_("ReconstructExistingMesh: All the vertices are aligned");
		 }
		//Move this vertex (i) to the 2d position in orderedvertices
		Exchange(orderedvertices[2], orderedvertices[i]);

		/*Reconstruct mesh beginning with 2 triangles*/
		BamgVertex *  v0=orderedvertices[0], *v1=orderedvertices[1];

		triangles[0](0) = NULL; // Infinite vertex
		triangles[0](1) = v0;
		triangles[0](2) = v1;

		triangles[1](0) = NULL;// Infinite vertex
		triangles[1](2) = v0;
		triangles[1](1) = v1;
		const int e0 = OppositeEdge[0];
		const int e1 = NextEdge[e0];
		const int e2 = PreviousEdge[e0];
		triangles[0].SetAdj2(e0, &triangles[1] ,e0);
		triangles[0].SetAdj2(e1, &triangles[1] ,e2);
		triangles[0].SetAdj2(e2, &triangles[1] ,e1);

		triangles[0].det = -1;  // boundary triangles
		triangles[1].det = -1;  // boundary triangles

		triangles[0].SetSingleVertexToTriangleConnectivity();
		triangles[1].SetSingleVertexToTriangleConnectivity();

		triangles[0].link=&triangles[1];
		triangles[1].link=&triangles[0];

		if (!quadtree) delete quadtree; //ReInitialise;
		quadtree = new BamgQuadtree(this,0);
		quadtree->Add(*v0);
		quadtree->Add(*v1);

		// vertices are added one by one
		long NbSwap=0;
		for (int icount=2; icount<nbvb; icount++) {
			BamgVertex *vi  = orderedvertices[icount];
			long long det3[3];
			Triangle *tcvi = TriangleFindFromCoord(vi->i,det3);
			quadtree->Add(*vi); 
			AddVertex(*vi,tcvi,det3);
			NbSwap += vi->Optim(1,1);
		}

		//enforce the boundary 
		AdjacentTriangle ta(0,0);
		long nbloss = 0,knbe=0;
		for ( i = 0; i < nbe; i++){
			if (st[i] >=0){ //edge alone => on border
				BamgVertex &a=edges[i][0], &b=edges[i][1];
				if (a.t && b.t){
					knbe++;
					if (ForceEdge(a,b,ta)<0) nbloss++;
				}
			}
		}
		if(nbloss) {
			_error_("we lost " << nbloss << " existing edges other " << knbe);
		}

		FindSubDomain(bamgopts,1);
		// remove all the hole 
		// remove all the good sub domain
		long krm =0;
		for (i=0;i<nbt;i++){
			if (triangles[i].link){ // remove triangles
				krm++;
				for (int j=0;j<3;j++){
					AdjacentTriangle ta =  triangles[i].Adj(j);
					Triangle &tta = *(Triangle*)ta;
					//if edge between remove and not remove 
					if(! tta.link){ 
						// change the link of ta;
						int ja = ta;
						BamgVertex *v0= ta.EdgeVertex(0);
						BamgVertex *v1= ta.EdgeVertex(1);
						long k =edge4->SortAndAdd(v0?GetId(v0):nbv,v1? GetId(v1):nbv);

						_assert_(st[k]>=0);
						tta.SetAdj2(ja,savetriangles + st[k] / 3,(int) (st[k]%3));
						ta.SetLock();
						st[k]=-2-st[k]; 
					}
				}
			}
		}
		long NbTfillHoll =0;
		for (i=0;i<nbt;i++){
			if (triangles[i].link) {
				triangles[i]=Triangle((BamgVertex *) NULL,(BamgVertex *) NULL,(BamgVertex *) NULL);
				triangles[i].color=-1;
			}
			else{
				triangles[i].color= savenbt+ NbTfillHoll++;
			}
		}
		_assert_(savenbt+NbTfillHoll<=savemaxnbt);

		// copy of the outside triangles in saveMesh 
		for (i=0;i<nbt;i++){
			if(triangles[i].color>=0) {
				savetriangles[savenbt]=triangles[i];
				savetriangles[savenbt].link=0;
				savenbt++;
			}
		}
		// gestion of the adj
		k =0;
		Triangle * tmax = triangles + nbt;
		for (i=0;i<savenbt;i++) { 
			Triangle & ti = savetriangles[i];
			for (int j=0;j<3;j++){
				Triangle * ta = ti.TriangleAdj(j);
				int aa = ti.NuEdgeTriangleAdj(j);
				int lck = ti.Locked(j);
				if (!ta) k++; // bug 
				else if ( ta >= triangles && ta < tmax){
					ta= savetriangles + ta->color;
					ti.SetAdj2(j,ta,aa);
					if(lck) ti.SetLocked(j);
				}
			}
		}

		// restore triangles;
		nbt=savenbt;
		maxnbt=savemaxnbt;
		delete [] triangles;
		delete [] subdomains;
		triangles = savetriangles;
		subdomains = savesubdomains;
		if (k) {
			_error_("number of triangles edges alone = " << k);
		}
		FindSubDomain(bamgopts);

		delete edge4;
		delete [] st;
		for (i=0;i<nbv;i++) quadtree->Add(vertices[i]);

		SetVertexFieldOn();

		/*Check requirements consistency*/
		for (i=0;i<nbe;i++){
			/*If the current mesh edge is on Geometry*/
			if(edges[i].GeomEdgeHook){
				for(int j=0;j<2;j++){
					/*Go through the edges adjacent to current edge (if on the same curve)*/
					if (!edges[i].adj[j]){
						/*The edge is on Geometry and does not have 2 adjacent edges... (not on a closed curve)*/
						/*Check that the 2 vertices are on geometry AND required*/
						if(!edges[i][j].GeomEdgeHook->IsRequiredVertex()){
							_printf_("ReconstructExistingMesh error message: problem with the edge number " << i+1 << ": [" << GetId(edges[i][0])+1 << " " << GetId(edges[i][1])+1 << "]\n");
							_printf_("This edge is on geometrical edge number " << Gh.GetId(edges[i].GeomEdgeHook)+1 << "\n");
							if (edges[i][j].GeomEdgeHook->OnGeomVertex())
							 _printf_("the vertex number " << GetId(edges[i][j])+1 << " of this edge is a geometric BamgVertex number " << Gh.GetId(edges[i][j].GeomEdgeHook->gv)+1 << "\n");
							else if (edges[i][j].GeomEdgeHook->OnGeomEdge())
							 _printf_("the vertex number " << GetId(edges[i][j])+1 << " of this edge is a geometric Edge number " << Gh.GetId(edges[i][j].GeomEdgeHook->ge)+1 << "\n");
							else
							 _printf_("Its pointer is " << edges[i][j].GeomEdgeHook << "\n");

							_printf_("This edge is on geometry and has no adjacent edge (open curve) and one of the tip is not required\n");
							_error_("See above (might be cryptic...)");
						}
					}
				}
			}
		}
	}
	/*}}}*/
	void Mesh::TrianglesRenumberBySubDomain(bool justcompress){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/ReNumberingTheTriangleBySubDomain)*/

		long *renu= new long[nbt];
		Triangle *t0,*t,*te=triangles+nbt;
		long k=0,it,i,j;

		for ( it=0;it<nbt;it++) 
		 renu[it]=-1; // outside triangle 
		for ( i=0;i<nbsubdomains;i++)
		  { 
			t=t0=subdomains[i].head;
			if (!t0){ // not empty sub domain
				_error_("!t0");
			}
			do { 
				long kt = GetId(t);
				if (kt<0 || kt >= nbt ){
					_error_("kt<0 || kt >= nbt");
				}
				if (renu[kt]!=-1){
					_error_("renu[kt]!=-1");
				}
				renu[kt]=k++;
			}
			while (t0 != (t=t->link));
		  }
		// take is same numbering if possible    
		if(justcompress)
		 for ( k=0,it=0;it<nbt;it++) 
		  if(renu[it] >=0 ) 
			renu[it]=k++;

		// put the outside triangles at the end
		for ( it=0;it<nbt;it++){
			if (renu[it]==-1) renu[it]=k++;
		}
		if (k != nbt){
			_error_("k != nbt");
		}
		// do the change on all the pointeur 
		for ( it=0;it<nbt;it++)
		 triangles[it].Renumbering(triangles,te,renu);

		for ( i=0;i<nbsubdomains;i++)
		 subdomains[i].head=triangles+renu[GetId(subdomains[i].head)];

		// move the Triangles  without a copy of the array 
		// be carefull not trivial code 
		for ( it=0;it<nbt;it++) // for all sub cycles of the permutation renu
		 if (renu[it] >= 0) // a new sub cycle
			{ 
			 i=it;
			 Triangle ti=triangles[i],tj;
			 while ( (j=renu[i]) >= 0) 
				{ // i is old, and j is new 
				 renu[i] = -1; // mark 
				 tj = triangles[j]; // save new
				 triangles[j]= ti; // new <- old
				 i=j;     // next 
				 ti = tj;
				}  
			}
		delete [] renu;

	}
	/*}}}*/
	void Mesh::SetIntCoor(const char * strfrom) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/SetIntCoor)*/

		/*Set integer coordinate for existing vertices*/

		//Get extrema coordinates of the existing vertices
		pmin =  vertices[0].r;
		pmax =  vertices[0].r;
		long i;
		for (i=0;i<nbv;i++) {
			pmin.x = Min(pmin.x,vertices[i].r.x);
			pmin.y = Min(pmin.y,vertices[i].r.y);
			pmax.x = Max(pmax.x,vertices[i].r.x);
			pmax.y = Max(pmax.y,vertices[i].r.y);
		}
		R2 DD = (pmax-pmin)*0.05;
		pmin = pmin-DD;
		pmax = pmax+DD; 

		//Compute coefIcoor
		long MaxICoord = 1073741823; //2^30 - 1 = =111...111 (29 times one)
		coefIcoor= (MaxICoord)/(Max(pmax.x-pmin.x,pmax.y-pmin.y));
		if (coefIcoor<=0){
			_error_("coefIcoor should be positive, a problem in the geometry is likely");
		}

		// generation of integer coord  
		for (i=0;i<nbv;i++) {
			vertices[i].i = R2ToI2(vertices[i].r);    
		}

		// computation of the det 
		int number_of_errors=0;
		for (i=0;i<nbt;i++) {
			BamgVertex* v0 = triangles[i](0);
			BamgVertex* v1 = triangles[i](1);
			BamgVertex* v2 = triangles[i](2);

			//If this is not a boundary triangle
			if (v0 && v1 && v2){

				/*Compute determinant*/
				triangles[i].det= det(v0->GetIntegerCoordinates(),v1->GetIntegerCoordinates(),v2->GetIntegerCoordinates());

				/*Check that determinant is positive*/
				if (triangles[i].det <=0){

					/*increase number_of_errors and print error only for the first 20 triangles*/
					number_of_errors++;
					if (number_of_errors<20){
						_printf_("Area of Triangle " << i+1 << " < 0 (det=" << triangles[i].det << ")\n");
					}
				}
			}

			//else, set as -1
			else triangles[i].det=-1;
		}

		if (number_of_errors) _error_("Fatal error: some triangles have negative areas, see above");
	}
	/*}}}*/
	void Mesh::SmoothingVertex(BamgOpts* bamgopts,int nbiter,double omega ) { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/SmoothingVertex)*/

		/*Intermediaries*/
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//  if quatree exist remove it end reconstruct
		if (quadtree) delete quadtree;
		quadtree=0;
		CreateSingleVertexToTriangleConnectivity();
		Triangle vide; // a triangle to mark the boundary vertex
		Triangle   ** tstart= new Triangle* [nbv];
		long i,j,k;
		//   attention si Background == Triangle alors on ne peut pas utiliser la rechech rapide 
		if ( this == & BTh)
		 for ( i=0;i<nbv;i++)
		  tstart[i]=vertices[i].t;     
		else 
		 for ( i=0;i<nbv;i++)
		  tstart[i]=0;
		for ( j=0;j<NbVerticesOnGeomVertex;j++ ) 
		 tstart[ GetId(VerticesOnGeomVertex[j].meshvertex)]=&vide;
		for ( k=0;k<NbVerticesOnGeomEdge;k++ ) 
		 tstart[ GetId(VerticesOnGeomEdge[k].meshvertex)]=&vide;
		if(verbose>2) _printf_("   SmoothingVertex: nb Iteration = " << nbiter << ", Omega=" << omega << "\n");
		for (k=0;k<nbiter;k++)
		  {
			long i,NbSwap =0;
			double delta =0;
			for ( i=0;i<nbv;i++)
			 if (tstart[i] != &vide) // not a boundary vertex 
			  delta=Max(delta,vertices[i].Smoothing(*this,BTh,tstart[i],omega));
			for ( i=0;i<nbv;i++)
			 if (tstart[i] != &vide) // not a boundary vertex 
			  NbSwap += vertices[i].Optim(1);
			if (verbose>3) _printf_("      move max = " << pow(delta,0.5) << ", iteration = " << k << ", nb of swap = " << NbSwap << "\n");
		  }

		delete [] tstart;
		if (quadtree) quadtree= new BamgQuadtree(this);
	}
	/*}}}*/
	void Mesh::SmoothMetric(BamgOpts* bamgopts,double raisonmax) { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/SmoothMetric)*/

		/*Intermediaries*/
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		if(raisonmax<1.1) return;
		if(verbose > 1) _printf_("   Mesh::SmoothMetric raisonmax = " << raisonmax << "\n");
		CreateSingleVertexToTriangleConnectivity();
		long i,j,kch,kk,ip;
		long *first_np_or_next_t0 = new long[nbv];
		long *first_np_or_next_t1 = new long[nbv];
		long Head0 =0,Head1=-1;
		double logseuil= log(raisonmax);

		for(i=0;i<nbv-1;i++)
		 first_np_or_next_t0[i]=i+1; 
		first_np_or_next_t0[nbv-1]=-1;// end;
		for(i=0;i<nbv;i++)
		 first_np_or_next_t1[i]=-1;
		kk=0;
		while(Head0>=0&& kk++<100){
			kch=0;
			for(i=Head0;i>=0;i=first_np_or_next_t0[ip=i],first_np_or_next_t0[ip]=-1) {
				//  pour tous les triangles autour du sommet s
				Triangle* t= vertices[i].t;
				if (!t){
					_error_("!t");
				}
				BamgVertex & vi = vertices[i];
				AdjacentTriangle ta(t,EdgesVertexTriangle[vertices[i].IndexInTriangle][0]);
				BamgVertex *pvj0 = ta.EdgeVertex(0);
				while (1) {
					ta=Previous(Adj(ta));
					if (vertices+i != ta.EdgeVertex(1)){
						_error_("vertices+i != ta.EdgeVertex(1)");
					}
					BamgVertex *pvj = (ta.EdgeVertex(0));
					BamgVertex & vj = *pvj;
					if(pvj){
						j= &vj-vertices;
						if (j<0 || j >= nbv){
							_error_("j<0 || j >= nbv");
						}
						R2 Aij = (R2) vj - (R2) vi;
						double ll =  Norme2(Aij);
						if (0) {  
							double hi = ll/vi.m.Length(Aij.x,Aij.y);
							double hj = ll/vj.m.Length(Aij.x,Aij.y);
							if(hi < hj)
							  {
								double dh=(hj-hi)/ll;
								if (dh>logseuil) {
									vj.m.IntersectWith(vi.m/(1 +logseuil*ll/hi));
									if(first_np_or_next_t1[j]<0)
									 kch++,first_np_or_next_t1[j]=Head1,Head1=j;
								}
							  }
						} 
						else
						  {
							double li = vi.m.Length(Aij.x,Aij.y);
							if( vj.m.IntersectWith(vi.m/(1 +logseuil*li)) )
							 if(first_np_or_next_t1[j]<0) // if the metrix change 
							  kch++,first_np_or_next_t1[j]=Head1,Head1=j;
						  }
					}
					if  ( &vj ==  pvj0 ) break;
				}
			}
			Head0 = Head1;
			Head1 = -1;
			Exchange(first_np_or_next_t0,first_np_or_next_t1);
		}
		if(verbose>2) _printf_("      number of iterations = " << kch << "\n"); 
		delete [] first_np_or_next_t0;
		delete [] first_np_or_next_t1;
	}
	/*}}}*/
	long  Mesh::SplitInternalEdgeWithBorderVertices(){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/SplitInternalEdgeWithBorderVertices)*/

		long NbSplitEdge=0;
		SetVertexFieldOn();  
		long it;
		long nbvold=nbv;
		long int verbose=2;
		for (it=0;it<nbt;it++){
			Triangle &t=triangles[it];
			if (t.link)
			 for (int j=0;j<3;j++)
			  if(!t.Locked(j) && !t.Hidden(j)){

				  Triangle *ptt = t.TriangleAdj(j);
				  Triangle &tt = *ptt;

				  if( ptt && tt.link && it < GetId(tt)) 
					 { // an internal edge 
					  BamgVertex &v0 = t[VerticesOfTriangularEdge[j][0]];
					  BamgVertex &v1 = t[VerticesOfTriangularEdge[j][1]];
					  if (v0.GeomEdgeHook && v1.GeomEdgeHook){
						  R2 P= ((R2) v0 + (R2) v1)*0.5;
						  if ( nbv<maxnbv) {
							  vertices[nbv].r = P;
							  vertices[nbv++].m = Metric(0.5,v0.m,0.5,v1.m);
							  vertices[nbv].ReferenceNumber=0;
						  }
						  NbSplitEdge++;
					  }
					 }
			  }
		}
		CreateSingleVertexToTriangleConnectivity();    
		if (nbvold!=nbv){
			long  iv = nbvold;
			long NbSwap = 0;
			long long det3[3];  
			for (int i=nbvold;i<nbv;i++) {// for all the new point
				BamgVertex & vi = vertices[i];
				vi.i = R2ToI2(vi.r);
				vi.r = I2ToR2(vi.i);

				// a good new point 
				vi.ReferenceNumber=0; 
				Triangle *tcvi = TriangleFindFromCoord(vi.i,det3);
				if (tcvi && !tcvi->link) {
					_printf_("problem inserting point in SplitInternalEdgeWithBorderVertices (tcvj && !tcvj->link)\n");
				}

				quadtree->Add(vi);
				if (!tcvi || tcvi->det<0){// internal
					_error_("!tcvi || tcvi->det < 0");
				}
				AddVertex(vi,tcvi,det3);
				NbSwap += vi.Optim(1);          
				iv++;
			}
			if (verbose>3) {
				_printf_("   number of points: " << iv << "\n");
				_printf_("   number of swap to  split internal edges with border vertices: " << NbSwap << "\n");
				nbv = iv;
			}
		}
		if (NbSplitEdge>nbv-nbvold) _printf_("WARNING: not enough vertices to split all internal edges, we lost " << NbSplitEdge - ( nbv-nbvold) << " edges. To fix this, increase 'maxnbv' (see bamg help)\n");
		if (verbose>2) _printf_("SplitInternalEdgeWithBorderVertices: Number of splited edge " << NbSplitEdge << "\n");

		return  NbSplitEdge;
	}
	/*}}}*/
	I2 Mesh::R2ToI2(const R2 & P) const {/*{{{*/
		return  I2( (int) (coefIcoor*(P.x-pmin.x)),(int) (coefIcoor*(P.y-pmin.y)) );
	}
	/*}}}*/
	R2 Mesh::I2ToR2(const I2 & P) const {/*{{{*/
		return  R2( (double) P.x/coefIcoor+pmin.x, (double) P.y/coefIcoor+pmin.y);
	}
	/*}}}*/
	Triangle * Mesh::TriangleFindFromCoord(const I2 & B,long long det3[3], Triangle *tstart){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/FindTriangleContening)*/

		Triangle * t=0;	
		int j,jp,jn,jj;
		int counter;

		/*Get starting triangle. Take tsart if provided*/
		if (tstart) t=tstart;

		/*Else find the closest Triangle using the quadtree*/
		else {

			/*Check that the quadtree does exist*/
			if (!quadtree) _error_("no starting triangle provided and no quadtree available");

			/*Call NearestVertex*/
			BamgVertex *a = quadtree->NearestVertex(B.x,B.y) ;

			/*Check output (Vertex a)*/
			if (!a)    _error_("problem while trying to find nearest vertex from a given point. No output found");
			if (!a->t) _error_("no triangle is associated to vertex number " << GetId(a)+1 << " (orphan?)");
			_assert_(a>=vertices && a<vertices+nbv);

			/*Get starting triangle*/
			t = a->t;
			_assert_(t>=triangles && t<triangles+nbt);
		}

		long long  detop ;

		/*initialize number of test triangle*/
		counter=0; 

		/*The initial triangle might be outside*/
		while (t->det < 0){ 

			/*Get a real vertex from this triangle (k0)*/
			int k0=(*t)(0)?(((*t)(1)?((*t)(2)?-1:2):1)):0;
			_assert_(k0>=0);// k0 the NULL vertex
			int k1=NextVertex[k0],k2=PreviousVertex[k0];
			det3[k0]=det(B,(*t)[k1],(*t)[k2]);
			det3[k1]=det3[k2]=-1;     
			if (det3[k0] > 0) // outside B 
			 return t; 
			t = t->TriangleAdj(OppositeEdge[k0]);
			counter++;
			_assert_(counter<2);
		}

		jj=0;
		detop = det(*(*t)(VerticesOfTriangularEdge[jj][0]),*(*t)(VerticesOfTriangularEdge[jj][1]),B);

		while(t->det>0){

			/*Increase counter*/
			if (++counter>=10000) _error_("Maximum number of iteration reached (threshold = " << counter << ").");

			j= OppositeVertex[jj];
			det3[j] = detop;  //det(*b,*s1,*s2);
			jn = NextVertex[j];
			jp = PreviousVertex[j];
			det3[jp]= det(*(*t)(j),*(*t)(jn),B);
			det3[jn] = t->det-det3[j] -det3[jp];

			// count the number k of  det3 <0
			int k=0,ii[3];
			if (det3[0] < 0 ) ii[k++]=0; 
			if (det3[1] < 0 ) ii[k++]=1;
			if (det3[2] < 0 ) ii[k++]=2;
			// 0 => ok
			// 1 => go in way 1
			// 2 => two way go in way 1 or 2 randomly

			if (k==0) break;
			if (k==2 && this->RandomNumber(1)) Exchange(ii[0],ii[1]);
			_assert_(k<3);
			AdjacentTriangle t1 = t->Adj(jj=ii[0]);
			if ((t1.det() < 0 ) && (k == 2))
			 t1 = t->Adj(jj=ii[1]);
			t=t1;
			j=t1;// for optimisation we now the -det[OppositeVertex[j]];
			detop = -det3[OppositeVertex[jj]];
			jj = j;
		}

		if (t->det<0) // outside triangle 
		 det3[0]=det3[1]=det3[2]=-1,det3[OppositeVertex[jj]]=detop;
		return t;
	}
	/*}}}*/
	void Mesh::TriangleIntNumbering(long* renumbering){/*{{{*/

		long num=0;
		for (int i=0;i<nbt;i++){
			if (triangles[i].det>0) renumbering[i]=num++;
			else renumbering[i]=-1;
		}
		return;   
	}
	/*}}}*/
	long  Mesh::TriangleReferenceList(long* reft) const {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/ConsRefTriangle)*/

		Triangle *t0,*t;
		long k=0, num;   

		//initialize all triangles as -1 (outside)
		for (int it=0;it<nbt;it++) reft[it]=-1;

		//loop over all subdomains
		for (int i=0;i<nbsubdomains;i++){ 

			//first triangle of the subdomain i
			t=t0=subdomains[i].head;

			//check that the subdomain is not empty
			if (!t0){ _error_("At least one subdomain is empty");}

			//loop
			do{
				k++;

				//get current triangle number
				num = GetId(t);

				//check that num is in [0 nbt[
				_assert_(num>=0 && num<nbt);

				//reft of this triangle is the subdomain number
				reft[num]=i;

			} while (t0 != (t=t->link));
			//stop when all triangles of subdomains have been tagged

		}
		return k;   
	}
	/*}}}*/
	void Mesh::Triangulate(double* x,double* y,int nods,BamgOpts* bamgopts){/*{{{*/

		int verbose=0;
		int i;
		Metric M1(1);

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		/*Initialize mesh*/
		Init(nods);//this resets nbv to 0
		nbv=nods;

		//Vertices
		if(verbose) _printf_("Reading vertices (" << nbv << ")\n");
		for(i=0;i<nbv;i++){
			vertices[i].r.x=x[i];
			vertices[i].r.y=y[i];
			vertices[i].ReferenceNumber=1;
			vertices[i].m=M1;
			vertices[i].color=0;
		}
		maxnbt=2*maxnbv-2; // for filling The Holes and quadrilaterals 

		/*Insert Vertices*/
		Insert(bamgopts);
	}
	/*}}}*/
	void Mesh::TriangulateFromGeom0(BamgOpts* bamgopts){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/GeomToTriangles0)*/
		/*Generate mesh from geometry*/

		/*Intermediaries*/
		int                i,k;
		int					 verbose=0;
		int                nbcurves    = 0;
		int                NbNewPoints,NbEdgeCurve;
		double             lcurve,lstep,s;
		const int          MaxSubEdge  = 10;

		R2          AB;
		GeomVertex *a, *b;
		BamgVertex *va,*vb;
		GeomEdge   *e;

		// add a ref to GH to make sure that it is not destroyed by mistake
		Gh.NbRef++;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		//build background mesh flag (1 if background, else 0)
		bool background=(&BTh != this);

		/*Build VerticesOnGeomVertex*/

		//Compute the number of geometrical vertices that we are going to use to mesh
		for (i=0;i<Gh.nbv;i++){
			if (Gh[i].Required()) NbVerticesOnGeomVertex++;
		}
		//allocate
		VerticesOnGeomVertex = new VertexOnGeom[NbVerticesOnGeomVertex];  
		if(NbVerticesOnGeomVertex >= maxnbv) _error_("too many vertices on geometry: " << NbVerticesOnGeomVertex << " >= " << maxnbv);
		_assert_(nbv==0);
		//Build VerticesOnGeomVertex
		for (i=0;i<Gh.nbv;i++){
			/* Add vertex only if required*/
			if (Gh[i].Required()) {//Gh  vertices Required

				//Add the vertex
				_assert_(nbv<maxnbv);
				vertices[nbv]=Gh[i];

				//Add pointer from geometry (Gh) to vertex from mesh (Th)
				Gh[i].MeshVertexHook=vertices+nbv;

				//Build VerticesOnGeomVertex for current point
				VerticesOnGeomVertex[nbv]=VertexOnGeom(*Gh[i].MeshVertexHook,Gh[i]);

				//nbv increment
				nbv++;
			}
		}

		/*Build VerticesOnGeomEdge*/

		//check that edges is still empty (Init)
		_assert_(!edges);

		/* Now we are going to create the first edges corresponding
		 * to the one present in the geometry provided.
		 * We proceed in 2 steps
		 *  -step 0: we count all the edges
		 *           we allocate the number of edges at the end of step 0
		 *  -step 1: the edges are created */
		for (int step=0;step<2;step++){

			//initialize number of edges and number of edges max
			long nbex=0;
			nbe=0;
			long NbVerticesOnGeomEdge0=NbVerticesOnGeomEdge;
			Gh.UnMarkEdges();	
			nbcurves=0;

			//go through the edges of the geometry
			for (i=0;i<Gh.nbe;i++){

				//ei = current Geometrical edge
				GeomEdge &ei=Gh.edges[i];   

				//loop over the two vertices of the edge ei
				for(int j=0;j<2;j++) {

					/*Take only required vertices (corner->beginning of a new curve)*/
					if (!ei.Mark() && ei[j].Required()){ 

						long  nbvend=0;
						Edge* PreviousNewEdge=NULL;
						lstep = -1;

						/*If Edge is required (do that only once for the 2 vertices)*/
						if(ei.Required()){
							if (j==0){
								//do not create internal points if required (take it as is)
								if(step==0) nbe++;
								else{ 
									e=&ei;
									a=ei(0);
									b=ei(1);

									//check that edges has been allocated
									_assert_(edges);
									edges[nbe].v[0]=a->MeshVertexHook;
									edges[nbe].v[1]=b->MeshVertexHook;;
									edges[nbe].ReferenceNumber = e->ReferenceNumber;
									edges[nbe].GeomEdgeHook = e;
									edges[nbe].adj[0] = 0;
									edges[nbe].adj[1] = 0;
									nbe++;
								}
							}
						}

						/*If Edge is not required: we are on a curve*/
						else {
							for (int kstep=0;kstep<=step;kstep++){
								//kstep=0, compute number of edges (discretize curve)
								//kstep=1  create the points and edge
								PreviousNewEdge=0;
								NbNewPoints=0;
								NbEdgeCurve=0;
								if (nbvend>=maxnbv) _error_("maximum number of vertices too low! Check the domain outline or increase maxnbv");
								lcurve =0;
								s = lstep; //-1 initially, then length of each sub edge

								/*reminder: i = edge number, j=[0;1] vertex index in edge*/
								k=j;            // k = vertex index in edge (0 or 1)
								e=&ei;          // e = reference of current edge
								a=ei(k);        // a = pointer toward the kth vertex of the current edge
								va = a->MeshVertexHook; // va = pointer toward mesh vertex associated
								e->SetMark();   // Mark edge

								/*Loop until we reach the end of the curve*/
								for(;;){ 
									k = 1-k;            // other vertx index of the curve
									b = (*e)(k);        // b = pointer toward the other vertex of the current edge
									AB= b->r - a->r;   // AB = vector of the current edge
									Metric MA = background ? BTh.MetricAt(a->r) :a->m ;  //Get metric associated to A
									Metric MB = background ? BTh.MetricAt(b->r) :b->m ;  //Get metric associated to B
									double ledge = (MA.Length(AB.x,AB.y) + MB.Length(AB.x,AB.y))/2;                  //Get edge length in metric

									/* We are now creating the mesh edges from the geometrical edge selected above.
									 * The edge will be divided according to the metric previously computed and cannot
									 * be divided more than 10 times (MaxSubEdge). */

									//By default, there is only one subedge that is the geometrical edge itself
									int NbSubEdge = 1;

									//initialize lSubEdge, holding the length of each subedge (cannot be higher than 10)
									double lSubEdge[MaxSubEdge];

									//Build Subedges according to the edge length
									if (ledge < 1.5){
										//if ledge < 1.5 (between one and 2), take the edge as is
										lSubEdge[0] = ledge;
									}
									//else, divide the edge
									else {
										//compute number of subedges (division of the edge), Maximum is 10
										NbSubEdge = Min( MaxSubEdge, (int) (ledge +0.5));
										/*Now, we are going to divide the edge according to the metric.
										 * Get segment by sement along the edge.
										 * Build lSubEdge, which holds the distance between the first vertex
										 * of the edge and the next point on the edge according to the 
										 * discretization (each SubEdge is AB)*/
										R2 A,B;
										A=a->r;
										Metric MAs=MA,MBs;
										ledge=0; 
										double x =0, xstep= 1./NbSubEdge;
										for (int kk=0; kk < NbSubEdge; kk++,A=B,MAs=MBs ) {
											x += xstep;
											B =  e->F(k ? x : 1-x);
											MBs= background ? BTh.MetricAt(B) : Metric(1-x,MA,x,MB);
											AB = A-B;
											lSubEdge[kk]=(ledge+=(MAs.Length(AB.x,AB.y)+MBs.Length(AB.x,AB.y))/2);
										}
									}

									double lcurveb = lcurve+ledge;

									/*Now, create corresponding points*/
									while(s>=lcurve && s<=lcurveb && nbv<nbvend){

										/*Schematic of current curve
										 *
										 *  a                   vb                  b          // vertex
										 *  0              ll0     ll1              ledge      // length from a
										 *  + --- + - ... - + --S-- + --- + - ... - +          // where is S
										 *  0              kk0     kk1              NbSubEdge  // Sub edge index
										 *
										 */

										double ss = s-lcurve;

										/*Find the SubEdge containing ss using Dichotomy*/
										int kk0=-1,kk1=NbSubEdge-1,kkk;
										double ll0=0,ll1=ledge,llk;
										while (kk1-kk0>1){
											if (ss < (llk=lSubEdge[kkk=(kk0+kk1)/2]))
											 kk1=kkk,ll1=llk;
											else
											 kk0=kkk,ll0=llk;
										}
										_assert_(kk1!=kk0);

										/*Curvilinear coordinate in [0 1] of ss in current edge*/
										// WARNING: This is what we would do
										// ssa = (ss-ll0)/(ll1-ll0);
										// aa = (kk0+ssa)/NbSubEdge
										// This is what Bamg does:
										double sbb = (ss-ll0)/(ll1-ll0);
										/*Curvilinear coordinate in [0 1] of ss in current curve*/
										double bb = (kk1+sbb)/NbSubEdge;
										double aa = 1-bb;

										// new vertex on edge
										vb = &vertices[nbv++];
										vb->m = Metric(aa,a->m,bb,b->m);
										vb->ReferenceNumber = e->ReferenceNumber;
										double abcisse = k ? bb : aa;
										vb->r =  e->F(abcisse);
										VerticesOnGeomEdge[NbVerticesOnGeomEdge++]= VertexOnGeom(*vb,*e,abcisse);        

										// to take into account the direction of the edge
										s += lstep;
										edges[nbe].v[0]=va;
										edges[nbe].v[1]=vb;
										edges[nbe].ReferenceNumber =e->ReferenceNumber;
										edges[nbe].GeomEdgeHook = e;
										edges[nbe].adj[0] = PreviousNewEdge;
										if(PreviousNewEdge) PreviousNewEdge->adj[1]=&edges[nbe];
										PreviousNewEdge=edges+nbe;
										nbe++;
										va = vb;
									}

									/*We just added one edge to the curve: Go to the next one*/
									lcurve = lcurveb;
									e->SetMark();
									a=b;

									/*If b is required, we are on a new curve->break*/
									if (b->Required()) break;
									int kprev=k;
									k = e->AdjVertexIndex[kprev];// next vertices
									e = e->Adj[kprev];
									_assert_(e);
								}// for(;;)
								vb = b->MeshVertexHook;

								/*Number of edges in the last disretized curve*/
								NbEdgeCurve = Max((long) (lcurve +0.5), (long) 1);
								/*Number of internal vertices in the last disretized curve*/
								NbNewPoints = NbEdgeCurve-1;
								if(!kstep){
									NbVerticesOnGeomEdge0 += NbNewPoints;
									nbcurves++;
								}
								nbvend=nbv+NbNewPoints; 
								lstep = lcurve / NbEdgeCurve; //approximately one
							}// end of curve --
							if (edges) { // last edges of the curves 
								edges[nbe].v[0]=va;
								edges[nbe].v[1]=vb;
								edges[nbe].ReferenceNumber = e->ReferenceNumber;
								edges[nbe].GeomEdgeHook = e;
								edges[nbe].adj[0] = PreviousNewEdge;
								edges[nbe].adj[1] = 0;
								if(PreviousNewEdge) PreviousNewEdge->adj[1] = & edges[nbe];
								nbe++;
							}
							else nbe += NbEdgeCurve;
						} // end on  curve ---
					}
				}
			} // for (i=0;i<nbe;i++)
			if(!step) {
				_assert_(!edges);
				_assert_(!VerticesOnGeomEdge);

				edges = new Edge[nbex=nbe];
				if(NbVerticesOnGeomEdge0) VerticesOnGeomEdge = new VertexOnGeom[NbVerticesOnGeomEdge0];

				// do the vertex on a geometrical vertex
				_assert_(VerticesOnGeomEdge || NbVerticesOnGeomEdge0==0);
				NbVerticesOnGeomEdge0 = NbVerticesOnGeomEdge;       
			}
			else{
				_assert_(NbVerticesOnGeomEdge==NbVerticesOnGeomEdge0);
			}
		}

		//Insert points inside existing triangles
		if (verbose>4) _printf_("      -- current number of vertices = " << nbv << "\n");
		if (verbose>3) _printf_("      Creating initial Constrained Delaunay Triangulation...\n");
		if (verbose>3) _printf_("         Inserting boundary points\n");
		Insert(bamgopts);

		//Force the boundary
		if (verbose>3) _printf_("         Forcing boundaries\n");
		ForceBoundary(bamgopts);

		//Extract SubDomains
		if (verbose>3) _printf_("         Extracting subdomains\n");
		FindSubDomain(bamgopts);

		if (verbose>3) _printf_("      Inserting internal points\n");
		NewPoints(*this,bamgopts,0) ;
		if (verbose>4) _printf_("      -- current number of vertices = " << nbv << "\n");
	}
	/*}}}*/
	void Mesh::TriangulateFromGeom1(BamgOpts* bamgopts,int KeepVertices){ /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/GeomToTriangles1)*/

		/*Intermediaries*/
		int verbose=0;

		/*Get options*/
		if(bamgopts) verbose=bamgopts->verbose;

		Gh.NbRef++;// add a ref to Gh

		/************************************************************************* 
		 * method in 2 steps
		 * 1 - compute the number of new edges to allocate
		 * 2 - construct the edges
		 * remark:
		 * in this part we suppose to have a background mesh with the same geometry 
		 * 
		 * To construct the discretization of the new mesh we have to 
		 * rediscretize the boundary of background Mesh 
		 * because we have only the pointeur from the background mesh to the geometry.
		 * We need the abcisse of the background mesh vertices on geometry
		 * so a vertex is 
		 * 0 on GeomVertex ;
		 * 1 on GeomEdge + abcisse
		 * 2 internal 
		 *************************************************************************/

		//Check that background mesh and current mesh do have the same geometry
		_assert_(&BTh.Gh==&Gh);
		BTh.NbRef++; // add a ref to BackGround Mesh

		//Initialize new mesh
		BTh.SetVertexFieldOn();
		int* bcurve = new int[Gh.nbcurves]; // 

		/* There are 2 ways to make the loop 
		 * 1) on the geometry 
		 * 2) on the background mesh
		 *  if you do the loop on geometry, we don't have the pointeur on background,
		 *  and if you do the loop in background we have the pointeur on geometry
		 * so do the walk on  background */

		NbVerticesOnGeomVertex=0;
		NbVerticesOnGeomEdge=0;

		/*STEP 1 copy of Required vertices*/

		int i; 
		for (i=0;i<Gh.nbv;i++) if (Gh[i].Required()) NbVerticesOnGeomVertex++;
		printf("\n");
		if(NbVerticesOnGeomVertex >= maxnbv){
			delete [] bcurve;
			_error_("too many vertices on geometry: " << NbVerticesOnGeomVertex << " >= " << maxnbv);
		}

		VerticesOnGeomVertex = new VertexOnGeom[  NbVerticesOnGeomVertex];
		VertexOnBThVertex    = new VertexOnVertex[NbVerticesOnGeomVertex];

		//At this point there is NO vertex but vertices should have been allocated by Init
		_assert_(vertices);
		for (i=0;i<Gh.nbv;i++){
			if (Gh[i].Required()) {//Gh vertices Required
				vertices[nbv]  =Gh[i];
				vertices[nbv].i=I2(0,0);
				Gh[i].MeshVertexHook = vertices + nbv;// save Geom -> Th
				VerticesOnGeomVertex[nbv]= VertexOnGeom(vertices[nbv],Gh[i]);
				nbv++;
			}
			else Gh[i].MeshVertexHook=0;
		} 
		for (i=0;i<BTh.NbVerticesOnGeomVertex;i++){ 
			VertexOnGeom &vog=BTh.VerticesOnGeomVertex[i];
			if (vog.IsRequiredVertex()){
				GeomVertex* gv=vog;
				BamgVertex *bv = vog;
				_assert_(gv->MeshVertexHook); // use of Geom -> Th
				VertexOnBThVertex[NbVertexOnBThVertex++]=VertexOnVertex(gv->MeshVertexHook,bv);
				gv->MeshVertexHook->m = bv->m; // for taking the metric of the background mesh
			}
		}
		_assert_(NbVertexOnBThVertex==NbVerticesOnGeomVertex); /*This might be due to MaxCornerAngle too small*/

		/*STEP 2: reseed boundary edges*/

		//  find the begining of the curve in BTh
		Gh.UnMarkEdges();	
		int bfind=0;
		for (int i=0;i<Gh.nbcurves;i++) bcurve[i]=-1; 

		/*Loop over the backgrounf mesh BTh edges*/
		for (int iedge=0;iedge<BTh.nbe;iedge++){      
			Edge &ei = BTh.edges[iedge];

			/*Loop over the 2 vertices of the current edge*/
			for(int je=0;je<2;je++){

				/* If one of the vertex is required we are in a new curve*/
				if (ei[je].GeomEdgeHook->IsRequiredVertex()){ 

					/*Get curve number*/
					int nc=ei.GeomEdgeHook->CurveNumber;

					//_printf_("Dealing with curve number " << nc << "\n");
					//_printf_("edge on geometry is same as GhCurve? " << (ei.GeomEdgeHook==Gh.curves[nc].FirstEdge || ei.GeomEdgeHook==Gh.curves[nc].LastEdge)?"yes":"no\n");
					//if(ei.GeomEdgeHook==Gh.curves[nc].FirstEdge || ei.GeomEdgeHook==Gh.curves[nc].LastEdge){
					//	_printf_("Do we have the right extremity? curve first vertex -> " << ((GeomVertex *)*ei[je].GeomEdgeHook==&(*Gh.curves[nc].FirstEdge)[Gh.curves[nc].FirstVertexIndex])?"yes":"no\n");
					//	_printf_("Do we have the right extremity? curve last  vertex -> " << ((GeomVertex *)*ei[je].GeomEdgeHook==&(*Gh.curves[nc].LastEdge)[Gh.curves[nc].LastVertexIndex])?"yes":"no\n");
					//}
					//BUG FIX from original bamg
					/*Check that we are on the same edge and right vertex (0 or 1) */
					if(ei.GeomEdgeHook==Gh.curves[nc].FirstEdge  && (GeomVertex *)*ei[je].GeomEdgeHook==&(*Gh.curves[nc].FirstEdge)[Gh.curves[nc].FirstVertexIndex]){
						bcurve[nc]=iedge*2+je;
						bfind++;	
					}
					else if ((ei.GeomEdgeHook==Gh.curves[nc].LastEdge  && (GeomVertex *)*ei[je].GeomEdgeHook==&(*Gh.curves[nc].LastEdge)[Gh.curves[nc].LastVertexIndex]) && bcurve[nc]==-1){
						bcurve[nc]=iedge*2+je;
						bfind++;	
					}
				}
			}
		} 
		if (bfind!=Gh.nbcurves){
			delete [] bcurve;
			_error_("problem generating number of curves (" << Gh.nbcurves << " found in the geometry but " << bfind << " curve found in the mesh)");
		}

		// method in 2 + 1 step 
		//  0.0) compute the length and the number of vertex to do allocation
		//  1.0) recompute the length
		//  1.1) compute the  vertex 

		long nbex=0,NbVerticesOnGeomEdgex=0;
		for (int step=0; step <2;step++){

			long NbOfNewPoints=0;
			long NbOfNewEdge=0;
			long iedge;
			Gh.UnMarkEdges();	
			double L=0;

			/*Go through all geometrical curve*/
			for (int icurve=0;icurve<Gh.nbcurves;icurve++){ 

				/*Get edge and vertex (index) of background mesh on this curve*/
				iedge=bcurve[icurve]/2;
				int jedge=bcurve[icurve]%2;

				/*Get edge of Bth with index iedge*/
				Edge &ei = BTh.edges[iedge];

				/*Initialize variables*/
				double Lstep=0;             // step between two points   (phase==1) 
				long NbCreatePointOnCurve=0;// Nb of new points on curve (phase==1) 

				/*Do phase 0 to step*/
				for(int phase=0;phase<=step;phase++){

					/*Current curve pointer*/
					Curve *curve= Gh.curves+icurve;

					/*Get index of current curve*/
					int icurveequi= Gh.GetId(curve);

					/*For phase 0, check that we are at the begining of the curve only*/
					if(phase==0 &&  icurveequi!=icurve)  continue;

					int   k0=jedge,k1;
					Edge* pe=  BTh.edges+iedge;
					int   iedgeequi=bcurve[icurveequi]/2;
					int   jedgeequi=bcurve[icurveequi]%2;

					int k0equi=jedgeequi,k1equi;		  
					Edge * peequi= BTh.edges+iedgeequi;
					GeomEdge *ongequi = peequi->GeomEdgeHook;

					double sNew=Lstep;// abscisse of the new points (phase==1) 
					L=0;// length of the curve
					long i=0;// index of new points on the curve
					GeomVertex * GA0 = *(*peequi)[k0equi].GeomEdgeHook;
					BamgVertex *A0;
					A0 = GA0->MeshVertexHook;  // the vertex in new mesh
					BamgVertex *A1;
					VertexOnGeom *GA1;
					Edge* PreviousNewEdge = 0;

					// New Curve phase 
					_assert_(A0-vertices>=0 && A0-vertices<nbv);
					if(ongequi->Required()){
						GeomVertex *GA1 = *(*peequi)[1-k0equi].GeomEdgeHook;
						A1 = GA1->MeshVertexHook;  //
					}       
					else {
						for(;;){
							Edge &ee=*pe; 
							Edge &eeequi=*peequi; 
							k1 = 1-k0; // next vertex of the edge 
							k1equi= 1 - k0equi;
							_assert_(pe && ee.GeomEdgeHook);
							ee.GeomEdgeHook->SetMark();
							BamgVertex & v0=ee[0], & v1=ee[1];
							R2 AB=(R2)v1-(R2)v0;
							double L0=L,LAB;
							LAB=LengthInterpole(v0.m,v1.m,AB);
							L+= LAB;

							if (phase){
								// computation of the new points for the given curve
								while ((i!=NbCreatePointOnCurve) && sNew<=L) { 

									//some checks
									_assert_(sNew>=L0);
									_assert_(LAB);
									_assert_(vertices && nbv<maxnbv);
									_assert_(edges && nbe<nbex);
									_assert_(VerticesOnGeomEdge && NbVerticesOnGeomEdge<NbVerticesOnGeomEdgex);

									// new vertex on edge
									A1=vertices+nbv++;
									GA1=VerticesOnGeomEdge+NbVerticesOnGeomEdge;
									Edge* e = edges + nbe++;
									double se= (sNew-L0)/LAB;
									if (se<0 || se>=1.000000001){
										_error_("Problem creating point on a boundary: se=" << se << " should be in [0 1]");
									}
									se = abscisseInterpole(v0.m,v1.m,AB,se,1);
									if (se<0 || se>1){
										_error_("Problem creating point on a boundary: se=" << se << " should be in [0 1]");
									}
									se = k1         ? se : 1. - se;
									se = k1==k1equi ? se : 1. - se;
									VertexOnBThEdge[NbVerticesOnGeomEdge++] = VertexOnEdge(A1,&eeequi,se); // save 
									ongequi=Gh.ProjectOnCurve(eeequi,se,*A1,*GA1); 
									A1->ReferenceNumber = eeequi.ReferenceNumber;
									e->GeomEdgeHook = ongequi;
									e->v[0]=A0;
									e->v[1]=A1;
									e->ReferenceNumber = eeequi.ReferenceNumber;
									e->adj[0]=PreviousNewEdge;

									if (PreviousNewEdge) PreviousNewEdge->adj[1]=e;
									PreviousNewEdge=e;
									A0=A1;
									sNew += Lstep;
									if (++i== NbCreatePointOnCurve) break;
								}
							}

							//some checks
							_assert_(ee.GeomEdgeHook->CurveNumber==ei.GeomEdgeHook->CurveNumber);
							if (ee[k1].GeomEdgeHook->IsRequiredVertex()) {
								_assert_(eeequi[k1equi].GeomEdgeHook->IsRequiredVertex());
								GeomVertex * GA1 = *eeequi[k1equi].GeomEdgeHook;
								A1=GA1->MeshVertexHook;// the vertex in new mesh
								_assert_(A1-vertices>=0 && A1-vertices<nbv);
								break;
							}
							if (!ee.adj[k1]) {
								_error_("adj edge " << BTh.GetId(ee) << ", nbe=" << nbe << ", Gh.vertices=" << Gh.vertices);
							}
							pe = ee.adj[k1]; // next edge
							k0 = pe->Intersection(ee); 
							peequi= eeequi.adj[k1equi];  // next edge
							k0equi=peequi->Intersection(eeequi);            
						}// for(;;) end of the curve
					}

					if (phase){ // construction of the last edge
						Edge* e=edges + nbe++;
						e->GeomEdgeHook  = ongequi;
						e->v[0]=A0;
						e->v[1]=A1;
						e->ReferenceNumber = peequi->ReferenceNumber;
						e->adj[0]=PreviousNewEdge;
						e->adj[1]=0;
						if (PreviousNewEdge) PreviousNewEdge->adj[1]=e;
						PreviousNewEdge = e;

						_assert_(i==NbCreatePointOnCurve);
					}

					if (!phase)  { // 
						long NbSegOnCurve = Max((long)(L+0.5),(long) 1);// nb of seg
						Lstep = L/NbSegOnCurve; 
						NbCreatePointOnCurve = NbSegOnCurve-1;
						NbOfNewEdge += NbSegOnCurve;
						NbOfNewPoints += NbCreatePointOnCurve;
					}
				}
			}//  end of curve loop 

			//Allocate memory
			if(step==0){
				if(nbv+NbOfNewPoints > maxnbv) {
					_error_("too many vertices on geometry: " << nbv+NbOfNewPoints << " >= " << maxnbv);
				}
				edges = new Edge[NbOfNewEdge];
				nbex = NbOfNewEdge;
				if(NbOfNewPoints) {
					VerticesOnGeomEdge    = new VertexOnGeom[NbOfNewPoints];
					NbVertexOnBThEdge     = NbOfNewPoints;
					VertexOnBThEdge       = new  VertexOnEdge[NbOfNewPoints];
					NbVerticesOnGeomEdgex = NbOfNewPoints;
				}
				NbOfNewPoints =0;
				NbOfNewEdge = 0;
			}
		}
		_assert_(nbe!=0);
		delete [] bcurve;

		//Insert points inside existing triangles
		if (verbose>4) _printf_("      -- current number of vertices = " << nbv << "\n");
		if (verbose>3) _printf_("      Creating initial Constrained Delaunay Triangulation...\n");
		if (verbose>3) _printf_("         Inserting boundary points\n");
		Insert(bamgopts);

		//Force the boundary
		if (verbose>3) _printf_("         Forcing boundaries\n");
		ForceBoundary(bamgopts);

		//Extract SubDomains
		if (verbose>3) _printf_("         Extracting subdomains\n");
		FindSubDomain(bamgopts);

		if (verbose>3) _printf_("      Inserting internal points\n");
		NewPoints(BTh,bamgopts,KeepVertices) ;
		if (verbose>4) _printf_("      -- current number of vertices = " << nbv << "\n");
	}
	/*}}}*/
	long  Mesh::RandomNumber(long max){/*{{{*/
		/*  Generate a random number from 0 to max-1 using linear congruential
		 *  random number generator*/

		this->randomseed = (this->randomseed * 1366l + 150889l) % 714025l;
		return this->randomseed / (714025l / max + 1);
	} /*}}}*/
	int Mesh::ForceEdge(BamgVertex &a, BamgVertex & b,AdjacentTriangle & taret)  { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/ForceEdge)*/

		int NbSwap =0;
		if (!a.t || !b.t){ // the 2 vertex is in a mesh
			_error_("!a.t || !b.t");
		}
		int k=0;
		taret=AdjacentTriangle(0,0); // erreur 

		AdjacentTriangle tta(a.t,EdgesVertexTriangle[a.IndexInTriangle][0]);
		BamgVertex   *v1, *v2 = tta.EdgeVertex(0),*vbegin =v2;
		// we turn around a in the  direct direction  

		long long det2 = v2 ? det(*v2,a,b): -1 , det1;
		if(v2) // normal case 
		 det2 = det(*v2,a,b);
		else { // no chance infini vertex try the next
			tta= Previous(Adj(tta));
			v2 = tta.EdgeVertex(0);
			vbegin =v2;
			if (!v2){
				_error_("!v2");
			}
			det2 = det(*v2,a,b);
		}

		while (v2 != &b) {
			AdjacentTriangle tc = Previous(Adj(tta));    
			v1 = v2; 
			v2 = tc.EdgeVertex(0);
			det1 = det2;
			det2 =  v2 ? det(*v2,a,b): det2; 

			if((det1 < 0) && (det2 >0)) { 
				// try to force the edge 
				BamgVertex * va = &a, *vb = &b;
				tc = Previous(tc);
				if (!v1 || !v2){
					_error_("!v1 || !v2");
				}
				long long detss = 0,l=0;
				while ((this->SwapForForcingEdge(  va,  vb, tc, detss, det1,det2,NbSwap)))
				 if(l++ > 10000000) {
					 _error_("Loop in forcing Egde, nb de swap=" << NbSwap << ", nb of try swap (" << l << ") too big");
				 }
				BamgVertex *aa = tc.EdgeVertex(0), *bb = tc.EdgeVertex(1);
				if (((aa == &a ) && (bb == &b)) ||((bb ==  &a ) && (aa == &b))){
					tc.SetLock();
					a.Optim(1,0);
					b.Optim(1,0);
					taret = tc;
					return NbSwap;
				}
				else 
				  {
					taret = tc;
					return -2; // error  boundary is crossing
				  }
			}
			tta = tc;
			k++;
			if (k>=2000){
				_error_("k>=2000");
			}
			if ( vbegin == v2 ) return -1;// error 
		}

		tta.SetLock();
		taret=tta;
		a.Optim(1,0);
		b.Optim(1,0);
		return NbSwap; 
	}
	/*}}}*/
	int Mesh::SwapForForcingEdge(BamgVertex   *  & pva ,BamgVertex  * &   pvb ,AdjacentTriangle & tt1,long long & dets1, long long & detsa,long long & detsb, int & NbSwap) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/SwapForForcingEdge)*/
		// l'arete ta coupe l'arete pva pvb
		// de cas apres le swap sa coupe toujours
		// on cherche l'arete suivante 
		// on suppose que detsa >0 et detsb <0
		// attention la routine echange pva et pvb 

		if(tt1.Locked()) return 0; // frontiere croise 

		AdjacentTriangle tt2 = Adj(tt1);
		Triangle *t1=tt1,*t2=tt2;// les 2 triangles adjacent
		short a1=tt1,a2=tt2;// les 2 numero de l arete dans les 2 triangles
		if ( a1<0 || a1>=3 ){
			_error_("a1<0 || a1>=3");
		}

		BamgVertex & sa= (* t1)[VerticesOfTriangularEdge[a1][0]];
		BamgVertex & s1= (*t1)[OppositeVertex[a1]];
		BamgVertex & s2= (*t2)[OppositeVertex[a2]];

		long long dets2 = det(*pva,*pvb,s2);
		long long det1=t1->det , det2=t2->det ;
		long long detT = det1+det2;
		if ((det1<=0 ) || (det2<=0)){
			_error_("(det1<=0 ) || (det2<=0)");
		}
		if ( (detsa>=0) || (detsb<=0) ){ // [a,b] cut infinite line va,bb
			_error_("(detsa>=0) || (detsb<=0)");
		}
		long long ndet1 = bamg::det(s1,sa,s2);
		long long ndet2 = detT - ndet1;

		int ToSwap =0; //pas de swap
		if ((ndet1 >0) && (ndet2 >0)) 
		  { // on peut swaper  
			if ((dets1 <=0 && dets2 <=0) || (dets2 >=0 && detsb >=0))
			 ToSwap =1; 
			else // swap alleatoire 
			 if (this->RandomNumber(1)) 
			  ToSwap =2; 
		  }
		if (ToSwap) NbSwap++,
		 bamg::swap(t1,a1,t2,a2,&s1,&s2,ndet1,ndet2);

		int ret=1;

		if (dets2 < 0) {// haut
			dets1 = ToSwap ? dets1 : detsa ;
			detsa = dets2; 
			tt1 =  Previous(tt2) ;}
		else if (dets2 > 0){// bas 
			dets1 = ToSwap ? dets1 : detsb ;
			detsb = dets2;
			//xxxx tt1 = ToSwap ? tt1 : Next(tt2);
			if(!ToSwap) tt1 =  Next(tt2);
		}
		else { // changement de direction 
			ret = -1;
			Exchange(pva,pvb);
			Exchange(detsa,detsb);
			Exchange(dets1,dets2);
			Exchange(tt1,tt2);
			dets1=-dets1;
			dets2=-dets2;
			detsa=-detsa;
			detsb=-detsb;

			if(ToSwap){
				if (dets2 < 0) {// haut
					dets1 = (ToSwap ? dets1 : detsa) ;
					detsa = dets2; 
					tt1 =  Previous(tt2) ;}
				else if(dets2 > 0){// bas 
					dets1 = (ToSwap ? dets1 : detsb) ;
					detsb =  dets2;
					if(!ToSwap) tt1 =  Next(tt2);
				}
				else {// on a fin ???
					tt1 = Next(tt2);
					ret =0;}
			}

		}
		return ret;
	}
	/*}}}*/

	/*Intermediary*/
	AdjacentTriangle CloseBoundaryEdge(I2 A,Triangle *t, double &a,double &b) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/CloseBoundaryEdge)*/

		int k=(*t)(0) ?  ((  (*t)(1) ? ( (*t)(2) ? -1 : 2) : 1  )) : 0;
		int dir=0;
		if (k<0){
			_error_("k<0");
		}
		int kkk=0;  
		long long IJ_IA,IJ_AJ;
		AdjacentTriangle edge(t,OppositeEdge[k]);          
		for (;;edge = dir >0 ? Next(Adj(Next(edge))) : Previous(Adj(Previous(edge)))) {  
			kkk++;
			if (kkk>=1000){
				_error_("kkk>=1000");
			}
			BamgVertex  &vI =  *edge.EdgeVertex(0);
			BamgVertex  &vJ =  *edge.EdgeVertex(1);
			I2 I=vI, J=vJ, IJ= J-I;
			IJ_IA = (IJ ,(A-I));
			if (IJ_IA<0) {
				if (dir>0) {a=1;b=0;return edge;}// change of signe => I
				else {dir=-1;
					continue;}};// go in direction i 
					IJ_AJ = (IJ ,(J-A));
					if (IJ_AJ<0) {
						if(dir<0)  {a=0;b=1;return edge;}            
						else {dir = 1;
							continue;}}// go in direction j
							double IJ2 = IJ_IA + IJ_AJ;
							if (IJ2==0){
								_error_("IJ2==0");
							}
							a= IJ_AJ/IJ2;
							b= IJ_IA/IJ2;
							return edge;
		} 
	}
	/*}}}*/
	void  swap(Triangle *t1,short a1, Triangle *t2,short a2, BamgVertex *s1,BamgVertex *s2,long long det1,long long det2){ /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/swap)*/
		// --------------------------------------------------------------
		// short a2=aa[a];// les 2 numero de l arete dans les 2 triangles
		//                               
		//               sb                     sb    
		//             / | \                   /   \                      !
		//         as1/  |  \                 /a2   \                     !
		//           /   |   \               /    t2 \                    !
		//       s1 /t1  | t2 \s2  -->   s1 /___as2___\s2                 !
		//          \  a1|a2  /             \   as1   /  
		//           \   |   /               \ t1    /   
		//            \  |  / as2             \   a1/    
		//             \ | /                   \   /     
		//              sa                       sa   
		//  -------------------------------------------------------------
		int as1 = NextEdge[a1];
		int as2 = NextEdge[a2];
		int ap1 = PreviousEdge[a1];
		int ap2 = PreviousEdge[a2];
		(*t1)(VerticesOfTriangularEdge[a1][1]) = s2 ; // avant sb
		(*t2)(VerticesOfTriangularEdge[a2][1]) = s1  ; // avant sa
		// mise a jour des 2 adjacences externes 
		AdjacentTriangle taas1 = t1->Adj(as1),
							  taas2 = t2->Adj(as2),
							  tas1(t1,as1), tas2(t2,as2),
							  ta1(t1,a1),ta2(t2,a2);
		// externe haut gauche
		taas1.SetAdj2(ta2, taas1.GetAllFlag_UnSwap());
		// externe bas droite
		taas2.SetAdj2(ta1, taas2.GetAllFlag_UnSwap());
		// remove the Mark  UnMarkSwap 
		t1->SetUnMarkUnSwap(ap1);
		t2->SetUnMarkUnSwap(ap2);
		// interne 
		tas1.SetAdj2(tas2);

		t1->det = det1;
		t2->det = det2;

		t1->SetSingleVertexToTriangleConnectivity();
		t2->SetSingleVertexToTriangleConnectivity();
	} // end swap 
	/*}}}*/
}
