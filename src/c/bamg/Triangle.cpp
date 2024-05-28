#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "./bamgobjects.h"
#include "../shared/shared.h"
#include "./det.h"

namespace bamg {

	/*Constructors/Destructors*/
	Triangle::Triangle(void){/*{{{*/

	}
	/*}}}*/
	Triangle::Triangle(Mesh *Th,long i,long j,long k) {/*{{{*/
		BamgVertex *v=Th->vertices;
		long nbv = Th->nbv;
		if (i<0 || j<0 || k<0){
			_error_("i<0 || j<0 || k<0");
		}
		if (i>=nbv || j>=nbv || k>=nbv){
			_error_("i>=nbv || j>=nbv || k>=nbv");
		}
		vertices[0]=v+i;
		vertices[1]=v+j;
		vertices[2]=v+k;
		adj[0]=adj[1]=adj[2]=0;
		AdjEdgeIndex[0]=AdjEdgeIndex[1]=AdjEdgeIndex[2]=0;
		det=0;
	}
	/*}}}*/
	Triangle::Triangle(BamgVertex *v0,BamgVertex *v1,BamgVertex *v2){/*{{{*/
		vertices[0]=v0;
		vertices[1]=v1;
		vertices[2]=v2;
		adj[0]=adj[1]=adj[2]=0;
		AdjEdgeIndex[0]=AdjEdgeIndex[1]=AdjEdgeIndex[2]=0;
		if (v0) det=0;
		else {
			det=-1;
			link=NULL;};  
	}
	/*}}}*/

	/*Methods*/
	AdjacentTriangle Triangle::Adj(int i)  const {/*{{{*/
		return AdjacentTriangle(adj[i],AdjEdgeIndex[i]&3);
	};/*}}}*/
	double Triangle::Length() const{/*{{{*/

		double l;

		/*Get three vertices A,B and C*/
		R2 A=*this->vertices[0];
		R2 B=*this->vertices[1];
		R2 C=*this->vertices[2];

		/*Compute edges*/
		R2 e1=B-A;
		R2 e2=C-A;
		R2 e3=B-C;

		/*Compute edge length*/
		l=Norme2(e1);
		l=max(l,Norme2(e2));
		l=max(l,Norme2(e3));

		return l;
	};/*}}}*/
	void Triangle::Echo(void){/*{{{*/

		int i;

		_printf_("Triangle:\n");
		_printf_("   vertices pointer towards three vertices\n");
		_printf_("      vertices[0] vertices[1] vertices[2] = " << vertices[0] << " " << vertices[1] << " " << vertices[2] << "\n");
		_printf_("   adj pointer towards three adjacent triangles\n");
		_printf_("      adj[0] adj[1] adj[2] = " << adj[0] << " " << adj[1] << " " << adj[2] << "\n");
		_printf_("   det (integer triangle determinant) = " << det << "\n");
		if (link){
			_printf_("   link (pointer toward duplicate triangle)= " << link << "\n");
		}
		else{
			_printf_("   color = " << color << "\n");
		}

		_printf_("\nThree vertices:\n");
		for(i=0;i<3;i++){
			if (vertices[i]){
				vertices[i]->Echo();
			}
			else{
				_printf_("   vertex " << i+1 << " does not exist\n");
			}
		}

		return;
	}
	/*}}}*/
	int    Triangle::GetAllflag(int a){/*{{{*/
		return AdjEdgeIndex[a] & 1020;
	}/*}}}*/
	int    Triangle::Hidden(int a)const {/*{{{*/
		return AdjEdgeIndex[a]&16;
	} /*}}}*/
	int    Triangle::Locked(int a)const {/*{{{*/
		return AdjEdgeIndex[a]&4;
	} /*}}}*/
	short  Triangle::NuEdgeTriangleAdj(int i) const {/*{{{*/
		/*Number of the  adjacent edge in adj tria (make sure it is between 0 and 2*/
		return AdjEdgeIndex[i&3]&3;
	}/*}}}*/
	long  Triangle::Optim(short i,int koption) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Optim)*/

		// turn around (positive direction)
		Triangle *t=this;
		long NbSwap =0;
		int  k = 0;
		int  j = OppositeEdge[i];
		int  jp= PreviousEdge[j];

		// initialize tp, jp the previous triangle & edge
		Triangle *tp=adj[jp];
		jp = AdjEdgeIndex[jp]&3;
		do {
			while (t->swap(j,koption)){
				if (k>=20000) _error_("k>=20000");
				NbSwap++;
				k++;
				t=  tp->adj[jp];      // set unchange t qnd j for previous triangles
				j=  NextEdge[tp->AdjEdgeIndex[jp]&3];
			}
			// end on this  Triangle 
			tp = t;
			jp = NextEdge[j];

			t=  tp->adj[jp];      // set unchange t qnd j for previous triangles
			j=  NextEdge[tp->AdjEdgeIndex[jp]&3];

		} while( t != this);
		return NbSwap;
	}
	/*}}}*/
	void  Triangle::Renumbering(Triangle *tb,Triangle *te, long *renu){/*{{{*/

		if (link  >=tb && link  <te) link  = tb + renu[link -tb];
		if (adj[0] >=tb && adj[0] <te) adj[0] = tb + renu[adj[0]-tb];
		if (adj[1] >=tb && adj[1] <te) adj[1] = tb + renu[adj[1]-tb];
		if (adj[2] >=tb && adj[2] <te) adj[2] = tb + renu[adj[2]-tb];    
	}/*}}}*/
	void Triangle::Renumbering(BamgVertex *vb,BamgVertex *ve, long *renu){/*{{{*/
		if (vertices[0] >=vb && vertices[0] <ve) vertices[0] = vb + renu[vertices[0]-vb];
		if (vertices[1] >=vb && vertices[1] <ve) vertices[1] = vb + renu[vertices[1]-vb];
		if (vertices[2] >=vb && vertices[2] <ve) vertices[2] = vb + renu[vertices[2]-vb];    
	}/*}}}*/
	void Triangle::Set(const Triangle & rec,const Mesh & Th ,Mesh & ThNew){ /*{{{*/
		*this = rec;
		if ( vertices[0] ) vertices[0] = ThNew.vertices +  Th.GetId(vertices[0]);
		if ( vertices[1] ) vertices[1] = ThNew.vertices +  Th.GetId(vertices[1]);
		if ( vertices[2] ) vertices[2] = ThNew.vertices +  Th.GetId(vertices[2]);
		if(adj[0]) adj[0] =  ThNew.triangles + Th.GetId(adj[0]);
		if(adj[1]) adj[1] =  ThNew.triangles + Th.GetId(adj[1]);
		if(adj[2]) adj[2] =  ThNew.triangles + Th.GetId(adj[2]);
		if (link  >= Th.triangles && link  < Th.triangles + Th.nbt)
		 link = ThNew.triangles + Th.GetId(link);
	}
	/*}}}*/
	void Triangle::SetAdjAdj(short a){/*{{{*/
		// Copy all the mark 
		a &= 3;
		Triangle *tt=adj[a];
		AdjEdgeIndex [a] &= 55; // remove MarkUnSwap
		short aatt = AdjEdgeIndex[a] & 3;
		if(tt){ 
			tt->adj[aatt]=this;
			tt->AdjEdgeIndex[aatt]=a + (AdjEdgeIndex[a] & 60 ) ;
		}
	}/*}}}*/
	void Triangle::SetAdj2(short a,Triangle *t,short aat){/*{{{*/
		/*For current triangle:
		 * - a is the index of the edge were the adjency is set (in [0 2])
		 * - t is the adjacent triangle
		 * - aat is the index of the same edge in the adjacent triangle*/
		adj[a]=t;
		AdjEdgeIndex[a]=aat;
		if(t){ //if t!=NULL add adjacent triangle to t (this)
			t->adj[aat]=this;
			t->AdjEdgeIndex[aat]=a;
		}
	}/*}}}*/
	void Triangle::SetHidden(int a){/*{{{*/
		//Get Adjacent Triangle number a
		Triangle* t = adj[a];
		//if it exist
		//C|=D -> C=(C|D) bitwise inclusive OR
		if(t) t->AdjEdgeIndex[AdjEdgeIndex[a] & 3] |=16;
		AdjEdgeIndex[a] |= 16;
	}/*}}}*/
	void Triangle::SetLocked(int a){/*{{{*/
		//mark the edge as on Boundary
		Triangle * t = adj[a];
		t->AdjEdgeIndex[AdjEdgeIndex[a] & 3] |=4;
		AdjEdgeIndex[a] |= 4;
	}/*}}}*/
	void Triangle::SetMarkUnSwap(int a){/*{{{*/
		Triangle * t = adj[a];
		t->AdjEdgeIndex[AdjEdgeIndex[a] & 3] |=8;
		AdjEdgeIndex[a] |=8 ;
	}/*}}}*/
	void Triangle::SetSingleVertexToTriangleConnectivity() { /*{{{*/
		if (vertices[0]) (vertices[0]->t=this,vertices[0]->IndexInTriangle=0);
		if (vertices[1]) (vertices[1]->t=this,vertices[1]->IndexInTriangle=1);
		if (vertices[2]) (vertices[2]->t=this,vertices[2]->IndexInTriangle=2);
	}/*}}}*/
	void Triangle::SetUnMarkUnSwap(int a){ /*{{{*/
		Triangle * t = adj[a];
		t->AdjEdgeIndex[AdjEdgeIndex[a] & 3] &=55; // 23 + 32 
		AdjEdgeIndex[a] &=55 ;
	}/*}}}*/
	Triangle* Triangle::TriangleAdj(int i) const {/*{{{*/
		return adj[i&3];
	}/*}}}*/
	int Triangle::swap(short a,int koption){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/swap)*/

		if(a/4 !=0) return 0;// arete lock or MarkUnSwap

		Triangle *t1=this,*t2=adj[a];// les 2 triangles adjacent
		short a1=a,a2=AdjEdgeIndex[a];// les 2 numero de l arete dans les 2 triangles
		if(a2/4 !=0) return 0; // arete lock or MarkUnSwap

		BamgVertex  *sa=t1->vertices[VerticesOfTriangularEdge[a1][0]];
		BamgVertex  *sb=t1->vertices[VerticesOfTriangularEdge[a1][1]];
		BamgVertex  *s1=t1->vertices[OppositeVertex[a1]];
		BamgVertex  *s2=t2->vertices[OppositeVertex[a2]];

		long long det1=t1->det , det2=t2->det ;
		long long detT = det1+det2;
		long long detA = Abs(det1) + Abs(det2);
		long long detMin = Min(det1,det2);

		int OnSwap = 0;       
		// si 2 triangle infini (bord) => detT = -2;
		if (sa == 0) {// les deux triangles sont frontieres
			det2=bamg::det(s2->i,sb->i,s1->i);
			OnSwap = det2 >0;}
		else if (sb == 0) { // les deux triangles sont frontieres
			det1=bamg::det(s1->i,sa->i,s2->i);
			OnSwap = det1 >0;}
		else if(( s1 != 0) && (s2 != 0) ) {
			det1 = bamg::det(s1->i,sa->i,s2->i);
			det2 = detT - det1;
			OnSwap = (Abs(det1) + Abs(det2)) < detA;

			long long detMinNew=Min(det1,det2);
			if (! OnSwap &&(detMinNew>0)) {
				OnSwap = detMin ==0;
				if (! OnSwap) {
					int  kopt = koption;
					while (1)
					 if(kopt) {
						 // critere de Delaunay pure isotrope
						 long long xb1 = sb->i.x - s1->i.x,
								  x21 = s2->i.x - s1->i.x,
								  yb1 = sb->i.y - s1->i.y,
								  y21 = s2->i.y - s1->i.y,
								  xba = sb->i.x - sa->i.x, 
								  x2a = s2->i.x - sa->i.x,
								  yba = sb->i.y - sa->i.y,
								  y2a = s2->i.y - sa->i.y;
						 double
							cosb12 =  double(xb1*x21 + yb1*y21),
									 cosba2 =  double(xba*x2a + yba*y2a) ,
									 sinb12 = double(det2),
									 sinba2 = double(t2->det);

						 // angle b12 > angle ba2 => cotg(angle b12) < cotg(angle ba2)
						 OnSwap =  ((double) cosb12 * (double)  sinba2) <  ((double) cosba2 * (double) sinb12);
						 break;
					 }
					 else {	
						 // critere de Delaunay anisotrope 
						 double som;
						 I2 AB=(I2) *sb - (I2) *sa;
						 I2 MAB2=((I2) *sb + (I2) *sa);
						 R2 MAB(MAB2.x*0.5,MAB2.y*0.5);
						 I2 A1=(I2) *s1 - (I2) *sa;
						 I2 D = (I2) * s1 - (I2) * sb ;
						 R2 S2(s2->i.x,s2->i.y);
						 R2 S1(s1->i.x,s1->i.y);
							{
							 Metric M=s1->m;
							 R2 ABo = M.Orthogonal(AB);
							 R2 A1o = M.Orthogonal(A1);
							 // (A+B)+ x ABo = (S1+B)/2+ y A1 
							 // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2
							 double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
							 double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
							 if (Abs(d) > dd*1.e-3) {
								 R2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
								 som  = M.Length(C.x-S2.x,C.y-S2.y) / M.Length(C.x-S1.x,C.y-S1.y);
							 } else 
								{kopt=1;continue;}

							}
							{
							 Metric M=s2->m;
							 R2 ABo = M.Orthogonal(AB);
							 R2 A1o = M.Orthogonal(A1);
							 // (A+B)+ x ABo = (S1+B)/2+ y A1 
							 // ABo x - A1o y =  (S1+B)/2-(A+B)/2 = (S1-B)/2 = D/2 
							 double dd = Abs(ABo.x*A1o.y)+Abs(ABo.y*A1o.x);
							 double d = (ABo.x*A1o.y - ABo.y*A1o.x)*2; // because D/2
							 if(Abs(d) > dd*1.e-3) {
								 R2 C(MAB+ABo*((D.x*A1o.y - D.y*A1o.x)/d));
								 som += M.Length(C.x-S2.x,C.y-S2.y) / M.Length(C.x-S1.x,C.y-S1.y);
							 } else 
								{kopt=1;continue;}
							}
						 OnSwap = som < 2;
						 break;
					 }

				} // OnSwap 
			} // (! OnSwap &&(det1 > 0) && (det2 > 0) )
		}
		if( OnSwap ) 
		 bamg::swap(t1,a1,t2,a2,s1,s2,det1,det2);
		else {
			t1->SetMarkUnSwap(a1);     
		}
		return OnSwap;
	}
	/*}}}*/

}
