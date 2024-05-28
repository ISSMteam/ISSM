#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "./bamgobjects.h"
#include "../shared/shared.h"
#include "./det.h"

namespace bamg {

	/*Constructors Destructors*/
	ListofIntersectionTriangles::ListofIntersectionTriangles(int n,int m)/*{{{*/
	  : MaxSize(n), Size(0), len(-1),state(-1),lIntTria(new IntersectionTriangles[n]) ,
	  NbSeg(0), MaxNbSeg(m), lSegsI(new SegInterpolation[m]){
	  }
	/*}}}*/
	ListofIntersectionTriangles::~ListofIntersectionTriangles(){/*{{{*/
		if (lIntTria) delete [] lIntTria,lIntTria=0;
		if (lSegsI) delete [] lSegsI,lSegsI=0;
	} 
	/*}}}*/

	/*Methods*/
	void ListofIntersectionTriangles::Init(void){/*{{{*/
		state=0;
		len=0;
		Size=0;
	}
	/*}}}*/
	double  ListofIntersectionTriangles::Length(){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Length)*/

		// computation of the length

		// check Size
		if (Size<=0){
			_error_("Size<=0");
		}

		Metric Mx,My;
		int ii,jj;
		R2 x,y,xy;

		SegInterpolation* SegI=lSegsI;
		lSegsI[NbSeg].last=Size;
		int EndSeg=Size;

		y = lIntTria[0].x;
		double sxy, s = 0;
		lIntTria[0].s =0;
		SegI->lBegin=s;

		for (jj=0,ii=1;ii<Size;jj=ii++) {  
			// seg jj,ii
			x  = y;
			y  = lIntTria[ii].x;
			xy = y-x;
			Mx = lIntTria[ii].m;
			My = lIntTria[jj].m;
			sxy= LengthInterpole(Mx,My,xy);
			s += sxy;
			lIntTria[ii].s = s;
			if (ii == EndSeg){
				SegI->lEnd=s;
				SegI++;
				EndSeg=SegI->last;
				SegI->lBegin=s;
			}
		}
		len = s;
		SegI->lEnd=s;

		return s;
	}
	/*}}}*/
	int  ListofIntersectionTriangles::NewItem(Triangle * tt,double d0,double d1,double d2) { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/NewItem)*/

		int n;
		R2 x(0,0);
		if ( d0) x =      (*tt)[0].r * d0;
		if ( d1) x = x +  (*tt)[1].r * d1;
		if ( d2) x = x +  (*tt)[2].r * d2;
		// newer add same point 
		if(!Size ||  Norme2_2(lIntTria[Size-1].x-x)) {
			if (Size==MaxSize) ReShape();
			lIntTria[Size].t=tt;
			lIntTria[Size].bary[0]=d0;
			lIntTria[Size].bary[1]=d1;
			lIntTria[Size].bary[2]=d2;
			lIntTria[Size].x = x;
			Metric m0,m1,m2;
			BamgVertex * v;
			if ((v=(*tt)(0))) m0    = v->m;
			if ((v=(*tt)(1))) m1    = v->m;
			if ((v=(*tt)(2))) m2    = v->m;
			lIntTria[Size].m =  Metric(lIntTria[Size].bary,m0,m1,m2);
			n=Size++;}
		else n=Size-1;
		return n;
	}
	/*}}}*/
	int ListofIntersectionTriangles::NewItem(R2 A,const Metric & mm) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/NewItem)*/

		int n;
		if(!Size ||  Norme2_2(lIntTria[Size-1].x-A)) {
			if (Size==MaxSize) ReShape();
			lIntTria[Size].t=0;
			lIntTria[Size].x=A;
			lIntTria[Size].m=mm;
			n=Size++;
		}
		else  n=Size-1;
		return  n; 
	}
	/*}}}*/
	long ListofIntersectionTriangles::NewPoints(BamgVertex* vertices,long &nbv,long maxnbv){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/NewPoints)*/

		//If length<1.5, do nothing
		double s=Length();
		if (s<1.5) return 0;

		const long nbvold=nbv;
		int ii = 1 ;
		R2 y,x;
		Metric My,Mx ;
		double sx =0,sy;
		int nbi=Max(2,(int) (s+0.5));
		double sint=s/nbi;
		double si  =sint;

		int EndSeg=Size;
		SegInterpolation* SegI=NULL;
		if (NbSeg) SegI=lSegsI,EndSeg=SegI->last;

		for (int k=1;k<nbi;k++){
			while ((ii < Size) && ( lIntTria[ii].s <= si )){
				if (ii++ == EndSeg){
					SegI++;
					EndSeg=SegI->last;
				}
			}

			int ii1=ii-1;
			x  =lIntTria[ii1].x;
			sx =lIntTria[ii1].s;
			Metric Mx=lIntTria[ii1].m;
			y  =lIntTria[ii].x;
			sy =lIntTria[ii].s;
			Metric My=lIntTria[ii].m;
			double lxy = sy-sx;
			double cy = abscisseInterpole(Mx,My,y-x,(si-sx)/lxy);

			R2 C;
			double cx = 1-cy;
			C = SegI ? SegI->F(si): x * cx + y *cy;
			//C.Echo();
			//x.Echo();
			//y.Echo();
			//_printf_("cx = " << cx << ", cy=" << cy << "\n");

			si += sint;
			if ( nbv<maxnbv) {
				vertices[nbv].r = C;
				vertices[nbv++].m = Metric(cx,lIntTria[ii-1].m,cy,lIntTria[ii].m);
			}
			else return nbv-nbvold;
		  }
		return nbv-nbvold;
	}
	/*}}}*/
	void ListofIntersectionTriangles::ReShape(){ /*{{{*/

		int newsize = MaxSize*2;
		IntersectionTriangles* nw = new IntersectionTriangles[newsize];
		_assert_(nw);

		// recopy
		for (int i=0;i<MaxSize;i++) nw[i] = lIntTria[i];       
		long int verbosity=0;
		if(verbosity>3) _printf_("   ListofIntersectionTriangles  ReShape Maxsize " << MaxSize << " -> " << MaxNbSeg << "\n");
		MaxSize = newsize; 
		delete [] lIntTria;// remove old
		lIntTria = nw; // copy pointer
	}
	/*}}}*/
	void ListofIntersectionTriangles::SplitEdge(Mesh & Bh, const R2 &A,const R2  &B,int nbegin) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/ListofIntersectionTriangles)*/

		Triangle *tbegin, *t;

		long long deta[3], deti,detj;
		double ba[3];
		int ifirst=-1,ilast;
		int i0,i1,i2;
		int ocut,i,j,k=-1;
		//  int OnAVertices =0;
		long long dt[3];
		I2 a = Bh.R2ToI2(A) ,b= Bh.R2ToI2(B);// compute  the Icoor a,b
		I2 vi,vj;  
		int iedge =-1;// not a edge

		if(nbegin)  {// optimisation 
			// we suppose  knowing the starting  triangle
			t=tbegin=lIntTria[ilast=(Size-1)].t;
			if (tbegin->det>=0) 
			 ifirst = ilast;}  
		else {// not optimisation 
			Init();
			t=tbegin = Bh.TriangleFindFromCoord(a,deta);
			if( t->det>=0)
			 ilast=NewItem(t,double(deta[0])/t->det,double(deta[1])/t->det,double(deta[2])/t->det);
			else 
			  {// find the nearest boundary edge  of the vertex A
				// find a edge or such normal projection a the edge IJ is on the edge
				//   <=> IJ.IA >=0 && IJ.AJ >=0
				ilast=ifirst;
				double ba,bb;
				AdjacentTriangle edge=CloseBoundaryEdge(a,t,ba,bb);
				BamgVertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
				NewItem(A,Metric(ba,v0,bb,v1));
				t=edge;
				// test if the point b is in the same side
				if (det(v0.i,v1.i,b)>=0) {
					AdjacentTriangle edge=CloseBoundaryEdge(a,t,ba,bb);
					BamgVertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
					NewItem(A,Metric(ba,v0,bb,v1));
					return;
				}
			  } // find the nearest boundary edge  of the vertex A
		} // end not optimisation 
		if (t->det<0) {  // outside departure
			while (t->det <0) { // intersection boundary edge and a,b,
				k=(*t)(0) ?  ((  (*t)(1) ? ( (*t)(2) ? -1 : 2) : 1  )) : 0;
				if (k<0){
					_error_("k<0");
				}
				ocut = OppositeEdge[k];
				i=VerticesOfTriangularEdge[ocut][0];
				j=VerticesOfTriangularEdge[ocut][1];
				vi=(*t)[i];
				vj=(*t)[j];
				deti = bamg::det(a,b,vi);
				detj = bamg::det(a,b,vj);
				if (deti>0) // go to  i direction on gamma
				 ocut = PreviousEdge[ocut];      
				else if (detj<=0) // go to j direction on gamma
				 ocut = NextEdge[ocut];         
				AdjacentTriangle tadj =t->Adj(ocut);
				t = tadj;
				iedge= tadj; 
				if (t == tbegin) { // 
					double ba,bb;
					AdjacentTriangle edge=CloseBoundaryEdge(a,t,ba,bb);
					BamgVertex & v0 = *edge.EdgeVertex(0), & v1 = *edge.EdgeVertex(1);
					NewItem(A,Metric(ba,v0,bb,v1));
					return;
				}
			} //  end while (t->det <0)
			// theoriticaly we have: deti =<0 and detj>0

			// computation of barycentric coor
			// test if the point b is on size on t
			// we revert vi,vj because vi,vj is def in Adj triangle
			if ( det(vi,vj,b)>=0) {
				t=tbegin;
				double ba,bb;
				AdjacentTriangle edge=CloseBoundaryEdge(b,t,ba,bb);
				NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
				return;
			}
			else
			  {
				k = OppositeVertex[iedge];
				i=VerticesOfTriangularEdge[iedge][0];
				j=VerticesOfTriangularEdge[iedge][1];
				double dij = detj-deti;
				if (i+j+k != 0 + 1 +2){
					_error_("i+j+k != 0 + 1 +2");
				}
				ba[j] =  detj/dij;
				ba[i] = -deti/dij;
				ba[k] = 0;
				ilast=NewItem(t,ba[0],ba[1],ba[2]); }
		}  //  outside departure

		// recherche the intersection of [a,b] with Bh Mesh.
		// we know  a triangle ta contening the vertex a
		// we have 2 case for intersection [a,b] with a edge [A,B] of Bh
		// 1) the intersection point is in ]A,B[
		// 2)                        is A or B
		// first version --- 
		for (;;) {
			//    t->Draw();
			if (iedge < 0) {
				i0 =0;i1=1;i2=2;
				dt[0] =bamg::det(a,b,(*t)[0]);
				dt[1] =bamg::det(a,b,(*t)[1]);
				dt[2] =bamg::det(a,b,(*t)[2]);}
			else {
				i2 = iedge;
				i0 = NextEdge[i2];
				i1 = NextEdge[i0]; 
				dt[VerticesOfTriangularEdge[iedge][0]] = detj;// we revert i,j because
				dt[VerticesOfTriangularEdge[iedge][1]] = deti;// we take the Triangle by the other side
				dt[iedge] = det(a,b,(*t)[OppositeVertex[iedge]]);}

				// so we have just to see the transition from - to + of the det0..2 on edge of t
				// because we are going from a to b
				if       ((dt[i=VerticesOfTriangularEdge[i0][0]] <  0) &&
							( dt[j=VerticesOfTriangularEdge[i0][1]] > 0))
				 ocut =i0;
				else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] <  0) &&
							(dt[j=VerticesOfTriangularEdge[i1][1]] >  0))
				 ocut =i1;
				else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] <  0) && 
							(dt[j=VerticesOfTriangularEdge[i2][1]] >  0))
				 ocut =i2;
				else if   ((dt[i=VerticesOfTriangularEdge[i0][0]] == 0) &&
							( dt[j=VerticesOfTriangularEdge[i0][1]] >  0))
				 ocut =i0;
				else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] == 0) &&
							(dt[j=VerticesOfTriangularEdge[i1][1]] >  0))
				 ocut =i1;
				else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] == 0) && 
							(dt[j=VerticesOfTriangularEdge[i2][1]] >  0))
				 ocut =i2;
				else if   ((dt[i=VerticesOfTriangularEdge[i0][0]] <  0) &&
							( dt[j=VerticesOfTriangularEdge[i0][1]] == 0))
				 ocut =i0;
				else  if ((dt[i=VerticesOfTriangularEdge[i1][0]] <  0) &&
							(dt[j=VerticesOfTriangularEdge[i1][1]] == 0))
				 ocut =i1;
				else  if ((dt[i=VerticesOfTriangularEdge[i2][0]] <  0) && 
							(dt[j=VerticesOfTriangularEdge[i2][1]] == 0))
				 ocut =i2;
				else { //  On a edge (2 zero)
					k =0;
					if (dt[0]) ocut=0,k++; 
					if (dt[1]) ocut=1,k++; 
					if (dt[2]) ocut=2,k++;
					if(k == 1) {
						if (dt[ocut] >0) // triangle upper AB
						 ocut = NextEdge[ocut];
						i= VerticesOfTriangularEdge[ocut][0];
						j= VerticesOfTriangularEdge[ocut][1];
					}
					else {
						_error_("Bug Split Edge");
					}
				}

				k = OppositeVertex[ocut];

				long long detbij = bamg::det((*t)[i],(*t)[j],b);

				if (detbij >= 0) { //we find the triangle contening b
					dt[0]=bamg::det((*t)[1],(*t)[2],b);
					dt[1]=bamg::det((*t)[2],(*t)[0],b);
					dt[2]=bamg::det((*t)[0],(*t)[1],b);
					double dd = t->det;
					NewItem(t,dt[0]/dd,dt[1]/dd,dt[2]/dd);      
					return ;}
				else { // next triangle by  adjacent by edge ocut 
					deti = dt[i];
					detj = dt[j];
					double dij = detj-deti;
					ba[i] =  detj/dij;
					ba[j] = -deti/dij;
					ba[3-i-j ] = 0;
					ilast=NewItem(t, ba[0],ba[1],ba[2]);      

					AdjacentTriangle ta =t->Adj(ocut);
					t = ta;
					iedge= ta; 
					if (t->det <= 0)  {
						double ba,bb;
						AdjacentTriangle edge=CloseBoundaryEdge(b,t,ba,bb);
						NewItem(B,Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1)));
						return;
					}
				}// we  go outside of omega 
		} // for(;;)
	}
	/*}}}*/

}
