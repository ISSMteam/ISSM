#include <cstdio>
#include <string.h>
#include <cmath>

#include "../shared/shared.h"

#include "GeomEdge.h"
#include "Geometry.h"

using namespace std;

namespace bamg {

	/*Constructor/Destructor*/

	/*Methods*/
	int    GeomEdge::Cracked() const  {/*{{{*/
		return type &1;  
	}/*}}}*/
	R2 GeomEdge::F(double theta) const{/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, MeshGeom.cpp/F)*/
		// parametrization of the curve edge

	   R2 A=v[0]->r,B=v[1]->r;
		double ca,cb,cta,ctb;

		//Check that theta is in [0 1]
		_assert_(theta>-1e-12 && theta<1+1e-12);

		if (TgA()){ 
			if (TgB()){ //Hermite interpolation
				cb =  theta*theta*(3-2*theta);
				ca =  1-cb;     
				cta = (1-theta)*(1-theta)*theta;
				ctb = (theta-1)*theta*theta ;
			}
			else {
				double t = theta;
				cb = t*t;
				ca = 1-cb;
				cta= t-cb;
				ctb=0;    
			}
		}
		else{
			if (TgB()){
				double t = 1-theta;
				ca = t*t;
				cb = 1-ca;
				ctb= -t+ca;
				cta=0;    
			}
			else { // lagrange P1
				ca =(1-theta);
				cb = theta;
				cta=ctb=0;
			}
		}
		return A*ca + B*cb + tg[0]*cta + tg[1]*ctb;
	  }
	/*}}}*/
	int    GeomEdge::Mark()    const  {/*{{{*/
		return type &16; 
	}/*}}}*/
	int    GeomEdge::Required()       {/*{{{*/
		return type &64; 
	}/*}}}*/
	void GeomEdge::Set(const GeomEdge & rec,const Geometry & Gh ,Geometry & GhNew){ /*{{{*/
		*this = rec;
		v[0] = GhNew.vertices + Gh.GetId(v[0]);    
		v[1] = GhNew.vertices + Gh.GetId(v[1]); 
		if (Adj[0]) Adj[0] =  GhNew.edges + Gh.GetId(Adj[0]);     
		if (Adj[1]) Adj[1] =  GhNew.edges + Gh.GetId(Adj[1]);     
	}
	/*}}}*/
	void   GeomEdge::SetCracked()     { /*{{{*/
		type |= 1;/*=>1st digit to 1*/
	}/*}}}*/
	void   GeomEdge::SetTgA()         { /*{{{*/
		type |=4; /*=>2d digit to 1*/
	}/*}}}*/
	void   GeomEdge::SetTgB()         { /*{{{*/
		type |=8; /*=> 3d digit to 1*/
	}/*}}}*/
	void   GeomEdge::SetMark()        { /*{{{*/
		type |=16;/*=> 4th digiy to 1*/
	}/*}}}*/
	void   GeomEdge::SetUnMark()      { /*{{{*/
		type &= 1007 /* 1023-16 = 000111110111 => 4th digit to 0*/;
	}/*}}}*/
	void   GeomEdge::SetRequired()    { /*{{{*/
		type |= 64;/*=>6th digit to 1*/ 
	}/*}}}*/
	int    GeomEdge::TgA()     const  {/*{{{*/
		return type &4;  
	}/*}}}*/
	int    GeomEdge::TgB()     const  {/*{{{*/
		return type &8;  
	}/*}}}*/
}
