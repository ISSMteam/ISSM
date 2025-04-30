/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, R2.h)*/
#ifndef _R2_H
#define _R2_H

#include <cstdio>
#include "../shared/shared.h"

namespace bamg {

	template <class R,class RR> class P2{

		  public:  

			  //fields
			  R x,y;

			  //functions
			  P2 () :x(0),y(0) {};
			  P2 (R a,R b)  :x(a),y(b)  {}
			  P2 (P2 A,P2 B) : x(B.x-A.x),y(B.y-A.y) {}
			  void Echo() const{
				  printf("Member of P2:\n");
				  std::cout<<"   x: "<<x<<std::endl;
				  std::cout<<"   y: "<<y<<std::endl;
			  }
			  //operators
			  RR       operator,(const P2<R,RR> & cc) const {return  (RR) x* (RR) cc.x+(RR) y* (RR) cc.y;} //scalar product
			  P2<R,RR> operator+(const P2<R,RR> & cc) const {return P2<R,RR>(x+cc.x,y+cc.y);}
			  P2<R,RR> operator-(const P2<R,RR> & cc) const {return P2<R,RR>(x-cc.x,y-cc.y);}
			  P2<R,RR> operator-()  const{return P2<R,RR>(-x,-y);}
			  P2<R,RR> operator*(R  cc) const {return P2<R,RR>(x*cc,y*cc);}
			  P2<R,RR> operator/(R  cc) const {return P2<R,RR>(x/cc,y/cc);}
			  P2<R,RR> operator+=(const  P2<R,RR> & cc) {x += cc.x;y += cc.y;return *this;}
			  P2<R,RR> operator/=(const  R r) {x /= r;y /= r;return *this;}
			  P2<R,RR> operator*=(const  R r) {x *= r;y *= r;return *this;}
			  P2<R,RR> operator-=(const  P2<R,RR> & cc) {x -= cc.x;y -= cc.y;return *this;}

	  };

	template <class R,class RR> class P2xP2{

		  public:

			  //fields
			  P2<R,RR> x,y; 

			  //functions
			  P2xP2 (): x(),y()  {}
			  P2xP2 (P2<R,RR> a,P2<R,RR> b): x(a),y(b) {}
			  P2xP2 (P2<R,RR> a,P2<R,RR> b,P2<R,RR> c ): x(b-a),y(c-a) {}
			  P2xP2 (R xx,R xy,R yx,R yy) :x(xx,xy),y(yx,yy) {}
			  void Echo(){
				  printf("Member of P2xP2:\n");
				  printf("   x.x: %g   x.y: %g\n",x.x,x.y);
				  printf("   y.x: %g   y.x: %g\n",y.x,y.y);
			  }
			  RR          det() const {return (RR) x.x* (RR) y.y - (RR) x.y * (RR) y.x;}
			  P2xP2<R,RR> inv()  const{
				  RR d = (*this).det(); 
				  return P2xP2<R,RR>((R)( y.y /d) ,(R)(-x.y/d),(R)( -y.x/d) ,(R)( x.x/d) );
			  };
			  P2xP2<R,RR> t()  {return P2xP2<R,RR>(x.x,y.x,x.y,y.y);} //transposer 
			  //Operators
			  P2<R,RR>     operator*(const P2<R,RR>& c) const {return P2<R,RR>(x.x*c.x + x.y*c.y, y.x*c.x + y.y*c.y);}
			  P2xP2<R,RR>  operator*(P2xP2<R,RR> c) const{
				  return  P2xP2<R,RR>(x.x*c.x.x + x.y*c.y.x,
							  x.x*c.x.y + x.y*c.y.y,
							  y.x*c.x.x + y.y*c.y.x,
							  y.x*c.x.y + y.y*c.y.y);
			  }
	  };  

	//inline functions
	template  <class R,class RR>  
	  inline RR Det(const P2<R,RR> x,const P2<R,RR> y) {
		  return (RR) x.x * (RR) y.y - (RR) x.y * (RR) y.x ;
	  } 
	template  <class R,class RR>  
	  inline RR Area2 (const P2<R,RR> a,const P2<R,RR> b,const P2<R,RR> c) {
		  return Det(b-a,c-a) ;
	  }
	template  <class R,class RR>  
	  inline R Norme1 (const P2<R,RR> x) {
		  return (Abs(x.x)+Abs(x.y)) ;
	  } 
	template  <class R,class RR>  
	  inline RR Norme2_2 (const P2<R,RR> x) {
		  return (RR)x.x*(RR)x.x + (RR)x.y*(RR)x.y ;
	  } 
	template  <class R,class RR>  
	  inline RR Norme2 (const P2<R,RR> x) {
		  return sqrt((RR)x.x*(RR)x.x + (RR)x.y*(RR)x.y) ;
	  } 
	template  <class R,class RR>  
	  inline P2<R,RR> Orthogonal (const P2<R,RR> x) {
		  return  P2<R,RR>(-x.y,x.x);
	  } 
}
#endif
