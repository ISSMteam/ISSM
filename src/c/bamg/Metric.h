#ifndef _METRIC_H
#define _METRIC_H

#include "./include.h"
#include "../shared/shared.h"
#include "R2.h"
#include <math.h>

namespace bamg {

	typedef P2<double,double>    D2;
	typedef P2xP2<double,double> D2xD2;

	class Metric;
	class EigenMetric;

	class Metric{

		public:

			//fields
			double a11,a21,a22;

			//friends
			friend class EigenMetric;

			//functions
			Metric():a11(0),a21(0),a22(0){};
			Metric(const EigenMetric&);
			Metric(double a);
			Metric(double a,double b,double c);
			Metric(double  a,const Metric ma,double  b,const Metric mb);
			Metric(const double  a[3],const Metric m0,const Metric m1,const Metric m2 );
			void        Echo();
			double      det() const;
			int         IntersectWith(const  Metric& M2);
			inline void Box(double &hx,double &hy) const;

			/*The following functions must remain the the header file because it is called before Metric
			 * is compiled by other classes*/
			R2 Orthogonal(const R2 x){ return R2(-(a21*x.x+a22*x.y),a11*x.x+a21*x.y); }
			R2 Orthogonal(const I2 x){ return R2(-(a21*x.x+a22*x.y),a11*x.x+a21*x.y); }
			double Length(double Ax,double Ay) const;

			//operators
			Metric operator*(double c) const {double c2=c*c;return  Metric(a11*c2,a21*c2,a22*c2);} 
			Metric operator/(double c) const {double c2=1/(c*c);return  Metric(a11*c2,a21*c2,a22*c2);} 
			operator D2xD2(){ return D2xD2(a11,a21,a21,a22);}
			//double  operator()(R2 x) const { return sqrt(x.x*x.x*a11+2*x.x*x.y*a21+x.y*x.y*a22);};        // length of x in metric sqrt(<Mx,x>) FIXME: replace by Length!
			double  operator()(R2 x,R2 y) const { return x.x*y.x*a11+(x.x*x.y+x.y*y.x)*a21+x.y*y.y*a22;};

	};

	class EigenMetric{
		public:

			//fields
			double lambda1,lambda2;
			double vx,vy;

			//friends
			friend  class Metric;

			//functions
			EigenMetric(const Metric& );
			EigenMetric(double r1,double r2,const D2& vp1);
			void   Echo();
			void   Abs();
			void   pow(double  p);
			void   Min(double  a);
			void   Max(double  a);
			void   Minh(double h);
			void   Maxh(double h);
			double hmin()   const;
			double hmax()   const;
			double lmax()   const;
			double lmin()   const;
			double Aniso2() const;
			inline void BoundAniso2(const double coef);

			//operators
			void operator *=(double coef){ lambda1*=coef;lambda2*=coef;}
	};

	class SaveMetricInterpole {
		friend double LengthInterpole(const Metric& Ma,const  Metric& Mb, R2 AB);
		friend double abscisseInterpole(const Metric& Ma ,const  Metric& Mb, R2 ,double s,int optim);
		public:
		int opt;
		double lab;
		double L[1024],S[1024];
	};

	extern SaveMetricInterpole  LastMetricInterpole; // for optimization 
	//Functions
	void  SimultaneousMatrixReduction( Metric M1,  Metric M2,D2xD2 &V);
	double LengthInterpole(const Metric& Ma,const  Metric& Mb, R2 AB);
	double abscisseInterpole(const Metric& Ma,const  Metric& Mb, R2 AB,double s,int optim=0);

	//inlines
	inline void  EigenMetric::BoundAniso2(const double coef){
		if (coef<=1.00000000001){
			if (lambda1 < lambda2)
			 lambda1 = bamg::Max(lambda1,lambda2*coef);
			else
			 lambda2 = bamg::Max(lambda2,lambda1*coef);
		}
		else{  //TO BE CHECKED
			if (lambda1 > lambda2)
			 lambda1 = bamg::Min(lambda1,lambda2*coef);
			else
			 lambda2 = bamg::Min(lambda2,lambda1*coef);
		}
	}
	inline Metric::Metric(const EigenMetric& M) {
		double v00=M.vx*M.vx;
		double v11=M.vy*M.vy;
		double v01=M.vx*M.vy;
		a11=v00*M.lambda1+v11*M.lambda2;
		a21=v01*(M.lambda1-M.lambda2);
		a22=v00*M.lambda2+v11*M.lambda1;
	}
	inline   void  Metric::Box(double &hx,double &hy) const {
		double d=  a11*a22-a21*a21;
		hx = sqrt(a22/d);
		hy = sqrt(a11/d);
	}
	inline double LengthInterpole(double la,double lb) {
		return ( Abs(la - lb) < 1.0e-6*Max3(la,lb,1.0e-20) ) ?  (la+lb)/2  : la*lb*log(la/lb)/(la-lb);
	}
	inline double abscisseInterpole(double la,double lb,double lab,double s){
		return ( Abs(la - lb) <1.0e-6*Max3(la,lb,1.0e-20))  ? s : (exp(s*lab*(la-lb)/(la*lb))-1)*lb/(la-lb);
	}
}
#endif
