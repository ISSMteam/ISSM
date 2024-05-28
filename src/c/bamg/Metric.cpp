#include <cstdio>
#include <string.h>
#include <cmath>

#include "Metric.h"
#include "../shared/shared.h"

using namespace std;

namespace bamg {

	SaveMetricInterpole  LastMetricInterpole;

	/*Constructor/Destructor*/
	Metric::Metric(double a){/*{{{*/

		/*Isotropic metric for a length of a as unit*/
		this->a11 = 1./(a*a);
		this->a21 = 0.;
		this->a22 = 1./(a*a);

	}/*}}}*/
	Metric::Metric(double a11_in,double a21_in,double a22_in){/*{{{*/

		this->a11 = a11_in;
		this->a21 = a21_in;
		this->a22 = a22_in;

	}/*}}}*/
	Metric::Metric(const double  a[3],const Metric m0,const Metric m1,const Metric m2 ){/*{{{*/

		/*Create metric using 3 inputs*/
		Metric mab(a[0]*m0.a11 + a[1]*m1.a11 + a[2]*m2.a11,
					a[0]*m0.a21 + a[1]*m1.a21 + a[2]*m2.a21,
					a[0]*m0.a22 + a[1]*m1.a22 + a[2]*m2.a22);

		/*Convert to eigen metric*/
		EigenMetric vab(mab);
		double v1x = + vab.vx;
		double v1y = + vab.vy;
		double v2x = - vab.vy;
		double v2y = + vab.vx;

		double h1=a[0] / m0.Length(v1x,v1y) + a[1] / m1.Length(v1x,v1y) + a[2] / m2.Length(v1x,v1y);
		double h2=a[0] / m0.Length(v2x,v2y) + a[1] / m1.Length(v2x,v2y) + a[2] / m2.Length(v2x,v2y);

		vab.lambda1 =  1. / (h1*h1);
		vab.lambda2 =  1. / (h2*h2);

		/*Build metric from vab*/
		double v00=vab.vx*vab.vx;
		double v11=vab.vy*vab.vy;
		double v01=vab.vx*vab.vy;
		this->a11=v00*vab.lambda1+v11*vab.lambda2;
		this->a21=v01*(vab.lambda1-vab.lambda2);
		this->a22=v00*vab.lambda2+v11*vab.lambda1;

	} /*}}}*/
	Metric::Metric(double  a,const Metric ma,double  b,const Metric mb) { /*{{{*/

		/*Compute metric (linear combination of ma and mb)*/
		Metric mab(a*ma.a11+b*mb.a11,a*ma.a21+b*mb.a21,a*ma.a22+b*mb.a22);

		/*Get Eigen values and vectors*/
		EigenMetric vab(mab);
		double v1x = + vab.vx;
		double v1y = + vab.vy;
		double v2x = - vab.vy;
		double v2y = + vab.vx;

		/*Modify eigen values (a+b=1)*/
		double h1 = a/ma.Length(v1x,v1y) + b/mb.Length(v1x,v1y);
		double h2 = a/ma.Length(v2x,v2y) + b/mb.Length(v2x,v2y);
		vab.lambda1 =  1./(h1*h1);
		vab.lambda2 =  1./(h2*h2);

		/*Build metric from vab*/
		double v00=vab.vx*vab.vx;
		double v11=vab.vy*vab.vy;
		double v01=vab.vx*vab.vy;
		this->a11=v00*vab.lambda1+v11*vab.lambda2;
		this->a21=v01*(vab.lambda1-vab.lambda2);
		this->a22=v00*vab.lambda2+v11*vab.lambda1;
	}
	/*}}}*/

	/*Methods*/
	double Metric::det() const {/*{{{*/
		return a11*a22-a21*a21;
	}  /*}}}*/
	void Metric::Echo(void){/*{{{*/

		_printf_("Metric:\n");
		_printf_("   [a11 a21 a22]: [" << a11 << " " << a21 << " " << a22 << "]\n");

		return;
	}
	/*}}}*/
	int Metric::IntersectWith(const Metric& M2) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/IntersectWith)*/

		/*Get a new metric from an existing metric (M1=this)
		 * and a new metric given in input M2 using a 
		 * Simultaneous Matrix Reduction:
		 * If M1 and M2 are 2 metrics, we must build N=M1^-1 M2 (Alauzet2003 p16) 
		 * the eigen vectors of N form a matrix P
		 * The new metric M = M1 inter M2 is then given by:
		 *
		 *      -T [ max(lambda1, mu1)          0         ]  -1				 
		 * M = P   [                                      ] P		 
		 *         [        0            max(lambda2, mu2)] 
		 *
		 * where lambdai and mui can be computed using Rayleigh formula: 
		 *    lambdai = vi' M1 vi
		 * with vi eigen vectors of N (columns of P)
		 */

		int         change=0;
		Metric &M1=*this;
		D2xD2       P;

		//Get P, eigen vectors of N=inv(M1) M2
		SimultaneousMatrixReduction(*this,M2,P);

		//extract the eigen vectors of P (columns)
		R2 v1(P.x.x,P.y.x);
		R2 v2(P.x.y,P.y.y);

		//compute lambdai mui
		double lambda1=M1(v1,v1);
		double lambda2=M1(v2,v2);
		double mu1=M2(v1,v1);
		double mu2=M2(v2,v2);

		//check where any change needs to be done on M1
		if ( lambda1 < mu1 )  change=1,lambda1=mu1;
		if ( lambda2 < mu2 )  change=1,lambda2=mu2; 

		//update M1 if necessary
		if (change) {
			D2xD2 invP(P.inv());
			D2xD2 D(lambda1,0,0,lambda2); 
			D2xD2 M(invP.t()*D*invP);
			a11=M.x.x;
			a21=0.5*(M.x.y+M.y.x);
			a22=M.y.y;
		}
		return change;
	}
	/*}}}*/
	double Metric::Length(double Ax,double Ay) const{/*{{{*/
		/*Length of A in the current metric*/
		return sqrt(Ax*Ax*a11+2*Ax*Ay*a21+Ay*Ay*a22);
	}
	/*}}}*/

	/*Intermediary*/
	double LengthInterpole(const Metric& Ma,const  Metric& Mb, R2 AB) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/LengthInterpole)*/

		double k=1./2.;
		int level=0;
		static int kkk=0;
		static  Metric Ms1[32],Ms2[32];
		static double lMs1[32],lMs2[32];
		static double K[32];
		double l=0,sss=0;
		Ms1[level]=Ma;
		Ms2[level]=Mb;
		double sa =  Ma.Length(AB.x,AB.y);
		double sb =  Mb.Length(AB.x,AB.y);
		lMs1[level]=sa;
		lMs2[level]=sb;
		K[level]=k;
		level++;
		int i=0;
		double * L= LastMetricInterpole.L, *S = LastMetricInterpole.S;
		double  sstop = 0.1; // Max(0.6,(sa+sb)/5000);
		while (level) {
			level--;
			Metric M1=Ms1[level];
			Metric M2=Ms2[level];
			k=K[level];
			double s1=  lMs1[level];
			double s2=  lMs2[level];

			double s= (s1+s2)*k;
			if( s > sstop   && level < 30 && i < 500-level ) {
				Metric Mi(0.5,M1,0.5,M2);
				double si = Mi.Length(AB.x,AB.y);
				if( Abs((s1+s2)-(si+si)) > s1*0.001) 
				  {
					k=k/2;
					// we begin by the end to walk in the correct direction from a to b
					// due to the stack 
					Ms1[level]=Mi;
					Ms2[level]=M2;
					lMs1[level]=si;
					lMs2[level]=s2;
					K[level]=k;
					level++;
					Ms1[level]=M1;
					Ms2[level]=Mi;
					lMs1[level]=s1;
					lMs2[level]=si;
					K[level]=k;
					level++;
				  }
				else
				 L[i]= l += s,S[i]=sss+=k,i++;
			}
			else 
			 L[i]= l += s,S[i]=sss+=k,i++;
		}
		// warning for optimisation S is in [0:0.5] not in [0:1]
		if (i>=512){
			_error_("i>=512");
		}
		LastMetricInterpole.lab=l;
		LastMetricInterpole.opt=i;
		if (i>200 && kkk++<10) _printf_("WARNING: LengthInterpole: ( i=" << i << " l=" << l << " sss=" << sss << " ) " << sstop << "\n"); 
		return l;
	}
	/*}}}*/
	void SimultaneousMatrixReduction( Metric M1,  Metric M2, D2xD2 &V) {/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/ReductionSimultanee)*/

		/*In this routine we must return a matrix V that is composed of the 
		 * eigen vectors of N=inv(M1) M2.
		 * Instead of looking at N directly, we are going to use the fact that
		 * M1 and M2 are symmetrical, positive definite. 
		 * The eigen values of N are given by solving
		 *    inv(M1) M2 V = lambda V
		 * which is equivalent to
		 *    M2 V = lambda M1 V
		 * and we will hence solve
		 *    (M2 - lambda M1) V = 0
		 */

		//M1 and M2 components
		double a11=M1.a11,a21=M1.a21,a22=M1.a22;
		double b11=M2.a11,b21=M2.a21,b22=M2.a22;

		/*To get the eigen values, we solve the following problem:
		 *    det(M2-lambda M1) = 0
		 *    (b11 - lambda a11)(b22-lambda a22) - (b21-lambda a21)^2
		 * and we have the following trinome:
		 *    a lambda^2 + b lambda + c =0
		 * with:
		 *    a = a11 a22 - a21 a21 (=det(M1))
		 *    b = -a11 b22 -b11 a22 + 2 b21 a21
		 *    c = b11 b22 - b21 b21 (=det(M2))
		 *    */
		const double a= a11*a22  - a21*a21;
		const double b=-a11*b22 - b11*a22+2*b21*a21;
		const double c=-b21*b21 + b11*b22;
		const double bb=b*b,ac=a*c;
		const double delta= bb-4*ac;

		// first, case of a double root if:
		//  - all the terms are very small (a??)
		//  - or : delta is very small
		if ( (bb + Abs(ac) < 1.0e-34 ) ||  (delta < 1.0e-6*bb) ){
			//all vectors are eigen vectors -> choose 1,0 and 0,1
			V= D2xD2(1,0,0,1);
		}

		//general case: two distinct roots: lambda1 and lambda2
		else {

			/*Compute eigen values*/
			const double delta2 = sqrt(delta);
			double lambda[2];
			lambda[0]= (-b - delta2)/(2*a);
			lambda[1]= (-b + delta2)/(2*a);

			/*compute eigen vectors*/
			double vp[2][2];
			double v0,v1,v2;
			double s0,s1;

			for(int i=0;i<2;i++){
				/*Now, one must find the eigen vectors. For that we use the 
				 * following property of the inner product
				 *    (Ax,b) = transp(b) Ax = transp(x) transp(A) b
				 *           = (transp(A) b ,x)
				 * Here we are dealing with A= M2 - lambda M1 which is symmetrical:
				 *    for all (x,y) in R2 
				 *       ((M2 - lambda M1)x,y)=((M2 - lambda M1)y,x)
				 * If y is in Ker(M2 - lambda M1):
				 *    for all x in R2
				 *       ((M2 - lambda M1)y,x)=0
				 * This shows that:
				 *    Ker(M2 - lambda M1) is orthogonal to Im(M2 - lambda M1)
				 * To find the eigen vectors, we only have to find two vectors
				 * of the image and take their perpendicular as long as they are
				 * not 0.
				 * To do that, we take (1,0) and (0,1) and take the larger norm*/

				//compute V = M2 - lambdai M1
				v0 = b11 - lambda[i]*a11;
				v1 = b21 - lambda[i]*a21;
				v2 = b22 - lambda[i]*a22;

				// compute s1=norm(V(1,0)) and s0=norm(V(0,1))
				s0 = v0*v0 + v1*v1;
				s1 = v1*v1 + v2*v2;

				//compute vp1 = (vp1x,vp1y)
				if(s1 < s0){
					s0=sqrt(s0);
					vp[0][i]=   v1/s0;
					vp[1][i]= - v0/s0;
				}
				else{
					s1=sqrt(s1);
					vp[0][i]=   v2/s1;
					vp[1][i]= - v1/s1;
				}
			}

			//compute V from vp
			V=D2xD2(vp[0][0],vp[0][1],vp[1][0],vp[1][1]);
		}
	}
	/*}}}*/
	double abscisseInterpole(const Metric& Ma,const  Metric& Mb, R2 AB,double s,int optim) { /*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Metric.cpp/abscisseInterpole)*/

		if(!optim)  LengthInterpole(Ma,Mb,AB);
		double l  = s* LastMetricInterpole.lab,r;
		int j=LastMetricInterpole.opt-1;

		double * L= LastMetricInterpole.L, *S = LastMetricInterpole.S;
		// warning for optimisation S is the abcisse in [0:0.5]
		// and L is le lenght 
		if(l<=L[0]){
			r=2*S[0]*l/L[0];
		}
		else if (l>=L[j]){
			r=1;
		}
		else{
			int i=0;
			while (j-i>1){
				int k;
				k= (i+j)/2;
				if(l<=L[k]){
					j=k;// l<=L[j] 
				}
				else{
					i=k; //  L[i]<l
				}
			};
			if (i==j){
				r = 2*S[i];
			}
			else{
				r =  2*(S[i]*(L[j]-l)+ S[j]*(l-L[i]))/(L[j]-L[i]);
			}
		}
		if (r>1 || r<0){
			_error_("r>1 || r<0");
		}
		return r ;
	}
	/*}}}*/

}
