#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "Metric.h"
#include "../shared/shared.h"

namespace bamg {

	/*Constructor*/
	EigenMetric::EigenMetric(const Metric& M){/*{{{*/
		/*From a metric (a11,a21,a22), get eigen values lambda1 and lambda2 and one eigen vector v*/

		/*Intermediaries*/
		double a11=M.a11,a21=M.a21,a22=M.a22;
		double normM;
		double delta,b;

		/*To get the eigen values, we must solve the following equation:
		 *     | a11 - lambda    a21        |
		 * det |                            | = 0
		 *     | a21             a22-lambda |
		 *
		 * We have to solve the following polynom:
		 *  lamda^2 + ( -a11 -a22)*lambda + (a11*a22-a21*a21) = 0*/

		/*Compute polynom determinant*/
		b=-a11-a22;
		delta=b*b - 4*(a11*a22-a21*a21);

		/*Compute norm of M to avoid round off errors*/
		normM=a11*a11 + a22*a22 + a21*a21;

		/*1: normM too small: eigen values = 0*/
		if(normM<1.e-30){
			lambda1=0;
			lambda2=0;
			vx=1;
			vy=0;
		}
		/*2: delta is small -> double root*/
		else if (delta < 1.e-5*normM){
			lambda1=-b/2;
			lambda2=-b/2;
			vx=1;
			vy=0;
		}
		/*3: general case -> two roots*/
		else{
			delta     = sqrt(delta);
			lambda1   = (-b-delta)/2.0;
			lambda2   = (-b+delta)/2.0;

			/*Now, one must find the eigen vectors. For that we use the following property of the inner product
			 *    <Ax,y> = <x,tAy>
			 * Here, M'(M-lambda*Id) is symmetrical, which gives:
			 *    ∀(x,y)∈R²xR² <M'x,y> = <M'y,x>
			 * And we have the following:
			 *    if y∈Ker(M'), ∀x∈R² <M'x,y> = <x,M'y> = 0
			 * We have shown that
			 *    Im(M') ⊥ Ker(M')
			 *
			 * To find the eigen vectors of M, we only have to find two vectors
			 * of the image of M' and take their perpendicular as long as they are
			 * not 0.
			 * To do that, we take the images (1,0) and (0,1):
			 *  x1 = (a11 - lambda)      x2 = a21
			 *  y1 = a21                 y2 = (a22-lambda)
			 *
			 * We take the vector that has the larger norm and take its perpendicular.*/

			double norm1 = (a11-lambda1)*(a11-lambda1) + a21*a21; 
			double norm2 = a21*a21 + (a22-lambda1)*(a22-lambda1);

			if (norm2<norm1){
				norm1=sqrt(norm1);
				vx = - a21/norm1;
				vy = (a11-lambda1)/norm1;
			}
			else{
				norm2=sqrt(norm2);
				vx = - (a22-lambda1)/norm2;
				vy = a21/norm2;
			}
		}

	}
	/*}}}*/
	EigenMetric::EigenMetric(double r1,double r2,const D2& vp1){/*{{{*/
		this->lambda1 = r1;
		this->lambda2 = r2;
		this->vx = vp1.x;
		this->vy = vp1.y;
	}/*}}}*/

	/*Methods*/
	void   EigenMetric::Abs(){/*{{{*/
		lambda1=bamg::Abs(lambda1),lambda2=bamg::Abs(lambda2);
	}/*}}}*/
	double EigenMetric::Aniso2() const  { /*{{{*/
		return lmax()/lmin();
	}/*}}}*/
	void EigenMetric::Echo(void){/*{{{*/

		_printf_("EigenMetric:\n");
		_printf_("   lambda1: " << lambda1 << "\n");
		_printf_("   lambda2: " << lambda2 << "\n");
		_printf_("   vx: " << vx << "\n");
		_printf_("   vy: " << vy << "\n");

		return;
	}
	/*}}}*/
	double EigenMetric::hmin() const {/*{{{*/
		return sqrt(1/bamg::Max3(lambda1,lambda2,1e-30));
	}/*}}}*/
	double EigenMetric::hmax() const {/*{{{*/
		return sqrt(1/bamg::Max(bamg::Min(lambda1,lambda2),1e-30));
	}/*}}}*/
	double EigenMetric::lmax() const {/*{{{*/
		return bamg::Max3(lambda1,lambda2,1e-30);
	}/*}}}*/
	double EigenMetric::lmin() const {/*{{{*/
		return bamg::Max(bamg::Min(lambda1,lambda2),1e-30);
	}/*}}}*/
	void   EigenMetric::Min(double a) { /*{{{*/
		lambda1=bamg::Min(a,lambda1); lambda2=bamg::Min(a,lambda2) ;
	}/*}}}*/
	void   EigenMetric::Max(double a) { /*{{{*/
		//change eigen values
		lambda1=bamg::Max(a,lambda1); lambda2=bamg::Max(a,lambda2) ;
	}/*}}}*/
	void   EigenMetric::Minh(double h) {/*{{{*/
		Min(1.0/(h*h));
	}/*}}}*/
	void   EigenMetric::Maxh(double h) {/*{{{*/
		//Call Max function
		Max(1.0/(h*h));
	}/*}}}*/
	void   EigenMetric::pow(double p){/*{{{*/
		lambda1=::pow(lambda1,p);lambda2=::pow(lambda2,p);
	}/*}}}*/

} 
