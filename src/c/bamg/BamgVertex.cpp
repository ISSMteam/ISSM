#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

#include "./bamgobjects.h"
#include "../shared/shared.h"
#include "./det.h"

namespace bamg {

	/*Constructor/Destructor*/
	BamgVertex::BamgVertex(){ /*{{{*/
		this->PreviousNumber = 0;
	}/*}}}*/

	/*Methods*/
	void BamgVertex::Echo(void){/*{{{*/

		_printf_("Vertex:\n");
		_printf_("  integer   coordinates i.x: " << i.x << ", i.y: " << i.y << "\n");
		_printf_("  Euclidean coordinates r.x: " << r.x << ", r.y: " << r.y << "\n");
		_printf_("  ReferenceNumber = " << ReferenceNumber << "\n");
		_printf_("  PreviousNumber  = " << PreviousNumber << "\n");
		m.Echo();

		return;
	}
	/*}}}*/
	int  BamgVertex::GetReferenceNumber() const { /*{{{*/
		return ReferenceNumber;
	}
	/*}}}*/
	void BamgVertex::MetricFromHessian(const double Hxx,const double Hyx, const double Hyy,const double smin,const double smax,const double s,double err,BamgOpts* bamgopts){/*{{{*/
		/*Compute Metric from Hessian*/

		/*get options*/
		double power=(bamgopts->power);
		double anisomax=(bamgopts->anisomax);
		double CutOff=bamgopts->cutoff;
		double hmin=(bamgopts->hmin);
		double hmax=(bamgopts->hmax);
		double coef=bamgopts->coeff;
		int    Metrictype=(bamgopts->Metrictype);

		/*Intermediary*/
		double ci;

		/*compute multiplicative coefficient depending on Metric Type (2/9 because it is 2d)*/

		//Absolute Error
		/*
		 *            2         1       
		 *Metric M = ---  ------------   Abs(Hessian)
		 *            9   err * coeff^2  
		 */
		if (Metrictype==0){
			ci= 2.0/9.0 * 1/(err*coef*coef);
		}

		//Relative Error
		/*
		 *            2         1            Abs(Hessian)
		 *Metric M = ---  ------------  ----------------------
		 *            9   err * coeff^2  max( |s| , cutoff*max(|s|) )
		 *
		 */
		else if (Metrictype==1){
			ci= 2.0/9.0 * 1/(err*coef*coef) * 1/Max( Abs(s), CutOff*(Max(Abs(smin),Abs(smax))));
		}

		//Rescaled absolute error
		/*
		 *            2         1            Abs(Hessian)
		 *Metric M = ---  ------------  ---------------------- 
		 *            9   err * coeff^2       (smax-smin)
		 */
		else if (Metrictype==2){
			ci= 2.0/9.0 * 1/(err*coef*coef) * 1/(smax-smin);
		}
		else{
			_error_("Metrictype " << Metrictype << " not supported yet (use 0,1 or 2(default))");
		}

		//initialize metric Miv with ci*H
		Metric Miv(Hxx*ci,Hyx*ci,Hyy*ci);

		//Get eigen values and vectors of Miv
		EigenMetric Vp(Miv);

		//move eigen valuse to their absolute values
		Vp.Abs();

		//Apply a power if requested by user
		if(power!=1.0) Vp.pow(power);

		//modify eigen values according to hmin and hmax
		Vp.Maxh(hmax);
		Vp.Minh(hmin);

		//Bound anisotropy by 1/(anisomax)^2
		Vp.BoundAniso2(1/(anisomax*anisomax));

		//rebuild Metric from Vp
		Metric MVp(Vp);

		//Apply Metric to vertex
		m.IntersectWith(MVp);

	}
	/*}}}*/
	long BamgVertex::Optim(int i,int koption){ /*{{{*/
		long ret=0;
		if ( t && (IndexInTriangle >= 0 ) && (IndexInTriangle <3) ){
			ret = t->Optim(IndexInTriangle,koption);
			if(i==0){
				t =0; // for no future optim
				IndexInTriangle= 0;
			}
		}
		return ret;
	}
	/*}}}*/
	double  BamgVertex::Smoothing(Mesh &Th,Mesh &BTh,Triangle* &tstart ,double omega){/*{{{*/
		/*Original code from Frederic Hecht <hecht@ann.jussieu.fr> (BAMG v1.01, Mesh2.cpp/Smoothing)*/

		BamgVertex* s=this;
		BamgVertex &vP = *s,vPsave=vP;

		Triangle* tbegin= t , *tria = t , *ttc;

		int k=0,kk=0,j = EdgesVertexTriangle[IndexInTriangle][0],jc;
		R2 P(s->r),PNew(0,0);
		do {
			k++; 

			if (!tria->Hidden(j)){
				BamgVertex &vQ = (*tria)[VerticesOfTriangularEdge[j][0]]; 

				R2 Q = vQ,QP(P-Q);
				double lQP = LengthInterpole(vP,vQ,QP);
				PNew += Q+QP/Max(lQP,1e-20);
				kk ++;
			}
			ttc =  tria->TriangleAdj(j);
			jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
			tria = ttc;
			j = NextEdge[jc];
			if (k>=2000){
				_error_("k>=2000 (Maximum number of iterations reached)");
			}
		} while ( tbegin != tria); 
		if (kk<4) return 0;
		PNew = PNew/(double)kk;
		R2 Xmove((PNew-P)*omega);
		PNew = P+Xmove;
		double delta=Norme2_2(Xmove); 

		long long deta[3];
		I2 IBTh  = BTh.R2ToI2(PNew);

		tstart=BTh.TriangleFindFromCoord(IBTh,deta,tstart);  

		if (tstart->det <0){ // outside
			double ba,bb;
			AdjacentTriangle edge= CloseBoundaryEdge(IBTh,tstart,ba,bb) ;
			tstart = edge;
			vP.m= Metric(ba,*edge.EdgeVertex(0),bb,*edge.EdgeVertex(1));
		}
		else { // inside
			double   aa[3];
			double s = deta[0]+deta[1]+deta[2];
			aa[0]=deta[0]/s;
			aa[1]=deta[1]/s;
			aa[2]=deta[2]/s;
			vP.m = Metric(aa,(*tstart)[0],(*tstart)[1],(*tstart)[2]);
		}

		// recompute the det of the triangle
		vP.r = PNew;

		vP.i = Th.R2ToI2(PNew);

		BamgVertex vPnew = vP;

		int ok=1;
		int loop=1;
		k=0;
		while (ok){
			ok =0;
			do {
				k++; 
				double detold = tria->det;
				tria->det =  bamg::det( (*tria)[0],(*tria)[1]  ,(*tria)[2]);
				if (loop) {
					if (tria->det<0) ok =1;			       
					else if ( (double)tria->det < detold/2 ) ok=1;
				}
				tria->SetUnMarkUnSwap(0);
				tria->SetUnMarkUnSwap(1);
				tria->SetUnMarkUnSwap(2);
				ttc =  tria->TriangleAdj(j);
				jc = NextEdge[tria->NuEdgeTriangleAdj(j)];
				tria = ttc;
				j = NextEdge[jc];
				if (k>=2000){
					_error_("k>=2000");
				}
			}while ( tbegin != tria); 

			if (ok && loop) vP=vPsave; // no move 
			loop=0;
		}
		return delta;
	}
	/*}}}*/
} 
