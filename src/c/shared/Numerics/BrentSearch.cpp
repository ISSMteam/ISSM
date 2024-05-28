/*!\file:  BrentSearch.cpp
 * \brief optimization algorithm
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <float.h>
#include <iomanip>
#include <cmath>

#include "../Exceptions/exceptions.h"
#include "../io/io.h"
#include "../MemOps/MemOps.h"
#include "./Verbosity.h"
#include "./OptPars.h"
#include "./types.h"
#include "./isnan.h"

void BrentSearch(IssmDouble** pJ,OptPars optpars,IssmDouble* X0,IssmDouble (*f)(IssmDouble*,void*),IssmDouble (*g)(IssmDouble**,IssmDouble*,void*),void* usr){

	/* This routine is optimizing a given function using Brent's method
	 * (Golden or parabolic procedure)*/

	/*Intermediary*/
	int        iter;
	IssmDouble si,gold,intervalgold,oldintervalgold;
	IssmDouble parab_num,parab_den,distance;
	IssmDouble fxmax,fxmin,fxbest;
	IssmDouble fx,fx1,fx2;
	IssmDouble x,x1,x2,xm,xbest;
	IssmDouble tol1,tol2,seps;
	IssmDouble tolerance = 1.e-4;

	/*Recover parameters:*/
	int         nsteps  = optpars.nsteps;
	int         nsize   = optpars.nsize;
	IssmDouble  xmin    = optpars.xmin;
	IssmDouble  xmax    = optpars.xmax;
	int        *maxiter = optpars.maxiter;
	IssmDouble *cm_jump = optpars.cm_jump;

	/*Initialize gradient and controls*/
	IssmDouble* G = NULL;
	IssmDouble* J = xNew<IssmDouble>(nsteps);
	IssmDouble* X = xNew<IssmDouble>(nsize);

	/*Header of printf*/
	_printf0_("\n");
	_printf0_("       x       |  Cost function f(x)  |  List of contributions\n");

	/*start iterations*/
	for(int n=0;n<nsteps;n++){

		/*Print iteration number*/
		_printf0_("====================== step "<< n+1 << "/" << nsteps <<" ===============================\n");

		/*Reset some variables*/
		iter = 0;
		xmin = 0.;
		xmax = 1.;
		bool loop = true;
		cout<<setprecision(5);

		/*Get current Gradient at xmin=0*/
		_printf0_(" x = "<<setw(9)<<xmin<<" | ");
		fxmin = (*g)(&G,X0,usr); if(xIsNan<IssmDouble>(fxmin)) _error_("Function evaluation returned NaN");

		/*Get f(xmax)*/
		_printf0_(" x = "<<setw(9)<<xmax<<" | ");
		for(int i=0;i<nsize;i++) X[i]=X0[i]+xmax*G[i];
		fxmax = (*f)(X,usr); if (xIsNan<IssmDouble>(fxmax)) _error_("Function evaluation returned NaN");
		//if(VerboseControl()) _printf0_("           N/A    "<<setw(12)<<xmax<<"  "<<setw(12)<<fxmax<<"           N/A         boundary\n");

		/*test if jump option activated and xmin==0*/
		if(!xIsNan<IssmDouble>(cm_jump[n]) && (xmin==0) && (fxmax/fxmin)<cm_jump[n]){
			for(int i=0;i<nsize;i++) X0[i]=X0[i]+xmax*G[i];
			xDelete<IssmDouble>(G);
			J[n]=fxmax;
			continue;
		}

		/*initialize optimization variables*/
		seps=sqrt(DBL_EPSILON);    //precision of a IssmDouble
		distance=0.0;              //new_x=old_x + distance
		gold=0.5*(3.0-sqrt(5.0));  //gold = 1 - golden ratio
		intervalgold=0.0;          //distance used by Golden procedure

		/*1: initialize the values of the 4 x needed (x1,x2,x,xbest)*/
		x1=xmin+gold*(xmax-xmin);
		x2=x1;
		xbest=x1;
		x=xbest;

		/*2: call the function to be evaluated*/
		_printf0_(" x = "<<setw(9)<<x<<" | ");
		for(int i=0;i<nsize;i++) X[i]=X0[i]+x*G[i];
		fxbest = (*f)(X,usr); if(xIsNan<IssmDouble>(fxbest)) _error_("Function evaluation returned NaN");
		iter++;

		/*3: update the other variables*/
		fx1=fxbest;
		fx2=fxbest;
		/*xm is always in the middle of a and b*/
		xm=0.5*(xmin+xmax);                           
		/*update tolerances*/
		tol1=seps*sqrt(pow(xbest,2))+tolerance/3.0;
		tol2=2.0*tol1;

		/*4: print result*/
		//if(VerboseControl())
		 //_printf0_("         "<<setw(5)<<iter<<"    "<<setw(12)<<xbest<<"  "<<setw(12)<<fxbest<<"  "<<setw(12)<<pow(pow(xbest-xm,2),0.5)<<"         initial\n");
		if (!xIsNan<IssmDouble>(cm_jump[n]) && (xmin==0) && ((fxbest/fxmin)<cm_jump[n])){
			//if(VerboseControl()) _printf0_("      optimization terminated: current x satisfies criteria 'cm_jump'=" << cm_jump[n] << "\n");
			loop=false;
		}

		while(loop){

			bool goldenflag=true;

			// Is a parabolic fit possible ?
			if (sqrt(pow(intervalgold,2))>tol1){

				// Yes, so fit parabola
				goldenflag=false;
				parab_num=(xbest-x1)*(xbest-x1)*(fxbest-fx2)-(xbest-x2)*(xbest-x2)*(fxbest-fx1);;
				parab_den=2.0*(xbest-x1)*(fxbest-fx2)-2.0*(xbest-x2)*(fxbest-fx1);

				//reverse p if necessary
				if(parab_den>0.0){ 
					parab_num=-parab_num;
				}
				parab_den=sqrt(pow(parab_den,2));
				oldintervalgold=intervalgold;
				intervalgold=distance;

				// Is the parabola acceptable (we use seps here because in some configuration parab_num==parab_den*(xmax-xbest)
				// and the result is not repeatable anymore
				if (( sqrt(pow(parab_num,2)) < sqrt(pow(0.5*parab_den*oldintervalgold,2))) &&
							(parab_num>parab_den*(xmin-xbest)+seps) && 
							(parab_num<parab_den*(xmax-xbest)-seps)){

					// Yes, parabolic interpolation step
					distance=parab_num/parab_den;
					x=xbest+distance;

					// f must not be evaluated too close to min_x or max_x
					if (((x-xmin)<tol2) || ((xmax-x)<tol2)){
						if ((xm-xbest)<0.0) si=-1;
						else                si=1;
						//compute new distance
						distance=tol1*si;
					}
				}
				else{
					// Not acceptable, must do a golden section step
					goldenflag=true;
				}
			}

			//Golden procedure
			if(goldenflag){
				// compute the new distance d
				if(xbest>=xm){
					intervalgold=xmin-xbest;    
				}
				else{ 
					intervalgold=xmax-xbest;  
				}
				distance=gold*intervalgold;
			}

			// The function must not be evaluated too close to xbest
			if(distance<0) si=-1;
			else           si=1;
			if(sqrt(pow(distance,2))>tol1) x=xbest+si*sqrt(pow(distance,2));
			else                           x=xbest+si*tol1;

			//evaluate function on x
			_printf0_(" x = "<<setw(9)<<x<<" | ");
			for(int i=0;i<nsize;i++) X[i]=X0[i]+x*G[i];
			fx = (*f)(X,usr); if(xIsNan<IssmDouble>(fx)) _error_("Function evaluation returned NaN");
			iter++;

			// Update a, b, xm, x1, x2, tol1, tol2
			if (fx<=fxbest){
				if (x>=xbest) xmin=xbest;
				else          xmax=xbest;
				x1=x2;    fx1=fx2;
				x2=xbest; fx2=fxbest;
				xbest=x;  fxbest=fx;
			}
			else{ // fx > fxbest
				if (x<xbest) xmin=x;
				else         xmax=x;
				if ((fx<=fx2) || (x2==xbest)){
					x1=x2; fx1=fx2;
					x2=x;  fx2=fx;
				}
				else if ( (fx <= fx1) || (x1 == xbest) || (x1 == x2) ){
					x1=x;  fx1=fx;
				}
			}
			xm = 0.5*(xmin+xmax);
			tol1=seps*pow(pow(xbest,2),0.5)+tolerance/3.0;
			tol2=2.0*tol1;
			//if(VerboseControl())
			// _printf0_("         "<<setw(5)<<iter<<"    "<<setw(12)<<x<<"  "<<setw(12)<<fx<<"  "<<setw(12)<<pow(pow(xbest-xm,2),0.5)<<
			//			 "         "<<(goldenflag?"golden":"parabolic")<<"\n");

			/*Stop the optimization?*/
			if (sqrt(pow(xbest-xm,2)) < (tol2-0.5*(xmax-xmin))){
				//if(VerboseControl()) _printf0_("      optimization terminated: current x satisfies criteria 'tolx'=" << tolerance << "\n");
				loop=false;
			}
			else if (iter>=maxiter[n]){
				//if(VerboseControl()) _printf0_("      exiting: Maximum number of iterations has been exceeded  ('maxiter'=" << maxiter[n] << ")\n");
				loop=false;
			}
			else if (!xIsNan<IssmDouble>(cm_jump[n]) && (xmin==0) && ((fxbest/fxmin)<cm_jump[n])){
				//if(VerboseControl()) _printf0_("      optimization terminated: current x satisfies criteria 'cm_jump'=" << cm_jump[n] << "\n");
				loop=false;
			}
			else{
				//continue
				loop=true;
			}
		}//end while

		//Now, check that the value on the boundaries are not better than current fxbest
		if (fxbest>fxmin){
			xbest=optpars.xmin; fxbest=fxmin;
		}
		if (fxbest>fxmax){
			xbest=optpars.xmax; fxbest=fxmax;
		}

		/*Assign output pointers: */
		for(int i=0;i<nsize;i++) X0[i]=X0[i]+xbest*G[i];
		xDelete<IssmDouble>(G);
		J[n]=fxbest;
	}

	/*return*/
	xDelete<IssmDouble>(X);
	*pJ=J;
}
