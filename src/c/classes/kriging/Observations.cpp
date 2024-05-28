/*
 * \file Observations.cpp
 * \brief: Implementation of Observations class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

#include "../Options/Options.h"
#include "./Observations.h"
#include "./Observation.h"
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"

#include "./Quadtree.h"
#include "./Covertree.h"
#include "./Variogram.h"
#include "../../toolkits/toolkits.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Observations::Observations(){/*{{{*/
	this->treetype  = 0;
	this->quadtree  = NULL;
	this->covertree = NULL;
	return;
}
/*}}}*/
Observations::Observations(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options){/*{{{*/

	/*Check that there are observations*/
	if(n<=0) _error_("No observation found");

	/*Get tree type (FIXME)*/
	double dtree = 0.;
	options->Get(&dtree,"treetype",1.);
	this->treetype = reCast<int>(dtree);
	switch(this->treetype){
		case 1:
			this->covertree = NULL;
			this->InitQuadtree(observations_list,x,y,n,options);
			break;
		case 2:
			this->quadtree = NULL;
			this->InitCovertree(observations_list,x,y,n,options);
			break;
		default:
			_error_("Tree type "<<this->treetype<<" not supported yet (1: quadtree, 2: covertree)");
	}
}
/*}}}*/
Observations::~Observations(){/*{{{*/
	switch(this->treetype){
		case 1:
			delete this->quadtree;
			break;
		case 2:
			delete this->covertree;
			break;
		default:
			_printf_("Tree type "<<this->treetype<<" not supported yet (1: quadtree, 2: covertree)");
	}
	return;
}
/*}}}*/

/*Initialize data structures*/
void Observations::InitCovertree(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options){/*{{{*/

	/*Intermediaries*/
	 IssmPDouble  minspacing,mintrimming,maxtrimming;

	/*Checks*/
	_assert_(n);

	/*Get trimming limits*/
	options->Get(&mintrimming,"mintrimming",-1.e+21);
	options->Get(&maxtrimming,"maxtrimming",+1.e+21);
	options->Get(&minspacing,"minspacing",0.01);
	if(minspacing<=0) _error_("minspacing must > 0");

	/*Get maximum distance between 2 points
	 *  maxDist should be the maximum distance that any two points
	 *  can have between each other. IE p.distance(q) < maxDist for all
	 *  p,q that you will ever try to insert. The cover tree may be invalid
	 *  if an inaccurate maxDist is given.*/
	IssmPDouble xmin = x[0];
	IssmPDouble xmax = x[0];
	IssmPDouble ymin = y[0];
	IssmPDouble ymax = y[0];
	for(int i=1;i<n;i++){
		if(x[i]<xmin) xmin=x[i];
		if(x[i]>xmax) xmax=x[i];
		if(y[i]<ymin) ymin=y[i];
		if(y[i]>ymax) ymax=y[i];
	}
	IssmPDouble maxDist = sqrt(pow(xmax-xmin,2)+pow(ymax-ymin,2));
	IssmPDouble base    = 2.;
	int         maxdepth = ceilf(log(maxDist)/log(base));

	 _printf0_("Generating covertree with a maximum depth " <<  maxdepth <<"... ");
    this->covertree=new Covertree(maxdepth);

    for(int i=0;i<n;i++){

		/*First check limits*/
		if(observations_list[i]>maxtrimming) continue;
		if(observations_list[i]<mintrimming) continue;

		/*Second, check that this observation is not too close from another one*/
		Observation newobs = Observation(x[i],y[i],observations_list[i]);
		if(i>0 && this->covertree->getRoot()){
			/*Get closest obs and see if it is too close*/
			std::vector<Observation> kNN=(this->covertree->kNearestNeighbors(newobs,1));
			Observation oldobs = (*kNN.begin());
			if(oldobs.distance(newobs)<minspacing) continue;
		}

		this->covertree->insert(newobs);
    }
	 _printf0_("done\n");
}
/*}}}*/
void Observations::InitQuadtree(IssmPDouble* observations_list,IssmPDouble* x,IssmPDouble* y,int n,Options* options){/*{{{*/

	/*Intermediaries*/
	int          i,maxdepth,level,counter,index;
	int          xi,yi;
	IssmPDouble  xmin,xmax,ymin,ymax;
	IssmPDouble  offset,minlength,minspacing,mintrimming,maxtrimming;
	Observation *observation = NULL;

	/*Checks*/
	_assert_(n);

	/*Get extrema*/
	xmin=x[0]; ymin=y[0];
	xmax=x[0]; ymax=y[0];
	for(i=1;i<n;i++){
		xmin=min(xmin,x[i]); ymin=min(ymin,y[i]);
		xmax=max(xmax,x[i]); ymax=max(ymax,y[i]);
	}
	offset=0.05*(xmax-xmin); xmin-=offset; xmax+=offset;
	offset=0.05*(ymax-ymin); ymin-=offset; ymax+=offset;

	/*Get trimming limits*/
	options->Get(&mintrimming,"mintrimming",-1.e+21);
	options->Get(&maxtrimming,"maxtrimming",+1.e+21);
	options->Get(&minspacing,"minspacing",0.01);
	if(minspacing<=0) _error_("minspacing must > 0");

	/*Get Minimum box size*/
	if(options->GetOption("boxlength")){
		options->Get(&minlength,"boxlength");
		if(minlength<=0)_error_("boxlength should be a positive number");
		maxdepth=reCast<int,IssmPDouble>(log(max(xmax-xmin,ymax-ymin)/minlength +1)/log(2.0));
	}
	else{
		maxdepth = 30;
		minlength=max(xmax-xmin,ymax-ymin)/IssmPDouble((1L<<maxdepth)-1);
	}

	/*Initialize Quadtree*/
	_printf0_("Generating quadtree with a maximum box size " << minlength << " (depth=" << maxdepth << ")... ");
	this->quadtree = new Quadtree(xmin,xmax,ymin,ymax,maxdepth);

	/*Add observations one by one*/
	counter = 0;
	for(i=0;i<n;i++){

		/*First check limits*/
		if(observations_list[i]>maxtrimming) continue;
		if(observations_list[i]<mintrimming) continue;

		/*Second, check that this observation is not too close from another one*/
		this->quadtree->ClosestObs(&index,x[i],y[i]);
		if(index>=0){
			observation=xDynamicCast<Observation*>(this->GetObjectByOffset(index));
			if(pow(observation->x-x[i],2)+pow(observation->y-y[i],2) < minspacing) continue;
		}

		this->quadtree->IntergerCoordinates(&xi,&yi,x[i],y[i]);
		this->quadtree->QuadtreeDepth2(&level,xi,yi);
		if((int)level <= maxdepth){
			observation = new Observation(x[i],y[i],xi,yi,counter++,observations_list[i]);
			this->quadtree->Add(observation);
			this->AddObject(observation);
		}
		else{
			/*We need to average with the current observations*/
			this->quadtree->AddAndAverage(x[i],y[i],observations_list[i]);
		}
	}
	_printf0_("done\n");
	_printf0_("Initial number of observations: " << n << "\n");
	_printf0_("  Final number of observations: " << this->quadtree->NbObs << "\n");
}
/*}}}*/

/*Methods*/
void Observations::ClosestObservation(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius){/*{{{*/

	switch(this->treetype){
		case 1:
			this->ClosestObservationQuadtree(px,py,pobs,x_interp,y_interp,radius);
			break;
		case 2:
			this->ClosestObservationCovertree(px,py,pobs,x_interp,y_interp,radius);
			break;
		default:
			_error_("Tree type "<<this->treetype<<" not supported yet (1: quadtree, 2: covertree)");
	}

}/*}}}*/
void Observations::ClosestObservationCovertree(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius){/*{{{*/

	IssmPDouble hmin  = UNDEF;

	if(this->covertree->getRoot()){
		/*Get closest obs and see if it is too close*/
		Observation newobs = Observation(x_interp,y_interp,0.);
		std::vector<Observation> kNN=(this->covertree->kNearestNeighbors(newobs,1));
		Observation observation = (*kNN.begin());
		hmin = observation.distance(newobs);
		if(hmin<=radius){
			*px   = observation.x;
			*py   = observation.y;
			*pobs = observation.value;
			return;
		}
	}

	*px   = UNDEF;
	*py   = UNDEF;
	*pobs = UNDEF;
}/*}}}*/
void Observations::ClosestObservationQuadtree(IssmPDouble *px,IssmPDouble *py,IssmPDouble *pobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius){/*{{{*/

	/*Output and Intermediaries*/
	int          nobs,i,index;
	IssmPDouble  hmin,h2,hmin2;
	int         *indices      = NULL;
	Observation *observation  = NULL;

	/*If radius is not provided or is 0, return all observations*/
	if(radius==0) radius=this->quadtree->root->length;

	/*For CPPcheck*/
	hmin = 2*radius;

	/*First, find closest point in Quadtree (fast but might not be the true closest obs)*/
	this->quadtree->ClosestObs(&index,x_interp,y_interp);
	if(index>=0){
		observation=xDynamicCast<Observation*>(this->GetObjectByOffset(index));
		hmin = sqrt((observation->x-x_interp)*(observation->x-x_interp) + (observation->y-y_interp)*(observation->y-y_interp));
		if(hmin<radius) radius=hmin;
	}

	/*Find all observations that are in radius*/
	this->quadtree->RangeSearch(&indices,&nobs,x_interp,y_interp,radius);
	for (i=0;i<nobs;i++){
		observation=xDynamicCast<Observation*>(this->GetObjectByOffset(indices[i]));
		h2 = (observation->x-x_interp)*(observation->x-x_interp) + (observation->y-y_interp)*(observation->y-y_interp);
		if(i==0){
			hmin2 = h2;
			index = indices[i];
		}
		else{
			if(h2<hmin2){
				hmin2 = h2;
				index = indices[i];
			}
		}
	}

	/*Assign output pointer*/
	if(nobs || hmin==radius){
		observation=xDynamicCast<Observation*>(this->GetObjectByOffset(index));
		*px   = observation->x;
		*py   = observation->y;
		*pobs = observation->value;
	}
	else{
		*px   = UNDEF;
		*py   = UNDEF;
		*pobs = UNDEF;
	}
	xDelete<int>(indices);

}/*}}}*/
void Observations::Distances(IssmPDouble* distances,IssmPDouble *x,IssmPDouble *y,int n,IssmPDouble radius){/*{{{*/

	IssmPDouble xi,yi,obs;

	for(int i=0;i<n;i++){
		this->ClosestObservation(&xi,&yi,&obs,x[i],y[i],radius);
		if(xi==UNDEF && yi==UNDEF){
		 distances[i]=UNDEF;
		}
		else{
		 distances[i]=sqrt( (x[i]-xi)*(x[i]-xi) + (y[i]-yi)*(y[i]-yi) );
		}
	}
}/*}}}*/
void Observations::ObservationList(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs){/*{{{*/

	/*Output and Intermediaries*/
	int          nobs;
	IssmPDouble *x            = NULL;
	IssmPDouble *y            = NULL;
	IssmPDouble *obs          = NULL;
	Observation *observation  = NULL;

	nobs = this->Size();

	if(nobs){
		x   = xNew<IssmPDouble>(nobs);
		y   = xNew<IssmPDouble>(nobs);
		obs = xNew<IssmPDouble>(nobs);
		int i=0;
		for(Object* & object: this->objects){
			observation=xDynamicCast<Observation*>(object);
			observation->WriteXYObs(&x[i],&y[i],&obs[i]);
			i++;
		}
	}

	/*Assign output pointer*/
	*px=x;
	*py=y;
	*pobs=obs;
	*pnobs=nobs;
}/*}}}*/
void Observations::ObservationList(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int maxdata){/*{{{*/

	switch(this->treetype){
		case 1:
			this->ObservationListQuadtree(px,py,pobs,pnobs,x_interp,y_interp,radius,maxdata);
			break;
		case 2:
			this->ObservationListCovertree(px,py,pobs,pnobs,x_interp,y_interp,radius,maxdata);
			break;
		default:
			_error_("Tree type "<<this->treetype<<" not supported yet (1: quadtree, 2: covertree)");
	}
}/*}}}*/
void Observations::ObservationListCovertree(double **px,double **py,double **pobs,int* pnobs,double x_interp,double y_interp,double radius,int maxdata){/*{{{*/

	double *x            = NULL;
	double *y            = NULL;
	double *obs          = NULL;
	Observation observation=Observation(x_interp,y_interp,0.);
	std::vector<Observation> kNN;

	kNN=(this->covertree->kNearestNeighbors(observation, maxdata));
	//cout << "kNN's size: " << kNN.size() << " (maxdata = " <<maxdata<<")"<<endl;

	//kNN is sort from closest to farthest neighbor
	//searches for the first neighbor that is out of radius
	//deletes and resizes the kNN vector
	vector<Observation>::iterator it;
	if(radius>0.){
		for (it = kNN.begin(); it != kNN.end(); ++it) {
			//(*it).print();
			//cout << "\n" << (*it).distance(observation) << endl;
			if ((*it).distance(observation) > radius) {
				break;
			}
		}
		kNN.erase(it, kNN.end());
	}

	/*Allocate vectors*/
	x   = new double[kNN.size()];
	y   = new double[kNN.size()];
	obs = new double[kNN.size()];

	/*Loop over all observations and fill in x, y and obs*/
	int i = 0;
	for(it = kNN.begin(); it != kNN.end(); ++it) {
		(*it).WriteXYObs((*it), &x[i], &y[i], &obs[i]);
		i++;
	}

	*px=x;
	*py=y;
	*pobs=obs;
	*pnobs = kNN.size();
}/*}}}*/
void Observations::ObservationListQuadtree(IssmPDouble **px,IssmPDouble **py,IssmPDouble **pobs,int* pnobs,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int maxdata){/*{{{*/

	/*Output and Intermediaries*/
	bool         stop;
	int          nobs,tempnobs,i,j,k,n,counter;
	IssmPDouble  h2,radius2;
	int         *indices      = NULL;
	int         *tempindices  = NULL;
	IssmPDouble *dists        = NULL;
	IssmPDouble *x            = NULL;
	IssmPDouble *y            = NULL;
	IssmPDouble *obs          = NULL;
	Observation *observation  = NULL;

	/*If radius is not provided or is 0, return all observations*/
	if(radius==0.) radius=this->quadtree->root->length*2.;

	/*Compute radius square*/
	radius2 = radius*radius;

	/*Find all observations that are in radius*/
	this->quadtree->RangeSearch(&tempindices,&tempnobs,x_interp,y_interp,radius);
	if(tempnobs){
		indices = xNew<int>(tempnobs);
		dists   = xNew<IssmPDouble>(tempnobs);
	}
	nobs = 0;
	for(i=0;i<tempnobs;i++){
		observation=xDynamicCast<Observation*>(this->GetObjectByOffset(tempindices[i]));
		h2 = (observation->x-x_interp)*(observation->x-x_interp) + (observation->y-y_interp)*(observation->y-y_interp);

		if(nobs==maxdata && h2>radius2) continue;
		if(nobs<maxdata){
			indices[nobs]   = tempindices[i];
			dists[nobs]     = h2;
			nobs++;
		}
		if(nobs==1) continue;

		/*Sort all dists up to now*/
		n=nobs-1;
		stop = false;
		for(k=0;k<n-1;k++){
			if(h2<dists[k]){
				counter=1;
				for(int jj=k;jj<n;jj++){
					j  = n-counter;
					dists[j+1]   = dists[j];
					indices[j+1] = indices[j];
					counter++;
				}
				dists[k]   = h2;
				indices[k] = tempindices[i];
				stop = true;
				break;
			}
			if(stop) break;
		}
	}
	xDelete<IssmPDouble>(dists);
	xDelete<int>(tempindices);

	if(nobs){
		/*Allocate vectors*/
		x   = xNew<IssmPDouble>(nobs);
		y   = xNew<IssmPDouble>(nobs);
		obs = xNew<IssmPDouble>(nobs);

		/*Loop over all observations and fill in x, y and obs*/
		for(i=0;i<nobs;i++){
			observation=xDynamicCast<Observation*>(this->GetObjectByOffset(indices[i]));
			observation->WriteXYObs(&x[i],&y[i],&obs[i]);
		}
	}

	/*Assign output pointer*/
	xDelete<int>(indices);
	*px=x;
	*py=y;
	*pobs=obs;
	*pnobs=nobs;
}/*}}}*/
void Observations::InterpolationIDW(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata,IssmPDouble power){/*{{{*/

	/*Intermediaries*/
	int         i,n_obs;
	IssmPDouble prediction;
	IssmPDouble numerator,denominator,h,weight;
	IssmPDouble *x   = NULL;
	IssmPDouble *y   = NULL;
	IssmPDouble *obs = NULL;

	/*Some checks*/
	_assert_(maxdata>0);
	_assert_(pprediction);
	_assert_(power>0);

	/*Get list of observations for current point*/
	this->ObservationList(&x,&y,&obs,&n_obs,x_interp,y_interp,radius,maxdata);

	/*If we have less observations than mindata, return UNDEF*/
	if(n_obs<mindata){
		prediction = UNDEF;
	}
	else{
		numerator   = 0.;
		denominator = 0.;
		for(i=0;i<n_obs;i++){
			h = sqrt( (x[i]-x_interp)*(x[i]-x_interp) + (y[i]-y_interp)*(y[i]-y_interp));
			if (h<0.0000001){
				numerator   = obs[i];
				denominator = 1.;
				break;
			}
			weight = 1./pow(h,power);
			numerator   += weight*obs[i];
			denominator += weight;
		}
		prediction = numerator/denominator;
	}

	/*clean-up*/
	*pprediction = prediction;
	xDelete<IssmPDouble>(x);
	xDelete<IssmPDouble>(y);
	xDelete<IssmPDouble>(obs);
}/*}}}*/
void Observations::InterpolationKriging(IssmPDouble *pprediction,IssmPDouble *perror,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata,Variogram* variogram){/*{{{*/

	/*Intermediaries*/
	int           i,j,n_obs;
	IssmPDouble   prediction,error;
	IssmPDouble  *x      = NULL;
	IssmPDouble  *y      = NULL;
	IssmPDouble  *obs    = NULL;
	IssmPDouble  *Lambda = NULL;

	/*Some checks*/
	_assert_(mindata>0 && maxdata>0);
	_assert_(pprediction && perror);

	/*Get list of observations for current point*/
	this->ObservationList(&x,&y,&obs,&n_obs,x_interp,y_interp,radius,maxdata);

	/*If we have less observations than mindata, return UNDEF*/
	if(n_obs<mindata){
		*pprediction = -999.0;
		*perror      = -999.0;
		return;
	}

	/*Allocate intermediary matrix and vectors*/
	IssmPDouble* A = xNew<IssmPDouble>((n_obs+1)*(n_obs+1));
	IssmPDouble* B = xNew<IssmPDouble>(n_obs+1);

	IssmPDouble unbias = variogram->Covariance(0.,0.);
	/*First: Create semivariogram matrix for observations*/
	for(i=0;i<n_obs;i++){
		//printf("%g %g ==> %g\n",x[i],y[i],sqrt(pow(x[i]-x_interp,2)+pow(y[i]-y_interp,2)));
		for(j=0;j<=i;j++){
			A[i*(n_obs+1)+j] = variogram->Covariance(x[i]-x[j],y[i]-y[j]);
			A[j*(n_obs+1)+i] = A[i*(n_obs+1)+j];
		}
		A[i*(n_obs+1)+n_obs] = unbias;
		//A[i*(n_obs+1)+n_obs] = 1.;
	}
	for(i=0;i<n_obs;i++) A[n_obs*(n_obs+1)+i]=unbias;
	//for(i=0;i<n_obs;i++) A[n_obs*(n_obs+1)+i]=1.;
	A[n_obs*(n_obs+1)+n_obs] = 0.;

	/*Get semivariogram vector associated to this location*/
	for(i=0;i<n_obs;i++) B[i] = variogram->Covariance(x[i]-x_interp,y[i]-y_interp);
	B[n_obs] = unbias;
	//B[n_obs] = 1.;

	/*Solve the three linear systems*/
#if _HAVE_GSL_
	DenseGslSolve(&Lambda,A,B,n_obs+1);    // Gamma^-1 Z
#else
	_error_("GSL is required");
#endif

	/*Compute predictor*/
	prediction = 0.;
	for(i=0;i<n_obs;i++) prediction += Lambda[i]*obs[i];

	/*Compute error (GSLIB p15 eq II.14)*/
	error = variogram->Covariance(0.,0.)*(1. - Lambda[n_obs]);;
	for(i=0;i<n_obs;i++) error += -Lambda[i]*B[i];

	/*clean-up*/
	*pprediction = prediction;
	*perror = error;
	xDelete<IssmPDouble>(x);
	xDelete<IssmPDouble>(y);
	xDelete<IssmPDouble>(obs);
	xDelete<IssmPDouble>(A);
	xDelete<IssmPDouble>(B);
	xDelete<IssmPDouble>(Lambda);
}/*}}}*/
void Observations::InterpolationNearestNeighbor(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius){/*{{{*/

	/*Intermediaries*/
	IssmPDouble x,y,obs;

	/*Get clostest observation*/
	this->ClosestObservation(&x,&y,&obs,x_interp,y_interp,radius);

	/*Assign output pointer*/
	*pprediction = obs;
}/*}}}*/
void Observations::InterpolationV4(IssmPDouble *pprediction,IssmPDouble x_interp,IssmPDouble y_interp,IssmPDouble radius,int mindata,int maxdata){/*{{{*/
	/* Reference:  David T. Sandwell, Biharmonic spline interpolation of GEOS-3
	 * and SEASAT altimeter data, Geophysical Research Letters, 2, 139-142,
	 * 1987.  Describes interpolation using value or gradient of value in any
	 * dimension.*/

	/*Intermediaries*/
	int         i,j,n_obs;
	IssmPDouble prediction,h;
	IssmPDouble *x       = NULL;
	IssmPDouble *y       = NULL;
	IssmPDouble *obs     = NULL;
	IssmPDouble *Green   = NULL;
	IssmPDouble *weights = NULL;
	IssmPDouble *g       = NULL;

	/*Some checks*/
	_assert_(maxdata>0);
	_assert_(pprediction);

	/*Get list of observations for current point*/
	this->ObservationList(&x,&y,&obs,&n_obs,x_interp,y_interp,radius,maxdata);

	/*If we have less observations than mindata, return UNDEF*/
	if(n_obs<mindata || n_obs<2){
		prediction = UNDEF;
	}
	else{

		/*Allocate intermediary matrix and vectors*/
		Green = xNew<IssmPDouble>(n_obs*n_obs);
		g     = xNew<IssmPDouble>(n_obs);

		/*First: distance vector*/
		for(i=0;i<n_obs;i++){
			h = sqrt( (x[i]-x_interp)*(x[i]-x_interp) + (y[i]-y_interp)*(y[i]-y_interp) );
			if(h>0){
				g[i] = h*h*(log(h)-1.);
			}
			else{
				g[i] = 0.;
			}
		}

		/*Build Green function matrix*/
		for(i=0;i<n_obs;i++){
			for(j=0;j<=i;j++){
				h = sqrt( (x[i]-x[j])*(x[i]-x[j]) + (y[i]-y[j])*(y[i]-y[j]) );
				if(h>0){
					Green[j*n_obs+i] = h*h*(log(h)-1.);
				}
				else{
					Green[j*n_obs+i] = 0.;
				}
				Green[i*n_obs+j] = Green[j*n_obs+i];
			}
			/*Zero diagonal (should be done already, but just in case)*/
			Green[i*n_obs+i] = 0.;
		}

		/*Compute weights*/
#if _HAVE_GSL_
		DenseGslSolve(&weights,Green,obs,n_obs); // Green^-1 obs
#else
		_error_("GSL is required");
#endif

		/*Interpolate*/
		prediction = 0;
		for(i=0;i<n_obs;i++) prediction += weights[i]*g[i];

	}

	/*clean-up*/
	*pprediction = prediction;
	xDelete<IssmPDouble>(x);
	xDelete<IssmPDouble>(y);
	xDelete<IssmPDouble>(obs);
	xDelete<IssmPDouble>(Green);
	xDelete<IssmPDouble>(g);
	xDelete<IssmPDouble>(weights);
}/*}}}*/
void Observations::QuadtreeColoring(IssmPDouble* A,IssmPDouble *x,IssmPDouble *y,int n){/*{{{*/

	if(this->treetype!=1) _error_("Tree type is not quadtree");
	int xi,yi,level;

	for(int i=0;i<n;i++){
		this->quadtree->IntergerCoordinates(&xi,&yi,x[i],y[i]);
		this->quadtree->QuadtreeDepth(&level,xi,yi);
		A[i]=(IssmPDouble)level;
	}

}/*}}}*/
void Observations::Variomap(IssmPDouble* gamma,IssmPDouble *x,int n){/*{{{*/

	/*Output and Intermediaries*/
	int          i,j,k;
	IssmPDouble  distance;
	Observation *observation1 = NULL;
	Observation *observation2 = NULL;

	IssmPDouble *counter = xNew<IssmPDouble>(n);
	for(j=0;j<n;j++) counter[j] = 0.0;
	for(j=0;j<n;j++) gamma[j]   = 0.0;

	for(Object* & object : this->objects){
		observation1=xDynamicCast<Observation*>(object);

		for(Object* & object : this->objects){
			observation2=xDynamicCast<Observation*>(object);

			distance=sqrt(pow(observation1->x - observation2->x,2) + pow(observation1->y - observation2->y,2));
			if(distance>x[n-1]) continue;

			int index = int(distance/(x[1]-x[0]));
			if(index>n-1) index = n-1;
			if(index<0)   index = 0;

			gamma[index]   += 1./2.*pow(observation1->value - observation2->value,2);
			counter[index] += 1.;
		}
	}

	/*Normalize semivariogram*/
	gamma[0]=0.;
	for(k=0;k<n;k++){
		if(counter[k]) gamma[k] = gamma[k]/counter[k];
	}

	/*Assign output pointer*/
	xDelete<IssmPDouble>(counter);
}/*}}}*/
