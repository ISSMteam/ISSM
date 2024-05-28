/*!\file: random
 * \brief random number generating functions
 */ 

/*Headers*/
/*{{{*/
#include <stdio.h>
#include <sys/types.h>
#include <math.h>
#include <float.h>    /*  DBL_EPSILON  */
#include <cstdarg>
#include <iostream>

#include "../Matrix/matrix.h"
#include "../Exceptions/exceptions.h"
#include "../MemOps/MemOps.h"
#include "../io/io.h"
#include "./randomgenerator.h"
/*}}}*/

void univariateNormal(IssmPDouble* prand, IssmPDouble mean, IssmPDouble sdev, int seed=-1) { /*{{{*/

	/*Seed the pseudo-random number generator*/
	rnd::linear_congruential_engine randomengine;
	randomengine.seed(seed);
	/*Normal distribution*/
	rnd::normal_distribution distriNormal(mean,sdev);
	/*Assign output pointer and cleanup*/
	*prand = distriNormal.generator(randomengine);
	randomengine.free_resources();
} /*}}}*/
void multivariateNormal(IssmDouble** prand, int dim, IssmDouble mean, IssmDouble* covariancematrix, int seed=-1) { /*{{{*/
   
	IssmPDouble* sampleStandardNormal    = xNew<IssmPDouble>(dim);
   IssmDouble* sampleMultivariateNormal = xNew<IssmDouble>(dim);
   IssmDouble* Lchol                    = xNewZeroInit<IssmDouble>(dim*dim);

	/*True randomness if seed<0, otherwise random seed is fixed at seed*/
	/*Seed the pseudo-random number generator, repeatedly calling univariateNormal does not ensure randomness*/
	rnd::linear_congruential_engine randomengine;
	randomengine.seed(seed);
	/*Normal distribution*/
	rnd::normal_distribution distriNormal(0.0,1.0);
	for(int i=0;i<dim;i++){
		sampleStandardNormal[i] = distriNormal.generator(randomengine);
		//_printf_("VV i sampleStandardNormal[i]: "<<i<<"  "<<sampleStandardNormal[i]<<'\n');
	}

	/*Cholsesky decomposition of the covariance matrix*/
	CholeskyRealPositiveDefinite(Lchol,covariancematrix,dim);
   
	/*Matrix by vector multiplication*/
	for(int i=0;i<dim;i++){ 
      /*Entry-by-entry multiplication along matrix row*/
      IssmDouble sum=0.;
      for(int j=0;j<dim;j++) sum += sampleStandardNormal[j]*Lchol[i*dim+j]; 
      sampleMultivariateNormal[i] = mean+sum;
   }

   /*Assign output pointer and cleanup*/
   *prand = sampleMultivariateNormal;
   xDelete<IssmPDouble>(sampleStandardNormal);
   xDelete<IssmDouble>(Lchol);
	randomengine.free_resources();
} /*}}}*/
void multivariateNormal(IssmDouble** prand, int dim, IssmDouble* mean, IssmDouble* covariancematrix, int seed=-1) { /*{{{*/
	
	IssmPDouble* sampleStandardNormal    = xNew<IssmPDouble>(dim);
	IssmDouble* sampleMultivariateNormal = xNew<IssmDouble>(dim);
	IssmDouble* Lchol                    = xNewZeroInit<IssmDouble>(dim*dim);
	
	/*True randomness if seed<0, otherwise random seed is fixed at seed*/
	/*Seed the pseudo-random number generator, repeatedly calling univariateNormal does not ensure randomness*/
	rnd::linear_congruential_engine randomengine;
	randomengine.seed(seed);
	/*Normal distribution*/
	rnd::normal_distribution distriNormal(0.0,1.0);
	for(int i=0;i<dim;i++){
		sampleStandardNormal[i] = distriNormal.generator(randomengine);
		//_printf_("VV i sampleStandardNormal[i]: "<<i<<"  "<<sampleStandardNormal[i]<<'\n');
	}

	/*Cholsesky decomposition of the covariance matrix*/
	CholeskyRealPositiveDefinite(Lchol,covariancematrix,dim);

	/*Matrix by vector multiplication*/
	for(int i=0;i<dim;i++){
		IssmDouble sum = 0.;
      for(int j=0;j<dim;j++) sum += sampleStandardNormal[j]*Lchol[i*dim+j]; 
      sampleMultivariateNormal[i] = mean[i]+sum;
	}
   
	/*Assign output pointer and cleanup*/
	*prand = sampleMultivariateNormal;
	xDelete<IssmPDouble>(sampleStandardNormal);
	xDelete<IssmDouble>(Lchol);
	randomengine.free_resources();
} /*}}}*/



