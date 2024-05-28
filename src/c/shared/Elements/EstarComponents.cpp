#include <math.h>
#include "../Numerics/types.h"
#include "../Exceptions/exceptions.h"
#include "./elements.h"
void EstarStrainrateQuantities(IssmDouble *pepsprime_norm, IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble vmag,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble* dvmag){/*{{{*/

	/*Intermediaries*/
	IssmDouble omega[3];
	IssmDouble nrsp[3],nrsp_norm;
	IssmDouble eps[3][3];
	IssmDouble epsprime[3],epsprime_norm;

	/*Get omega, correction for rigid body rotation*/
	EstarOmega(&omega[0],vx,vy,vz,vmag,dvx,dvy,dvz,dvmag);

	/*Non-rotating shear plane*/
	nrsp[0] =  vy*omega[2] - vz*omega[1];
	nrsp[1] =  vz*omega[0] - vx*omega[2];
	nrsp[2] =  vx*omega[1] - vy*omega[0];

	/*Normalize*/
	nrsp_norm = sqrt(nrsp[0]*nrsp[0] + nrsp[1]*nrsp[1] + nrsp[2]*nrsp[2]);
	if(nrsp_norm==0){
		nrsp[0] = 0.;
		nrsp[1] = 0.;
		nrsp[2] = 0.;
	}
	else{
		nrsp[0] =nrsp[0]/nrsp_norm;
		nrsp[1] =nrsp[1]/nrsp_norm;
		nrsp[2] =nrsp[2]/nrsp_norm;
	}

	/*Build strain rate tensor*/
	eps[0][0] = dvx[0];             eps[0][1] = .5*(dvx[1]+dvy[0]); eps[0][2] = .5*(dvx[2]+dvz[0]);
	eps[1][0] = .5*(dvx[1]+dvy[0]); eps[1][1] = dvy[1];             eps[1][2] = .5*(dvy[2]+dvz[1]);
	eps[2][0] = .5*(dvx[2]+dvz[0]); eps[2][1] = .5*(dvy[2]+dvz[1]); eps[2][2] = dvz[2];

	/*Compute the shear strain rate on the non rotating shear plane*/
	epsprime[0]=0.;
	epsprime[1]=0.;
	epsprime[2]=0.;
	/*term #1: eps'.n */
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			epsprime[i] += eps[i][j]*nrsp[j];
		}
	}
	/*term #2: ((eps'.n).n)n */
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				epsprime[j] += -nrsp[i]*eps[i][k]*nrsp[k]*nrsp[j];
			}
		}
	}
	/*term #3: ((eps'.n).omega)omega */
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				epsprime[j] += -nrsp[i]*eps[i][k]*omega[k]*omega[j];
			}
		}
	}

	/*Get norm of epsprime*/
	epsprime_norm = sqrt(epsprime[0]*epsprime[0] + epsprime[1]*epsprime[1] + epsprime[2]*epsprime[2]);

	/*Assign output pointers*/
	*pepsprime_norm=epsprime_norm;
}/*}}}*/
void EstarOmega(IssmDouble* omega,IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble vmag,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz, IssmDouble* dvmag){/*{{{*/

	/*Intermediaries*/
	IssmDouble omega_norm;
	IssmDouble omega_rigid[3];

	/*Create vorticity vector*/
	_assert_(dvx && dvy && dvz && dvmag);
	if(vmag<1e-12)vmag=1e-12;

	/*Create vorticity vector, corrected for rigid body rotation
	 * \overline{\omega} =\omega - \omega_rigid
	 *                   =\nabla\times{\bf v} - 2*U*\kappa*\hat{b};
	 *                   =\nabla\times{\bf v} - (2*{\bf v}\times(({\bf v}\cdot\nabla)*{\bf v}))/U^2
	 * check the magnitude of the second term -- if it is small, then the two
	 * vorticities (omega and first term in omega) are approx. equal
	 *
	 * */
	omega_rigid[0] = 2*(vy*(vx*dvz[0]+vy*dvz[1]+vz*dvz[2]) - vz*(vx*dvy[0]+vy*dvy[1]+vz*dvy[2]))/(vmag*vmag);
	omega_rigid[1] = 2*(vz*(vx*dvx[0]+vy*dvx[1]+vz*dvx[2]) - vx*(vx*dvz[0]+vy*dvz[1]+vz*dvz[2]))/(vmag*vmag);
	omega_rigid[2] = 2*(vx*(vx*dvy[0]+vy*dvy[1]+vz*dvy[2]) - vy*(vx*dvx[0]+vy*dvx[1]+vz*dvx[2]))/(vmag*vmag);

	omega[0] = (dvz[1] - dvy[2]) - omega_rigid[0];
	omega[1] = (dvx[2] - dvz[0]) - omega_rigid[1];
	omega[2] = (dvy[0] - dvx[1]) - omega_rigid[2];

	/*Take out vorticity component aligned with the velocity*/
	IssmDouble wdotv = vx/vmag*omega[0] + vy/vmag*omega[1] + vz/vmag*omega[2];
	omega[0] = omega[0] - wdotv*vx/vmag;
	omega[1] = omega[1] - wdotv*vy/vmag;
	omega[2] = omega[2] - wdotv*vz/vmag;

	/*Normalize*/
	omega_norm = sqrt(omega[0]*omega[0] + omega[1]*omega[1] + omega[2]*omega[2]);
	if(omega_norm==0){
		omega[0] = 0.;
		omega[1] = 0.;
		omega[2] = 0.;
	}
	else{
		omega[0] =omega[0]/omega_norm;
		omega[1] =omega[1]/omega_norm;
		omega[2] =omega[2]/omega_norm;
	}

}/*}}}*/
IssmDouble EstarLambdaS(IssmDouble epseff, IssmDouble epsprime_norm){/*{{{*/
	IssmDouble lambdas;

	_assert_(epsprime_norm>=0.); 
	if(epseff==0.){
		lambdas=0.;
	}
	else{
		lambdas=sqrt(epsprime_norm*epsprime_norm/(epseff*epseff));
	}
	return lambdas; 
}/*}}}*/
