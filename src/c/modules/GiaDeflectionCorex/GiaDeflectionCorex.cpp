/*!\file GiaDeflectionCorex
 * \brief: GIA solution from Erik Ivins. 
 * Compute deflection wi from a single disk of radius re, load history hes for 
 * numtimes time steps. 
 */

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./GiaDeflectionCorex.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"

#ifdef _HAVE_AD_
void GiaDeflectionCorex( IssmDouble* pwi, IssmDouble* pdwidt, GiaDeflectionCoreArgs* arguments){
	_error_("Not compiled with AD as this function requires Fortran");
}
#else
/*External blocks: {{{*/
struct blockp{
	double pset[7];
};

struct blocko{
	double rhoi;
};

struct blockrad{
	double distrad; 
};

struct blocks{
	double aswokm_w; 
	double aswokm_dwdt; 
};

extern "C" { 
	int distme_(int* pNtime,int* pNtimp,int* pNtimm,double* time,double* bi,double* dmi,double* zhload);

	int what0_(int* piedge,int* pNtimp,int* pNtimm,double* time,double* bi,double* dmi);
	extern struct blockp blockp_;
	extern struct blocko blocko_;
	extern struct blockrad blockrad_;
	extern struct blocks blocks_;
}

/*}}}*/

void GiaDeflectionCorex( IssmDouble* pwi, IssmDouble* pdwidt, GiaDeflectionCoreArgs* arguments){

	/*Recover material parameters and loading history: see GiaDeflectionCoreArgs for more details*/
	IssmDouble  ri                        = arguments->ri;                        //radial distance from center of disk to vertex i
	IssmDouble  re                        = arguments->re;                        //radius of disk
	IssmDouble *hes                       = arguments->hes;                       //loading history (in ice thickness)
	IssmDouble *times                     = arguments->times;                     //loading history times
	int         numtimes                  = arguments->numtimes;                  //loading history length
	IssmDouble  currenttime               = arguments->currenttime;
	IssmDouble  lithosphere_shear_modulus = arguments->lithosphere_shear_modulus;
	IssmDouble  lithosphere_density       = arguments->lithosphere_density;
	IssmDouble  mantle_shear_modulus      = arguments->mantle_shear_modulus;
	IssmDouble  mantle_viscosity          = arguments->mantle_viscosity;
	IssmDouble  mantle_density            = arguments->mantle_density;
	IssmDouble  lithosphere_thickness     = arguments->lithosphere_thickness;
	IssmDouble  rho_ice                   = arguments->rho_ice;
	int         disk_id                   = arguments->idisk;
	int         iedge                     = arguments->iedge;
	IssmDouble  yts                       = arguments->yts;

	/*Modify inputs to match naruse code: */
	int Ntime=numtimes; // number of times with load history
	int Ntimm=Ntime-1; // Ntime-1 : for slope/y-cept of load segments
	int Ntimp=Ntime+1; // Ntime+1 : for evaluation time

	/*Prepare block inputs for fortran distme and what0 routines of the naruse code:*/
	/*Now, let's set pset from the data that we got in input to GiaDeflectionCorex: */
	blockp_.pset[0]=reCast<IssmPDouble>(lithosphere_thickness);
	blockp_.pset[1]=reCast<IssmPDouble>(mantle_viscosity);
	blockp_.pset[2]=reCast<IssmPDouble>(lithosphere_shear_modulus);
	blockp_.pset[3]=reCast<IssmPDouble>(mantle_shear_modulus);
	blockp_.pset[4]=reCast<IssmPDouble>(lithosphere_density);
	blockp_.pset[5]=reCast<IssmPDouble>(mantle_density);
	blockp_.pset[6]=reCast<IssmPDouble>(re);
	blocko_.rhoi=reCast<IssmPDouble>(rho_ice); 
   blockrad_.distrad=reCast<IssmPDouble>(ri)/1000.0; // in km

	/*loading history: */
	IssmPDouble* blocky_zhload=xNew<IssmPDouble>(Ntime);
	for(int i=0;i<Ntime;i++) blocky_zhload[i]=reCast<IssmPDouble>(hes[i]);

	/*times in kyr: */
	IssmPDouble* blockt_time=xNew<IssmPDouble>(Ntimp);
	for(int i=0;i<Ntimp;i++){
		blockt_time[i]=reCast<IssmPDouble>(times[i]/1000.0/yts); 
		if(i==numtimes-1) blockt_time[i]=reCast<IssmPDouble>(times[numtimes-1]/1000.0/yts); // final loading time, same as evaluation time
		if(i==numtimes)   blockt_time[i]=reCast<IssmPDouble>(times[numtimes-1]/1000.0/yts);   // evaluation time
	}

	IssmPDouble* blockt_bi=xNew<IssmPDouble>(Ntimm);
	IssmPDouble* blockt_dmi=xNew<IssmPDouble>(Ntimm);

	/*Call distme driver: */
	distme_(&Ntime,&Ntimp,&Ntimm,blockt_time,blockt_bi,blockt_dmi,blocky_zhload); 

	/*Call what0 driver: */
	what0_(&iedge,&Ntimp,&Ntimm,blockt_time,blockt_bi,blockt_dmi); 

	/*output solution: */
	*pwi=reCast<IssmDouble>(blocks_.aswokm_w);
	*pdwidt=reCast<IssmDouble>(blocks_.aswokm_dwdt);

	/*Free resources: */
	xDelete<IssmPDouble>(blockt_time);
	xDelete<IssmPDouble>(blockt_bi);
	xDelete<IssmPDouble>(blockt_dmi);
	xDelete<IssmPDouble>(blocky_zhload);

}
#endif
