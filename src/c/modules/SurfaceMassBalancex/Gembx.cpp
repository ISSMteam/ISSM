/*!\file GEMB module from Alex Gardner.
 * \brief: calculates SMB 
 */

#include "./SurfaceMassBalancex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../modules.h"
#include "../../classes/Inputs/TransientInput.h"

const double Pi = 3.141592653589793;
const double CtoK = 273.15;             // Kelvin to Celcius conversion/ice melt. point T in K
const double dts = 86400.0;              // Number of seconds in a day

/* Tolerances have to be defined for if loops involving some specific values
	like densitiy of ice or melting point temp because values are "rounded"
	(e.g rho ice = 909.99... instead of 910.0) and we don't always go into the right loop */
const double Ttol = 1e-10;
const double Dtol = 1e-11;
const double Gdntol = 1e-10;
const double Wtol = 1e-13;
const double Ptol = 1e-6;

const double CI = 2102.0;                       // heat capacity of snow/ice (J kg-1 k-1)
const double LF = 0.3345e6;             // latent heat of fusion (J kg-1)
const double LV = 2.495e6;               // latent heat of vaporization (J kg-1)
const double LS = 2.8295e6;             // latent heat of sublimation (J kg-1)
const double SB = 5.67e-8;                // Stefan-Boltzmann constant (W m-2 K-4)
const double CA = 1005.0;                    // heat capacity of air (J kg-1 K-1)
const double R = 8.314;                      // gas constant (J mol-1 K-1)

const double Delflag=-99999;

void Gembx(FemModel* femmodel){  /*{{{*/

	int        count=0;
	int        steps=0;
	IssmDouble time,dt,finaltime;
	IssmDouble timeclim=0.0;
	IssmDouble t,smb_dt;
   IssmDouble delta;
	bool       linear_interp=true;

	femmodel->parameters->FindParam(&linear_interp,TimesteppingInterpForcingEnum); /*is interpolation requested*/
	femmodel->parameters->FindParam(&time,TimeEnum);                        /*transient core time at which we run the smb core*/
   femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);          /*transient core time step*/
	femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
   femmodel->parameters->FindParam(&smb_dt,SmbDtEnum);                     /*time period for the smb solution,  usually smaller than the glaciological dt*/

	//before starting loop, realize that the transient core runs this smb_core at time = time +deltaT.
	//go back to time - deltaT:
	time-=dt;

	IssmDouble timeinputs = time;

	steps=0;
	for (t=time;t<=time+dt-smb_dt;t=t+smb_dt){
		steps=steps+1;
	}

	/*Start loop: */
	count=1;
	for (t=time;t<=time+dt-smb_dt;t=t+smb_dt){

		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);

			timeclim=time;
			if (linear_interp) timeinputs = t-time+timeclim;
			else timeinputs = t-time+timeclim+smb_dt/2;
			element->SmbGemb(timeinputs,count,steps);
		}
		count=count+1;
	}

} /*}}}*/
void GembgridInitialize(IssmDouble** pdz, int* psize, IssmDouble zTop, IssmDouble dzTop, IssmDouble zMax, IssmDouble zY){ /*{{{*/

	/* This file sets up the initial grid spacing and total grid depth.  The
	grid structure is set as constant grid length 'dzTop' for the top
	'zTop' meters of the model grid. Bellow 'zTop' the gid length increases
	linearly with depth */

	/*intermediary:*/
	IssmDouble dgpTop=0.0;
	int gpTop=0;
	int gpBottom=0;
	int i=0;
	IssmDouble gp0=0.0;
	IssmDouble z0=0.0;
	IssmDouble* dzT=NULL;
	IssmDouble* dzB=NULL;

	/*output: */
	IssmDouble* dz=NULL;

	//----------------------Calculate Grid Lengths------------------------------
	//calculate number of top grid points
	dgpTop = zTop/dzTop;

	//check to see if the top grid cell structure length (dzTop) goes evenly 
	//into specified top structure depth (zTop). Also make sure top grid cell
	//structure length (dzTop) is greater than 5 cm
	#ifndef _HAVE_AD_  //avoid the round operation check!
	if (dgpTop != round(dgpTop)){ 
		_error_("top grid cell structure length does not go evenly into specified top structure depth, adjust dzTop or zTop\n");
	}
	#endif
	if(dzTop < 0.05-Dtol){
		_printf_("initial top grid cell length (dzTop) is < 0.05 m\n");
	}
	gpTop=reCast<int,IssmDouble>(dgpTop);

	//initialize top grid depth vector
	dzT = xNew<IssmDouble>(gpTop); 
	for (i=0;i<gpTop-Dtol;i++)dzT[i]=dzTop;

	//bottom grid cell depth = x*zY^(cells from to structure)
	//figure out the number of grid points in the bottom vector (not known a priori)
	gp0 = dzTop;
	z0 = zTop;
	gpBottom = 0;
	while (zMax > z0+Dtol){
		gp0= gp0 * zY;
		z0 = z0 + gp0;
		gpBottom++;
	}
	//initialize bottom vectors
	dzB = xNewZeroInit<IssmDouble>(gpBottom);
	gp0 = dzTop;
	z0 = zTop;
	for(i=0;i<gpBottom;i++){
		gp0=gp0*zY;
		dzB[i]=gp0;
	}

	//combine top and bottom dz vectors
	dz = xNew<IssmDouble>(gpTop+gpBottom);
	for(i=0;i<gpTop-Dtol;i++){
		dz[i]=dzT[i];
	}
	for(i=0;i<gpBottom;i++){
		dz[gpTop+i]=dzB[i];
	}

	/*Free resouces:*/
	xDelete(dzT);
	xDelete(dzB);

	//---------NEED TO IMPLEMENT A PROPER GRID STRECHING ALGORITHM------------

	/*assign ouput pointers: */
	*pdz=dz; 
	*psize=gpTop+gpBottom;
} /*}}}*/ 
IssmDouble Marbouty(IssmDouble T, IssmDouble d, IssmDouble dT){ /*{{{*/

	// calculates grain growth according to Fig. 9 of Marbouty, 1980
	// ------NO GRAIN GROWTH FOR d > 400 kg m-3 because H is set to zero------
	// ---------------this is a major limitation of the model-------------------

	// initialize
	IssmDouble F = 0.0, H=0.0, G=0.0;
	const IssmDouble E = 0.09;        //[mm d-1] model time growth constant E
	// convert T from K to degC
	T = T - CtoK;
	// convert dT from degC/m to degC/cm
	dT = dT/100.0;

	// temperature coefficient F
	if(T> -6.0+Ttol) F =  0.7 + ((T/-6.0) * 0.3);
	if(T<= -6.0+Ttol && T> -22.0+Ttol) F =  1.0 - ((T+6.0)/-16.0 * 0.8);
	if(T<= -22.0+Ttol && T> -40.0+Ttol) F =  0.2 - ((T+22.0)/-18.0 * 0.2);

	// density coefficient F
	if(d< 150.0-Dtol) H=1.0;

	if(d>= 150.0-Dtol && d <400.0-Dtol) H = 1.0 - ((d-150.0)/250.0);

	// temperature gradient coefficient G
	if(dT >= 0.16-Ttol && dT < 0.25-Ttol) G = ((dT - 0.16)/0.09) * 0.1;
	if(dT >= 0.25-Ttol && dT < 0.4-Ttol)  G = 0.1 + (((dT - 0.25)/0.15) * 0.57);
	if(dT >= 0.4-Ttol && dT < 0.5-Ttol)  G = 0.67 + (((dT - 0.4)/0.1) * 0.23);
	if(dT >= 0.5-Ttol && dT < 0.7-Ttol)  G = 0.9 + (((dT - 0.5)/0.2) * 0.1);
	if(dT >= 0.7-Ttol)  G = 1.0;

	// grouped coefficient Q
	return F*H*G*E;

} /*}}}*/
IssmDouble gardnerAlb(IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, IssmDouble dPHC, int m){ /*{{{*/
	//gardnerAlb(S1, c1, SZA, t, z1, S2, c2)
	//This is an implementation of the snow and ice broadband albedo
	//  parameterization developed by Alex Gardner.
	//Created By: Alex S. Gardner, Jet Propulsion Laboratory [alex.s.gardner@jpl.nasa.gov]
	//  Last Modified: June, 2014
	//Full Reference: Gardner, A. S., and Sharp, M. J.: A review of snow and
	//  ice albedo and the development of a new physically based broadband albedo
	//  parameterization, J. Geophys. Res., 115, F01009, 10.1029/2009jf001444,
	//  2010.

	//INPUTS
	// ONE LAYER
	//  - S1    : specific surface area of the snow or ice [cm^2 g-1]
	//  - c1    : concentration of light absorbing carbon  [ppm1]
	//  - SZA   : solar zenith angle of the incident radiation [deg]
	//  - t     : cloud optical thickness
	// TWO LAYER
	//  - z1    : depth of snow suface layer [mm w.e.]
	//  - S2    : specific surface area of bottom ice layer [cm^2 g-1]
	//  - c2    : concentration of light absorbing carbon of bottom ice
	//             layer [ppm1]
	IssmDouble c1=clabSnow;
	IssmDouble c2=clabIce;
	IssmDouble t=COT;
	IssmDouble a=0.0;

	//Single layer albedo parameterization
	//convert effective radius to specific surface area [cm2 g-1]
	IssmDouble S1 = 3.0 / (0.091 * re[0]);

	//effective solar zenith angle
	IssmDouble x = min(pow(t/(3.0*cos(Pi*SZA/180.0)),0.5), 1.0);
	IssmDouble u = 0.64*x + (1.0-x)*cos(Pi*SZA/180.0);

	// pure snow albedo
	IssmDouble as = 1.48 - pow(S1,-0.07);

	//change in pure snow albedo due to soot loading
	IssmDouble dac = max(0.04 - as, pow(-c1,0.55)/(0.16 + 0.6*pow(S1,0.5) + (1.8*pow(c1,0.6))*pow(S1,-0.25)));

	//Two layer albedo parameterization
	//  do two layer calculation if there is more than 1 layer
	IssmDouble z1=0.0;
	int lice=0;
	for(int l=0;(l<m && d[l]<dPHC-Dtol);l++){
		z1=z1+dz[l]*d[l]; //mm
		lice=l+1;
	}
	if (m>0 & lice<m & z1 > Dtol){
		// determine albedo values for bottom layer
		IssmDouble S2 = 3.0 / (0.091 * re[lice]);

		// pure snow albedo
		IssmDouble as2 = 1.48 - pow(S2,-0.07);

		// change in pure snow albedo due to soot loading
		IssmDouble dac2 = max(0.04 - as2, pow(-c2,0.55)/(0.16 + 0.6*pow(S2,0.5) + (1.8*pow(c2,0.6))*pow(S2,-0.25)));

		// determine the effective change due to finite depth and soot loading
		IssmDouble A = min(1.0, (2.1 * pow(z1,1.35*(1.0-as) - 0.1*c1 - 0.13)));

		dac =  (as2 + dac2 - as) + A*((as + dac) - (as2 + dac2));
	}

	// change in albedo due to solar zenith angle
	IssmDouble dasz = 0.53*as*(1.0 - (as + dac))*pow(1.0 - u,1.2);

	// change in albedo due to cloud (apart from change in diffuse fraction)
	IssmDouble dat = (0.1*t*pow(as + dac,1.3)) / (pow(1.0 + 1.5*t,as));

	// Broadband albedo
	a = as + dac + dasz + dat;

	return a;
}   /*}}}*/
void grainGrowth(IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, IssmDouble* T,IssmDouble* dz,IssmDouble* d, IssmDouble* W,IssmDouble smb_dt,int m,int aIdx,int sid){ /*{{{*/

	/*Created by: Alex S. Gardner, University of Alberta

	 *Description*: models the effective snow grain size

	 *Reference:*
	 DENDRITIC SNOW METAMORPHISM:
	 Brun, E., P. David, M. Sudul, and G. Brunot, 1992: A numerical model to
	 simulate snow-cover stratigraphy for operational avalanche forecasting.
	 Journal of Glaciology, 38, 13-22.

	 NONDENDRITIC SNOW METAMORPHISM:
	 Dry snow metamorphism:
	 Marbouty, D., 1980: An experimental study of temperature-gradient
	 metamorphism. Journal of Glaciology, 26, 303-312.

	 WET SNOW METAMORPHISM:
	 Brun, E., 1989: Investigation on wet-snow metamorphism in respect of
	 liquid-water content. Annals of Glaciology, 13, 22-26.

	 INPUTS
	 * T: grid cell temperature [K]
	 * dz: grid cell depth [m]
	 * d: grid cell density [kg m-3]
	 * W: water content [kg/m^2]
	 * re: effective grain size [mm]
	 * gdn: grain dentricity
	 * gsp: grain sphericity
	 * dt: time step of input data [s]

	 OUTPUTS
	 * re: effective grain size [mm]
	 * gdn: grain dentricity
	 * gsp: grain sphericity*/

	/*intermediary: */
	IssmDouble  dt=0.0;
	IssmDouble  Ti=0.0;
	IssmDouble* gsz=NULL;
	IssmDouble* dT=NULL;
	IssmDouble* zGPC=NULL;
	IssmDouble* lwc=NULL;
	IssmDouble  Q=0.0;

	/*output: */
	IssmDouble* re=NULL;
	IssmDouble* gdn=NULL;
	IssmDouble* gsp=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   grain growth module\n");

	/*Recover pointers: */
	re=*pre;
	gdn=*pgdn;
	gsp=*pgsp;

	/*only when aIdx = 1 or 2 do we run grainGrowth: */
	if(aIdx!=1 && aIdx!=2){
		/*come out as we came in: */
		return;
	}

	/*Figure out grain size from effective grain radius: */
	gsz=xNew<IssmDouble>(m); for(int i=0;i<m;i++)gsz[i]=re[i]*2.0;

	/*Convert dt from seconds to day: */
	dt=smb_dt/dts;

	/*Determine liquid-water content in percentage: */
	lwc=xNew<IssmDouble>(m); for(int i=0;i<m;i++)lwc[i]= W[i] / (d[i]*dz[i])*100.0;

	//set maximum water content by mass to 9 percent (Brun, 1980)
	for(int i=0;i<m;i++)if(lwc[i]>9.0+Wtol) lwc[i]=9.0;

	/* Calculate temperature gradiant across grid cells 
	 * Returns the average gradient across the upper and lower grid cell */

	//initialize
	dT=xNewZeroInit<IssmDouble>(m); 

	//depth of grid point center from surface
	zGPC=xNewZeroInit<IssmDouble>(m);  
	for(int i=0;i<m;i++){
		for (int j=0;j<=i;j++) zGPC[i]+=dz[j];
		zGPC[i]-=(dz[i]/2.0);
	}

	// Take forward differences on left and right edges
	if(m>2){
		dT[0] = (T[2] - T[0])/(zGPC[2]-zGPC[0]);
		dT[m-1] = (T[m-1] - T[m-3])/(zGPC[m-1]-zGPC[m-3]);
	}
	else if(m>1){
		dT[0] = (T[1] - T[0])/(zGPC[1]-zGPC[0]);
		dT[m-1] = (T[m-1] - T[m-2])/(zGPC[m-1]-zGPC[m-2]);
	}

	//Take centered differences on interior points
	for(int i=1;i<m-1;i++) dT[i] = (T[i+1]-T[i-1])/(zGPC[i+1]-zGPC[i-1]);

	// take absolute value of temperature gradient
	for(int i=0;i<m;i++)dT[i]=fabs(dT[i]);

	/*Snow metamorphism. Depends on value of snow dendricity and wetness of the snowpack: */
	for(int i=0;i<m;i++){

		// T for this layer
		Ti = T[i];

		if (gdn[i]>0.0+Gdntol){

			if(W[i]<=0.0+Wtol){
				//_printf_("Dendritic dry snow metamorphism\n");
				//index for dentricity > 0 and T gradients < 5 degC m-1 and >= 5 degC m-1
				if(fabs(dT[i])<=5.0+Ttol){
					//determine coefficients
					IssmDouble A = - 2e8 * exp(-6e3 / Ti) * dt;
					IssmDouble B = 1e9 * exp(-6e3 / Ti) * dt;
					//new dentricity and sphericity for dT < 5 degC m-1
					gdn[i] += A;
					gsp[i] += B;
				}
				else{
					// new dendricity and sphericity for dT >= 5 degC m-1

					//determine coefficients
					IssmDouble C = (-2e8 * exp(-6e3 / Ti) * dt) * pow(fabs(dT[i]),.4);
					gdn[i] +=C;
					gsp[i] +=C;
				}
			}
			else{
				// wet snow metamorphism
				//_printf_("Dendritic wet snow metamorphism\n");

				//determine coefficient
				IssmDouble D = (1.0/16.0) * pow(lwc[i],3.0) * dt;

				// new dendricity for wet snow
				gdn[i] -= D;
				// new sphericity for wet snow
				gsp[i] += D;
			}
			// dendricity and sphericity can not be > 1 or < 0
         if (gdn[i]<=0.0+Gdntol)gdn[i]=0.0;
         if (gsp[i]<=0.0+Gdntol)gsp[i]=0.0;
         if (gdn[i]>=1.0-Gdntol)gdn[i]=1.0;
         if (gsp[i]>=1.0-Gdntol)gsp[i]=1.0;

         // determine new grain size (mm)
			gsz[i] = max(1e-1*(gdn[i]/.99+(1.0-1.0*gdn[i]/.99)*(gsp[i]/.99*3.0+(1.0-gsp[i]/.99)*4.0)),Gdntol*2.0);

		}
		else{

			//When wet-snow grains (class 6) are submitted to a
			// temperature gradient higher than 5 degC m-1, their sphericity
			// decreases according to Equations (4). When sphericity
			// reaches 0, their size increases according to the functions
			// determined by Marbouty. (Brun et al., 1992)
			if(gsp[i]>0.0+Gdntol && gsp[i]<1.0-Gdntol){

				IssmDouble F = 0.0;

				if (fabs(dT[i])>5.0+Ttol){
					F = (-2e8 * exp(-6e3 / Ti) * dt) * pow(fabs(dT[i]),.4);
				}
				else if (W[i]>0.0+Wtol){
					F = (1.0/16.0) * pow(lwc[i],3.0) * dt;
				}
				else{
					F = 1e9 * exp(-6e3 / Ti) * dt;
				}
				gsp[i] +=F;

			}
			if (gsp[i]<=0.0+Gdntol)gsp[i]=0.0;
			if (gsp[i]>=1.0-Gdntol)gsp[i]=1.0;

			/*Dry snow metamorphism (Marbouty, 1980) grouped model coefficents
			 *from Marbouty, 1980: Figure 9*/
			if(W[i]<=0.0+Wtol || (gsp[i]<=0.0+Gdntol && fabs(dT[i])>5.0+Ttol)){
				//_printf_("Nondendritic snow metamorphism\n");
				Q = Marbouty(Ti,d[i],dT[i]);

				// calculate grain growth
				gsz[i] += (Q*dt);
			}
			//Wet snow metamorphism (Brun, 1989)
			else{
				//_printf_("Nondendritic wet snow metamorphism\n");
				//wet rate of change coefficient
				IssmDouble E = (1.28e-8 + 4.22e-10 * pow(lwc[i],3.0))* (dt *dts);   // [mm^3 s^-1]

				// calculate change in grain volume and convert to grain size
				gsz[i] = 2.0 * pow(3.0/(Pi * 4.0)*((4.0/ 3.0)*Pi*pow(gsz[i]/2.0,3.0) + E),1.0/3.0);
			}

			// grains with sphericity == 1 can not have grain sizes > 2 mm (Brun, 1992)
			if (fabs(gsp[i]-1.0)<Wtol && gsz[i]>2.0-Wtol) gsz[i]=2.0;

			// grains with sphericity == 0 can not have grain sizes > 5 mm (Brun, 1992)
			if (fabs(gsp[i]-1.0)>=Wtol && gsz[i]>5.0-Wtol) gsz[i]=5.0;
		}

		//convert grain size back to effective grain radius:
		re[i]=gsz[i]/2.0;
	}

	/*Free resources:*/
	xDelete<IssmDouble>(gsz);
	xDelete<IssmDouble>(dT);
	xDelete<IssmDouble>(zGPC);
	xDelete<IssmDouble>(lwc);

	/*Assign output pointers:*/
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;

}  /*}}}*/
void albedo(IssmDouble** pa, IssmDouble** padiff, int aIdx, IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble cldFrac, IssmDouble aIce, IssmDouble aSnow, IssmDouble aValue, IssmDouble adThresh, IssmDouble* TK, IssmDouble* W, IssmDouble P, IssmDouble EC, IssmDouble Msurf, IssmDouble clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, IssmDouble t0wet, IssmDouble t0dry, IssmDouble K, IssmDouble dt, IssmDouble dIce, int m,int sid) { /*{{{*/

	//// Calculates Snow, firn and ice albedo as a function of:
	//   0 : direct input from aValue parameter
	//   1 : effective grain radius (Gardner & Sharp, 2009)
	//   2 : effective grain radius (Brun et al., 1992, Lefebre et al., 2003)
	//   3 : density and cloud amount (Greuell & Konzelmann, 1994)
	//   4 : exponential time decay & wetness (Bougamont & Bamber, 2005)

	//// Inputs
	// aIdx      = albedo method to use

	// Method 0
	//  aValue   = direct input value for albedo, override all changes to albedo

	// adThresh
	//  Apply below method to all areas with densities below this value, 
	//  or else apply direct input value, allowing albedo to be altered.  
	//  Default value is rho water (1023 kg m-3).

	// Methods 1 & 2
	//   re      = surface effective grain radius [mm]
	// Method 1, optional
	//  clabSnow = concentration of light absorbing carbon  [ppm1], default 0
	//  SZA      = solar zenith angle of the incident radiation [deg], default 0
	//  COT      = cloud optical thickness, default 0
	//  For TWO LAYER
	//  clabIce  = concentration of light absorbing carbon of first ice layer [ppm1], default 0

	// Method 3
	//   d       = snow surface density [kg m-3]
	//   n       = cloud amount
	//   aIce    = albedo of ice
	//   aSnow   = albedo of fresh snow

	// Method 4
	//   aIce    = albedo of ice
	//   aSnow   = albedo of fresh snow
	//   a       = grid cell albedo from prevous time step;
	//   T       = grid cell temperature [k]
	//   W       = pore water [kg]
	//   P       = precipitation [mm w.e.] or [kg m-3]
	//   EC      = surface evaporation (-) condensation (+) [kg m-2]
	//   t0wet   = time scale for wet snow (15-21.9) [d]
	//   t0dry   = warm snow timescale [15] [d]
	//   K       = time scale temperature coef. (7) [d]
	//   dt      = time step of input data [s]

	//// Output
	//   a       = grid cell albedo 

	//// Usage 
	// Method 1
	// a = albedo(1, 0.1); 

	// Method 4
	// a = albedo(4, [], [], [], 0.48, 0.85, [0.8 0.5 ... 0.48], ...
	//   [273 272.5 ... 265], [0 0.001 ... 0], 0, 0.01, 15, 15, 7, 3600)

	/*output: */
	IssmDouble* a=NULL;
	IssmDouble* adiff=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   albedo module\n");

	/*Recover pointers: */
	a=*pa;
	adiff=*padiff;

	//some constants:
	const IssmDouble dSnow = 300.0;   // density of fresh snow [kg m-3]       
	const IssmDouble dPHC = 830.0;  //Pore closeoff density
	const IssmDouble ai_max = 0.58;  //maximum ice albedo, from Lefebre,2003
	const IssmDouble ai_min = aIce;  //minimum ice albedo
	const IssmDouble as_min = 0.65;  //minimum snow albedo, from Alexander 2014

	if(aIdx==0 || (adThresh - d[0])<Dtol){
		a[0] = aValue;
	}
	else{
		if(aIdx==1){ 
			//function of effective grain radius
			// clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, int m
			a[0]=gardnerAlb(re, dz, d, clabSnow, clabIce, SZA, COT, dPHC, m);
			adiff[0]=gardnerAlb(re, dz, d, clabSnow, clabIce, 50.0, COT, dPHC, m);
		}
		else if(aIdx==2){

			// Spectral fractions  (Lefebre et al., 2003)
			// [0.3-0.8um 0.8-1.5um 1.5-2.8um]

			IssmDouble sF[3] = {0.606, 0.301, 0.093};

			// convert effective radius to grain size in meters
			IssmDouble gsz = (re[0] * 2.0) / 1000.0;

			// spectral range:
			// 0.3 - 0.8um
			IssmDouble a0 = min(0.98, 0.95 - 1.58 *pow(gsz,0.5));
			// 0.8 - 1.5um
			IssmDouble a1 = max(0., 0.95 - 15.4 *pow(gsz,0.5));
			// 1.5 - 2.8um
			IssmDouble a2 = max(0.127, 0.88 + 346.3*gsz - 32.31*pow(gsz,0.5));

			// broadband surface albedo
			a[0] = sF[0]*a0 + sF[1]*a1 + sF[2]*a2;

		}
		else if(aIdx==3){

			// a as a function of density

			// calculate albedo
			a[0] = aIce + (d[0] - dIce)*(aSnow - aIce) / (dSnow - dIce) + (0.05 * (cldFrac - 0.5));
		}
		else if(aIdx==4){

			// exponential time decay & wetness

			// change in albedo with time:
			//   (d_a) = (a - a_old)/(t0)
			// where: t0 = timescale for albedo decay

			dt = dt / dts;    // convert from [s] to [d]

			// initialize variables
			IssmDouble* t0=xNew<IssmDouble>(m);
			IssmDouble* T=xNew<IssmDouble>(m);
			IssmDouble* t0warm=xNew<IssmDouble>(m);
			IssmDouble* d_a=xNew<IssmDouble>(m);

			// specify constants
			// a_wet = 0.15;        // water albedo (0.15)
			// a_new = aSnow        // new snow albedo (0.64 - 0.89)
			// a_old = aIce;        // old snow/ice albedo (0.27-0.53)
			// t0_wet = t0wet;      // time scale for wet snow (15-21.9) [d]
			// t0_dry = t0dry;      // warm snow timescale [15] [d]
			// K = 7                // time scale temperature coef. (7) [d]
			// W0 = 300;            // 200 - 600 [mm]
			const IssmDouble z_snow = 15.0;            // 16 - 32 [mm]

			// determine timescale for albedo decay
			for(int i=0;i<m;i++)if(W[i]>0.0+Wtol)t0[i]=t0wet; // wet snow timescale
			for(int i=0;i<m;i++)T[i]=TK[i] - CtoK; // change T from K to degC
			for(int i=0;i<m;i++) t0warm[i]= fabs(T[i]) * K + t0dry; //// 'warm' snow timescale
			for(int i=0;i<m;i++)if(fabs(W[i])<Wtol && T[i]>=-10.0-Ttol)t0[i]= t0warm[i];
			for(int i=0;i<m;i++)if(T[i]<-10.0-Ttol) t0[i] =  10.0 * K + t0dry; // 'cold' snow timescale

			// calculate new albedo
			for(int i=0;i<m;i++)d_a[i] = (a[i] - aIce) / t0[i] * dt;           // change in albedo
			for(int i=0;i<m;i++)a[i] -= d_a[i];                            // new albedo

			// modification of albedo due to thin layer of snow or solid
			// condensation (deposition) at the surface 

			// check if condensation occurs & if it is deposited in solid phase
			if ( EC > 0.0 + Dtol && T[0] < 0.0-Ttol) P = P + (EC/dSnow) * 1000.0;  // add cond to precip [mm]

			a[0] = aSnow - (aSnow - a[0]) * exp(-P/z_snow);

			//----------THIS NEEDS TO BE IMPLEMENTED AT A LATER DATE------------
			// modification of albedo due to thin layer of water on the surface
			// a_surf = a_wet - (a_wet - a_surf) * exp(-W_surf/W0);

			/*Free resources:*/
			xDelete<IssmDouble>(t0);
			xDelete<IssmDouble>(T);
			xDelete<IssmDouble>(t0warm);
			xDelete<IssmDouble>(d_a);

		}
		else _error_("albedo method switch should range from 0 to 4!");

		//If we do not have fresh snow
		if (aIdx<3 && aIdx>0 && (adThresh - d[0])>=Dtol){
			// In a snow layer < 10cm, account for mix of ice and snow,
			// after P. Alexander et al., 2014
			IssmDouble depthsnow=0.0;
			IssmDouble aice=0.0;
			int lice=0;
			for(int l=0;(l<m && d[l]<dPHC-Dtol);l++){
				depthsnow=depthsnow+dz[l];
				lice=l+1;
			}
			if (depthsnow<=0.1+Dtol && lice<m && d[lice]>=dPHC-Dtol){
				aice = ai_max + (as_min - ai_max)*(d[lice]-dIce)/(dPHC-dIce);
				a[0]= aice + max(a[0]-aice,0.0)*(depthsnow/0.1);
			}

			if (d[0]>=dPHC-Dtol){
				if (d[0]<dIce-Dtol){ //For continuity of albedo in firn i.e. P. Alexander et al., 2014

					//ai=ai_max + (as_min - ai_max)*(dI-dIce)/(dPHC-dIce);
					//dPHC is pore close off (830 kg m^-3)
					//dI is density of the upper firn layer

					a[0] = ai_max + (as_min - ai_max)*(d[0]-dIce)/(dPHC-dIce);

				}
				else{ //surface layer is density of ice

					//When density is > dIce (typically 910 kg m^-3, 920 is used by Alexander in MAR),
					//ai=ai_min + (ai_max - ai_min)*e^(-1*(Msw(t)/K))
					//K is a scale factor (set to 200 kg m^-2)
					//Msw(t) is the time-dependent accumulated amount of excessive surface meltwater
					//  before run-off in kg m^-2 (melt per GEMB timestep, i.e. 3 hourly)
					IssmDouble M = Msurf+W[0];
					a[0]=max(ai_min + (ai_max - ai_min)*exp(-1.0*(M/200.0)), ai_min);

				}
			}
		}
	}

	// Check for erroneous values
	if (a[0] > 1.0+Ttol) _printf_("albedo > 1.0\n");
	else if (a[0] < 0.0-Dtol) _printf_("albedo is negative\n");
	else if (xIsNan(a[0])) _error_("albedo == NAN\n");

	/*Assign output pointers:*/
	*pa=a;
	*padiff=adiff;

}  /*}}}*/
void thermo(IssmDouble* pshf, IssmDouble* plhf, IssmDouble* pEC, IssmDouble** pT, IssmDouble* pulwrf, IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble* swf, IssmDouble dlwrf, IssmDouble Ta, IssmDouble V, IssmDouble eAir, IssmDouble pAir, int tcIdx, int eIdx, IssmDouble teValue, IssmDouble dulwrfValue, IssmDouble teThresh, IssmDouble Ws, IssmDouble dt0, IssmDouble dzMin, int m, IssmDouble Vz, IssmDouble Tz, IssmDouble thermo_scaling, IssmDouble dIce, int sid, bool isconstrainsurfaceT, bool isdeltaLWup) { /*{{{*/

	/* ENGLACIAL THERMODYNAMICS*/

	/* Description: 
	   computes new temperature profile accounting for energy absorption and 
	   thermal diffusion.*/

	// INPUTS
	//  T: grid cell temperature [k]
	//  dz: grid cell depth [m]
	//  d: grid cell density [kg m-3]
	//  swf: shortwave radiation fluxes [W m-2]
	//  dlwrf: downward longwave radiation fluxes [W m-2]
	//  Ta: 2 m air temperature
	//  V:  wind velocity [m s-1]
	//  eAir: screen level vapor pressure [Pa]
	//  Ws: surface water content [kg]
	//  dt0: time step of input data [s]
	//  elev: surface elevation [m a.s.l.] 
	//  Vz: air temperature height above surface [m]
	//  Tz: wind height above surface [m]
	//  thermo_scaling: scaling factor to multiply the thermal diffusion timestep (delta t) 

	// OUTPUTS
	// T: grid cell temperature [k]
	// EC: evaporation/condensation [kg]
	// ulwrf: upward longwave radiation flux [W m-2]

	/*intermediary: */
	IssmDouble* K = NULL;
	IssmDouble* KU = NULL;
	IssmDouble* KD = NULL;
	IssmDouble* KP = NULL;
	IssmDouble* Au = NULL;
	IssmDouble* Ad = NULL;
	IssmDouble* Ap = NULL;
	IssmDouble* Nu = NULL;
	IssmDouble* Nd = NULL;
	IssmDouble* Np = NULL;
	IssmDouble* dzU = NULL;
	IssmDouble* dzD = NULL;
	IssmDouble* sw = NULL;
	IssmDouble* dT_sw = NULL;
	IssmDouble* T0 = NULL;
	IssmDouble* Tu = NULL;
	IssmDouble* Td = NULL;

	IssmDouble z0=0.0;	
	IssmDouble zT=0.0;
	IssmDouble zQ=0.0;
	IssmDouble zratio=1.0;
	IssmDouble dt=0.0;
	IssmDouble max_fdt=0.0;
	IssmDouble Ts=0.0;
	IssmDouble L=0.0;
	IssmDouble eS=0.0;
	IssmDouble Ri=0.0;
	IssmDouble coefM=0.0;
	IssmDouble coefH=0.0;
	IssmDouble coefHT=0.0;
	IssmDouble coefHQ=0.0;
	IssmDouble An_num=0.0;
	IssmDouble An_den_T=0.0;
	IssmDouble An_den_Q=0.0;
	IssmDouble An=0.0;
	IssmDouble C=0.0;
	IssmDouble shf=0.0;
	IssmDouble ds=0.0;
	IssmDouble dAir=0.0;
	IssmDouble TCs=0.0;
	IssmDouble lhf=0.0;
	IssmDouble EC_day=0.0;
	IssmDouble lhf_cum=0.0;
	IssmDouble shf_cum=0.0;
	IssmDouble dT_turb=0.0;
	IssmDouble turb=0.0;
	IssmDouble ulw=0.0;
	IssmDouble dT_ulw=0.0;
	IssmDouble dlw=0.0;
	IssmDouble dT_dlw=0.0;

	/*outputs:*/
	IssmDouble EC=0.0;
	IssmDouble* T=*pT;
	IssmDouble ulwrf=0.0;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   thermal module\n");

	ds = d[0];      // density of top grid cell

	// calculated air density [kg/m3]
	dAir = 0.029 * pAir /(R * Ta);

	// thermal capacity of top grid cell [J/k]
	TCs = d[0]*dz[0]*CI; 

	//initialize Evaporation - Condenstation 
	EC = 0.0;

	// check if all SW applied to surface or distributed throught subsurface
	// swIdx = length(swf) > 1

	// SURFACE ROUGHNESS (Bougamont, 2006)
	// wind/temperature surface roughness height [m]
	if (ds < dIce-Dtol && Ws < Wtol) z0 = 0.00012;        // 0.12 mm for dry snow
	else if (ds >= dIce-Dtol) z0 = 0.0032;             // 3.2 mm for ice
	else z0 = 0.0013;                            // 1.3 mm for wet snow

	//zT and zQ are percentage of z0 (Foken 2008)
	zratio=10.0;
	zT=z0/zratio;
	zQ=z0/zratio;

	// if V = 0 goes to infinity therfore if V = 0 change
	if(V<0.01-Dtol)V=0.01;

	// Bulk-transfer coefficient for turbulent fluxes
	An =  pow(0.4,2); // Bulk-transfer coefficient
	C = An*V;  // shf & lhf common coefficient
	An_den_T = (log(Tz/zT)*log(Vz/z0)); 
	An_den_Q = (log(Tz/zQ)*log(Vz/z0)); 

	// THERMAL CONDUCTIVITY (Sturm, 1997: J. Glaciology)
	// calculate new K profile [W m-1 K-1]

	// initialize conductivity
	K= xNewZeroInit<IssmDouble>(m);

	// for snow and firn (density < 910 kg m-3) (Sturm et al, 1997) or (Calonne et al., 2011)
	if (tcIdx == 2){
		for(int i=0;i<m;i++) if(d[i]<dIce-Dtol) K[i] = 0.024 - 1.23e-4 * d[i] + 2.5e-6 * (pow(d[i],2));
	}
	else{ //default (Sturm et al, 1997)
		for(int i=0;i<m;i++) if(d[i]<dIce-Dtol) K[i] = 0.138 - 1.01e-3 * d[i] + 3.233e-6 * (pow(d[i],2));
	}

	// for ice (density >= 910 kg m-3)
	for(int i=0;i<m;i++) if(d[i]>=dIce-Dtol) K[i] = 9.828 * exp(-5.7e-3*T[i]);

	// THERMAL DIFFUSION COEFFICIENTS

	// A discretization scheme which truncates the Taylor-Series expansion
	// after the 3rd term is used. See Patankar 1980, Ch. 3&4

	// discretized heat equation:

	//                 Tp = (Au*Tuo+ Ad*Tdo+ (Ap-Au-Ad)Tpo+ S) / Ap

	// where neighbor coefficients Au, Ap, & Ad are

	//                   Au = [dz_u/2KU + dz/2KP]^-1
	//                   Ad = [dz_d/2KD + dz/2KP]^-1
	//                   Ap = d*CI*dz/Dt 

	// and u & d represent grid points up and down from the center grid point 
	// point p and o identifies previous time step values. S is a source term.

	// u, d, and p conductivities
	KU = xNew<IssmDouble>(m);
	KD = xNew<IssmDouble>(m);
	KP = xNew<IssmDouble>(m);

	KU[0] = Delflag; //Thermal conductivity of air = 0.025 W/m/K
	KD[m-1] = Delflag;
	for(int i=1;i<m;i++) KU[i]= K[i-1];
	for(int i=0;i<m-1;i++) KD[i] = K[i+1];
	for(int i=0;i<m;i++) KP[i] = K[i];

	// determine u, d & p cell widths
	dzU = xNew<IssmDouble>(m);
	dzD = xNew<IssmDouble>(m);
	dzU[0]=Delflag;
	dzD[m-1]=Delflag;

	for(int i=1;i<m;i++) dzU[i]= dz[i-1];
	for(int i=0;i<m-1;i++) dzD[i] = dz[i+1];

	// determine minimum acceptable delta t (diffusion number > 1/2) [s]
	// NS: 2.16.18 divided dt by scaling factor, default set to 1/11 for stability
	dt=1e12; 
	for(int i=0;i<m;i++) dt = min(dt,CI * pow(dz[i],2) * d[i]  / (3. * K[i]) * thermo_scaling);

	// smallest possible even integer of 60 min where diffusion number > 1/2
	// must go evenly into one hour or the data frequency if it is smaller

	// all integer factors of the number of second in a day (86400 [s])
	int f[45] = {1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 18, 20, 24, 25, 30, 36, 40, 45, 48, 50, 60,
    72, 75, 80, 90, 100, 120, 144, 150, 180, 200, 225, 240, 300, 360, 400, 450, 600, 720, 900, 1200, 1800, 3600};

	// return the min integer factor that is < dt
	max_fdt=f[0];
	bool maxfound=false;
	for(int i=0;i<45;i++){
		if (f[i]<dt-Dtol){
			if (f[i]>=max_fdt-Dtol){
				max_fdt=f[i];
				maxfound=true;
			}
		}
	}
	dt=max_fdt;
	if (maxfound==false){
		if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0){
			_printf0_(" WARNING: calculated timestep for thermal loop is < 1 second.\n");
		}
	}

	// determine mean (harmonic mean) of K/dz for u, d, & p
	Au = xNew<IssmDouble>(m);
	Ad = xNew<IssmDouble>(m);
	Ap = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++){
		Au[i] = pow((dzU[i]/2.0/KU[i] + dz[i]/2.0/KP[i]),-1.0);
		Ad[i] = pow((dzD[i]/2.0/KD[i] + dz[i]/2.0/KP[i]),-1.0);
		Ap[i] = (d[i]*dz[i]*CI)/dt;
	}

	// create "neighbor" coefficient matrix
	Nu = xNewZeroInit<IssmDouble>(m);
	Nd = xNewZeroInit<IssmDouble>(m);
	Np = xNewZeroInit<IssmDouble>(m);
	for(int i=0;i<m;i++){
		Nu[i] = Au[i] / Ap[i];
		Nd[i] = Ad[i] / Ap[i];
		Np[i]= 1.0 - Nu[i] - Nd[i];
	}

	// specify boundary conditions: constant flux at bottom
	Nu[m-1] = 0.0;
	Np[m-1] = 1.0;

	// zero flux at surface
	Np[0] = 1.0 - Nd[0];

	// Create neighbor arrays for diffusion calculations instead of a tridiagonal matrix
	Nu[0] = 0.0;
	Nd[m-1] = 0.0;

	/* RADIATIVE FLUXES*/

	// energy supplied by shortwave radiation [J]
	sw = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) sw[i]= swf[i] * dt;

	// temperature change due to SW
	dT_sw = xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) dT_sw[i]= sw[i] / (CI * d[i] * dz[i]);

	// Upward longwave radiation flux is calculated from the snow surface
	// temperature which is set equal to the average temperature of the
	// top grid cells.

	// energy supplied by downward longwave radiation to the top grid cell [J]
	dlw = dlwrf * dt;

	// temperature change due to dlw_surf
	dT_dlw = dlw / TCs;

	// PREALLOCATE ARRAYS BEFORE LOOP FOR IMPROVED PERFORMANCE
	T0 = xNewZeroInit<IssmDouble>(m+2);
	Tu=xNew<IssmDouble>(m);
	Td=xNew<IssmDouble>(m);

	/* CALCULATE ENERGY SOURCES AND DIFFUSION FOR EVERY TIME STEP [dt]*/
	for (IssmDouble i=1.0;i<=dt0;i+=dt){

		// PART OF ENERGY CONSERVATION CHECK
		// store initial temperature
		//T_init = T;

		// calculate temperature of snow surface (Ts)
		Ts = T[0];
		Ts = min(CtoK,Ts);    // don't allow Ts to exceed 273.15 K (0 degC)

		//TURBULENT HEAT FLUX

		// Monin-Obukhov Stability Correction
		// Reference:
		// Ohmura, A., 1982: Climate and Energy-Balance on the Arctic Tundra.
		// Journal of Climatology, 2, 65-84.

		// calculate the Bulk Richardson Number (Ri)
		Ri = pow(100000./pAir,0.286)*(2.0*9.81*(Ta - Ts)) / (Tz*(Ta + Ts)* pow(V/(Vz),2.0));

		IssmDouble PhiM;
		IssmDouble PhiH;

		// calculate Monin-Obukhov stability factors 'coefM' and 'coefH'
		if (false){
			// do not allow Ri to exceed 0.16
			Ri = min(Ri, 0.16); //Ohmura, 1982

			// calculate momentum 'coefM' stability factor
			if (Ri > 0.0+Ttol){
				// if stable
				coefM = 1.0/(1.0-5.2*Ri);
			}
			else {
				coefM =pow (1.0-18.0*Ri,-0.25);
			}

			// calculate heat/wind 'coef_H' stability factor
			if (Ri <= -0.03+Ttol) coefH = coefM/1.3;
			else coefH = coefM;

         coefHT = coefH*An_den_T;
         coefHQ = coefH*An_den_Q;

		}
		else if(false){

			// do not allow Ri to exceed 0.19
			Ri = min(Ri, 0.19); //Ohmura, 1982

			// calculate momentum 'coefM' stability factor
			if (Ri > 0.0+Ttol){
				// if stable
				//coefM = pow(1.0-5.0*Ri,2.0); //Fitzpatrick et al., 2017, from Brock et al., 2010
				coefM=1.0+5.3*min((Ri/(1.0-5.0*Ri)),0.5);
				coefH=1.0+8.0*min((Ri/(1.0-5.0*Ri)),0.5);
			}
			else {
				//coefM =pow(1.0-16.0*max(Ri,-1.0),0.75); //Fitzpatrick et al., 2017, from Brock et al., 2010
				coefM=pow(1.0-19.0*max(Ri/1.5,-2.0),-0.25);
				coefH=0.95*pow(1.0-11.6*max(Ri/1.5,-2.0),-0.5);
			}

			coefHT = coefH*An_den_T;
			coefHQ = coefH*An_den_Q;

		}
      else if (false){
         // Greuell and Konzelman, 1994
         // calculate momentum 'coefM' stability factor

         if (Ri > 0.0+Ttol){
            // if stable
            coefM=1.0+15.0*Ri*pow(1.0+Ri,1./2.);
            coefH=1.0;
         }
         else {
            coefM=pow(1.0-15.0*Ri/(1.0+75.0*pow(0.4/log(Tz/zT),2)*pow(Tz/zT*fabs(Ri),1./2.)),-1);
            coefH=1.0;
         }

         coefHT = coefH*An_den_T;
         coefHQ = coefH*An_den_Q;

      }
      else {
         IssmDouble a1=1.0;
         IssmDouble b1=2.0/3.0;
         IssmDouble c1=5.0;
         IssmDouble d1=0.35;
         IssmDouble PhiMz=0.0;
         IssmDouble PhiHz=0.0;
         IssmDouble PhiMz0=0.0;
         IssmDouble PhiHzT=0.0;
         IssmDouble PhiHzQ=0.0;
         IssmDouble zL=0.0;
         IssmDouble zLT=0.0;
         IssmDouble zLM=0.0;

         if (Ri > 0.0+Ttol){
            // if stable

            if(Ri < 0.2-Ttol){
               zL = Ri/(1.0-5.0*Ri);
            }
            else{
               zL=Ri;
            }
            //zL = min(zL, 0.5); //Sjoblom, 2014
            zLM=max(zL/Vz*z0,1e-3);
            zLT=max(zL/Tz*zT,1e-3);

            //Ding et al. 2020, from Beljaars and Holtslag (1991)
            PhiMz=-1.*(a1*zL + b1*(zL-c1/d1)*exp(-1.*d1*zL) + b1*c1/d1);
            PhiHz=-1.*(pow(1.+2.*a1*zL/3.,1.5) + b1*(zL-c1/d1)*exp(-1.*d1*zL) + b1*c1/d1 - 1.0);
            PhiMz0=-1.*(a1*zLM + b1*(zLM-c1/d1)*exp(-1.*d1*zLM) + b1*c1/d1);
            PhiHzT=-1.*(pow(1.+2.*a1*zLT/3.,1.5) + b1*(zLT-c1/d1)*exp(-1.*d1*zLT) + b1*c1/d1 - 1.0);

            PhiHzQ=PhiHzT;
         }
         else {
            IssmDouble xm;
            IssmDouble xh;
            IssmDouble xmT;
            IssmDouble xmM;

            zL = Ri/1.5; //max(Ri, -0.5+Ttol)/1.5; //Hogstrom (1996)
            //zL = max(zL, -2.0); //Sjoblom, 2014
            zLM=min(zL/Vz*z0,-1e-3);
            zLT=min(zL/Tz*zT,-1e-3);

            if (true){ //Sjoblom, 2014
               xm=pow(1.0-19.0*zL,-0.25);
               PhiMz=2.0*log((1.+xm)/2.0) + log((1.+pow(xm,2))/2.0) - 2.*atan(xm) + Pi/2.;

               xh=0.95*pow(1.0-11.6*zL,-0.5);
               PhiHz=2.0*log((1.0+pow(xh,2))/2.0);
            }
            else{ //Ding et al., 2020
               xm=pow(1.0-16*zL,0.25);
               xmM=pow(1.0-16*zLM,0.25);
               xmT=pow(1.0-16*zLT,0.25);
               PhiMz=2.0*log((1.+xm)/2.0) + log((1.+pow(xm,2))/2.0) - 2.0*atan(xm) + Pi/2.0;
               PhiMz0=2.0*log((1.+xmM)/2.0) + log((1.+pow(xmM,2))/2.0) - 2.0*atan(xmM) + Pi/2.0;

               PhiHz=2.0*log((1.+pow(xm,2))/2.0);
               PhiHzT=2.0*log((1.+pow(xmT,2))/2.0);

               PhiHzQ=PhiHzT;
            }
         }

         PhiM=PhiMz;
         PhiH=PhiHz;
         coefM = log(Vz/z0) - PhiMz + PhiMz0; //Ding et al., 2019
         coefHT = log(Tz/zT) - PhiHz + PhiHzT; //Sjoblom, 2014, after Foken 2008
         coefHQ = log(Tz/zQ) - PhiHz + PhiHzQ; //Sjoblom, 2014, after Foken 2008

      }

      //// Sensible Heat
      // calculate the sensible heat flux [W m-2](Patterson, 1998)
      shf = dAir * C * CA * (Ta - Ts) * pow(100000./pAir,0.286);

      // adjust using Monin-Obukhov stability theory
      shf = shf/(coefM*coefHT);

      //// Latent Heat
      // determine if snow pack is melting & calcualte surface vapour pressure over ice or liquid water
      if (Ts >= CtoK-Ttol){
         L = LV; //for liquid water at 273.15 k to vapor

         //for liquid surface (assume liquid on surface when Ts == 0 deg C)
         // Wright (1997), US Meteorological Handbook from Murphy and Koop, 2005 Appendix A
         //eS = 611.21 * exp(17.502 * (Ts - CtoK) / (240.97 + Ts - CtoK));
         // Murray 1967, https://cran.r-project.org/web/packages/humidity/vignettes/humidity-measures.html
         eS = 610.78 * exp(17.2693882 * (Ts - CtoK - 0.01) / (Ts - 35.86));
      }
      else{
         L = LS; // latent heat of sublimation

         // for an ice surface Murphy and Koop, 2005 [Equation 7]
         //eS = exp(9.550426 - 5723.265/Ts + 3.53068 * log(Ts) - 0.00728332 * Ts);
         // for an ice surface Ding et al., 2019 after Bolton, 1980
         eS = 610.78 * exp(21.8745584 * (Ts - CtoK - 0.01) / (Ts - 7.66));
      }

      // Latent heat flux [W m-2]
      lhf = C * L * (eAir - eS) / (461.9*(Ta+Ts)/2.0);

      // adjust using Monin-Obukhov stability theory (if lhf '+' then there is energy and mass gained at the surface,
      // if '-' then there is mass and energy loss at the surface.
      lhf = lhf/(coefM*coefHQ);

		//mass loss (-)/acreation(+) due to evaporation/condensation [kg]
		EC_day = lhf * dts / L;

		// temperature change due turbulent fluxes
		turb = (shf + lhf)* dt;
		dT_turb = turb  / TCs;

		// upward longwave contribution
		IssmDouble deltaULW=0.0;
		IssmDouble emissivity=1.0;
		//If user wants to set a upward long wave bias
		if(isdeltaLWup) deltaULW = dulwrfValue;
		//If user wants to directly set emissivity, or grain radius is larger than the
		// threshold, or eIdx is 2 and we have wet snow or ice, use prescribed emissivity
		if(eIdx==0 || (teThresh - re[0])<Gdntol || (eIdx==2 && z0>(0.001+Gdntol))) emissivity = teValue;
		ulw = - (SB * pow(Ts,4.0)* emissivity + deltaULW) * dt; 
		ulwrf = ulwrf - ulw/dt0;

		dT_ulw = ulw / TCs;

		// new grid point temperature

		//SW penetrates surface
		if(!isconstrainsurfaceT){
			for(int j=0;j<m;j++) T[j] = T[j] + dT_sw[j];
			T[0] = T[0] + dT_dlw + dT_ulw + dT_turb;
		}

		// temperature diffusion
		for(int j=0;j<m;j++) T0[1+j]=T[j];
		T0[0]=Ta;
		T0[m+1]=T[m-1];
		for(int j=0;j<m;j++) Tu[j] = T0[j];
		for(int j=0;j<m;j++) Td[j] = T0[2+j];
		for(int j=0;j<m;j++) T[j] = (Np[j] * T[j]) + (Nu[j] * Tu[j]) + (Nd[j] * Td[j]);

		// calculate cumulative evaporation (+)/condensation(-)
		if(!isconstrainsurfaceT) EC = EC + (EC_day/dts)*dt;
		lhf_cum=lhf_cum+lhf*dt/dt0;
		shf_cum=shf_cum+shf*dt/dt0;

		/* CHECK FOR ENERGY (E) CONSERVATION [UNITS: J]
		//energy flux across lower boundary (energy supplied by underling ice)
		base_flux = Ad(-1)*(T_init()-T_init(-1)) * dt;

		E_used = sum((T - T_init) * (d*dz*CI));
		E_sup = ((sum(swf)  * dt) + dlw + ulw + turb + base_flux);

		E_diff = E_used - E_sup;

		if fabs(E_diff) > 1E-6 || isnan(E_diff)
		disp(T(1))
		_error_("energy not conserved in thermodynamics equations");
		*/
	}

	/*Free resources:*/
	xDelete<IssmDouble>(K);
	xDelete<IssmDouble>(KU);
	xDelete<IssmDouble>(KD);
	xDelete<IssmDouble>(KP);
	xDelete<IssmDouble>(Au);
	xDelete<IssmDouble>(Ad);
	xDelete<IssmDouble>(Ap);
	xDelete<IssmDouble>(Nu);
	xDelete<IssmDouble>(Nd);
	xDelete<IssmDouble>(Np);
	xDelete<IssmDouble>(dzU);
	xDelete<IssmDouble>(dzD);
	xDelete<IssmDouble>(sw);
	xDelete<IssmDouble>(dT_sw);
	xDelete<IssmDouble>(T0);
	xDelete<IssmDouble>(Tu);
	xDelete<IssmDouble>(Td);

	/*Assign output pointers:*/
	*pEC=EC;
	*plhf=lhf_cum;
	*pshf=shf_cum;
	*pT=T;
	*pulwrf=ulwrf;

}  /*}}}*/
void shortwave(IssmDouble** pswf, int swIdx, int aIdx, IssmDouble dsw, IssmDouble dswdiff, IssmDouble as, IssmDouble asdiff, IssmDouble* d, IssmDouble* dz, IssmDouble* re, IssmDouble dIce, int m, int sid){ /*{{{*/

	// DISTRIBUTES ABSORBED SHORTWAVE RADIATION WITHIN SNOW/ICE

	// swIdx = 0 : all absorbed SW energy is assigned to the top grid cell

	// swIdx = 1 : absorbed SW is distributed with depth as a function of:
	//   default   : snow density (taken from Bassford, 2002)
	//   if aIdx=2 : grain size in 3 spectral bands (Brun et al., 1992)

	// Inputs
	//   swIdx   = shortwave allowed to penetrate surface (0 = No, 1 = Yes)
	//   aIdx    = method for calculating albedo (1-4)
	//   dsw     = downward shortwave radiative flux [w m-2]
	//   dswdiff  = downward shortwave diffuse radiative flux [w m-2]
	//   as      = surface albedo
	//   asdiff  = surface albedo for diffuse radiation
	//   d       = grid cell density [kg m-3]
	//   dz      = grid cell depth [m]
	//   re      = grid cell effective grain radius [mm]

	// Outputs
	//   swf     = absorbed shortwave radiation [W m-2]
	//

	/*outputs: */
	IssmDouble* swf=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   shortwave module\n");

	/*Initialize and allocate: */
	swf=xNewZeroInit<IssmDouble>(m);

	// SHORTWAVE FUNCTION
	if (swIdx == 0 || (dIce - d[0])<Dtol) {// all sw radation is absorbed in by the top grid cell

		// calculate surface shortwave radiation fluxes [W m-2]
		if (aIdx == 1){
			swf[0] = (1.0 - as) * max(0.0,(dsw - dswdiff)) +  (1.0 - asdiff) * dswdiff;
		}
		else{
			swf[0] = (1.0 - as) * dsw;
		}
	}
	else{ // sw radation is absorbed at depth within the glacier

		if (aIdx == 2){    // function of effective radius (3 spectral bands)

			IssmDouble * gsz=NULL;
			IssmDouble * B1_cum=NULL;
			IssmDouble * B2_cum=NULL;
			IssmDouble* h =NULL;
			IssmDouble* B1 =NULL;
			IssmDouble* B2 =NULL;
			IssmDouble* exp1 = NULL;
			IssmDouble* exp2 = NULL;
			IssmDouble*  Qs1 = NULL;
			IssmDouble*  Qs2 = NULL;

			// convert effective radius [mm] to grain size [m]
			gsz=xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) gsz[i]= (re[i] * 2.0) / 1000.0;

			// Spectral fractions [0.3-0.8um 0.8-1.5um 1.5-2.8um]
			// (Lefebre et al., 2003)
			IssmDouble sF[3] = {0.606, 0.301, 0.093};

			// initialize variables
			B1_cum=xNew<IssmDouble>(m+1);
			B2_cum=xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++){
				B1_cum[i]=1.0;
				B2_cum[i]=1.0;
			}

			// spectral albedos:
			// 0.3 - 0.8um
			IssmDouble a0 = min(0.98, 0.95 - 1.58 *pow(gsz[0],0.5));
			// 0.8 - 1.5um
			IssmDouble a1 = max(0.0, 0.95 - 15.4 *pow(gsz[0],0.5));
			// 1.5 - 2.8um
			IssmDouble a2 = max(0.127, 0.88 + 346.3*gsz[0] - 32.31*pow(gsz[0],0.5));

			// separate net shortwave radiative flux into spectral ranges
			IssmDouble swfS[3];
			swfS[0] = (sF[0] * dsw) * (1.0 - a0);
			swfS[1] = (sF[1] * dsw) * (1.0 - a1);
			swfS[2] = (sF[2] * dsw) * (1.0 - a2);

			// absorption coefficient for spectral range:
			h =xNew<IssmDouble>(m);
			B1 =xNew<IssmDouble>(m);
			B2 =xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) h[i]= d[i] /(pow(gsz[i],0.5));
			for(int i=0;i<m;i++) B1[i] = 0.0192 * h[i];                 // 0.3 - 0.8um
			for(int i=0;i<m;i++) B2[i]= 0.1098 * h[i];                 // 0.8 - 1.5um
			// B3 = +inf                     // 1.5 - 2.8um

			// cumulative extinction factors
			exp1 = xNew<IssmDouble>(m); 
			exp2 = xNew<IssmDouble>(m); 
			for(int i=0;i<m;i++) exp1[i]=exp(-B1[i]*dz[i]);
			for(int i=0;i<m;i++) exp2[i]=exp(-B2[i]*dz[i]);

			for(int i=0;i<m;i++){
				IssmDouble cum1=exp1[0];
				IssmDouble cum2=exp2[0];
				for(int j=1;j<=i;j++){
					cum1 = cum1*exp1[j];
					cum2 = cum2*exp2[j];
				}
				B1_cum[i+1]=cum1;
				B2_cum[i+1]=cum2;
			}

			// flux across grid cell boundaries
			Qs1 = xNew<IssmDouble>(m+1);
			Qs2 = xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++){
				Qs1[i] = swfS[0] * B1_cum[i];
				Qs2[i] = swfS[1] * B2_cum[i];
			}

			// net energy flux to each grid cell
			for(int i=0;i<m;i++) swf[i]= (Qs1[i]-Qs1[i+1]) + (Qs2[i]-Qs2[i+1]);

			// add flux absorbed at surface
			swf[0] = swf[0]+ swfS[2];

			/*Free resources: */
			xDelete<IssmDouble>(gsz);
			xDelete<IssmDouble>(B1_cum);
			xDelete<IssmDouble>(B2_cum);
			xDelete<IssmDouble>(h);
			xDelete<IssmDouble>(B1);
			xDelete<IssmDouble>(B2);
			xDelete<IssmDouble>(exp1);
			xDelete<IssmDouble>(exp2);
			xDelete<IssmDouble>(Qs1);
			xDelete<IssmDouble>(Qs2);

		}
		else{  //function of grid cell density

			/*intermediary: */
			IssmDouble* B_cum = NULL;
			IssmDouble* exp_B = NULL;
			IssmDouble* Qs    = NULL;
			IssmDouble* B    = NULL;

			// fraction of sw radiation absorbed in top grid cell (wavelength > 0.8um)
			IssmDouble SWs = 0.36;

			// SWs and SWss coefficients need to be better constranted. Greuell
			// and Konzelmann 1994 used SWs = 0.36 and SWss = 0.64 as this the
			// the % of SW radiation with wavelengths > and < 800 nm
			// respectively.  This, however, may not account for the fact that
			// the albedo of wavelengths > 800 nm has a much lower albedo.

			// calculate surface shortwave radiation fluxes [W m-2]
			IssmDouble swf_s = SWs * (1.0 - as) * dsw;

			// calculate surface shortwave radiation fluxes [W m-2]
			IssmDouble swf_ss = (1.0-SWs) * (1.0 - as) * dsw;

			// SW allowed to penetrate into snowpack
			IssmDouble Bs = 10.0;    // snow SW extinction coefficient [m-1] (Bassford,2006)
			IssmDouble Bi = 1.3;   // ice SW extinction coefficient [m-1] (Bassford,2006)

			// calculate extinction coefficient B [m-1] vector
			B=xNew<IssmDouble>(m);
			for(int i=0;i<m;i++) B[i] = Bs + (300.0 - d[i]) * ((Bs - Bi)/(dIce - 300.0));

			// cumulative extinction factor
			B_cum = xNew<IssmDouble>(m+1);
			exp_B = xNew<IssmDouble>(m);
			for(int i=0;i<m;i++)exp_B[i]=exp(-B[i]*dz[i]);

			B_cum[0]=1.0;
			for(int i=0;i<m;i++){
				IssmDouble cum_B=exp_B[0];
				for(int j=1;j<=i;j++) cum_B=cum_B*exp_B[j];
				B_cum[i+1]=  cum_B;
			}

			// flux across grid cell boundaries
			Qs=xNew<IssmDouble>(m+1);
			for(int i=0;i<m+1;i++) Qs[i] = swf_ss * B_cum[i];

			// net energy flux to each grid cell
			for(int i=0;i<m;i++) swf[i] = (Qs[i]-Qs[i+1]);

			// add flux absorbed at surface
			swf[0] += swf_s;

			/*Free resources:*/
			xDelete<IssmDouble>(B_cum);
			xDelete<IssmDouble>(exp_B);
			xDelete<IssmDouble>(Qs);
			xDelete<IssmDouble>(B);
		}
	}
	/*Assign output pointers: */
	*pswf=swf;

} /*}}}*/ 
void accumulation(IssmDouble** pT, IssmDouble** pdz, IssmDouble** pd, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, IssmDouble* pRa, int* pm, int aIdx, int dsnowIdx, IssmDouble Tmean, IssmDouble T_air, IssmDouble P, IssmDouble dzMin, IssmDouble aSnow, IssmDouble C, IssmDouble V, IssmDouble Vmean, IssmDouble dIce, int sid){ /*{{{*/

	// Adds precipitation and deposition to the model grid

	// Author: Alex Gardner, University of Alberta
	// Date last modified: JAN, 2008

	/* Description:
	   adjusts the properties of the top grid cell to account for accumulation
	   T_air & T = Air and top grid cell temperatures [K]
		Tmean =  average surface temperature [K]
		Vmean =  average wind velocity [m s-1]
		V =  wind velocity [m s-1]
		C =  average accumulation rate [kg m-2 yr-1]
	   dz = topgrid cell length [m]
	   d = density of top grid gell [kg m-3]
		P = precipitation [mm w.e.] or [kg m-3]
		Ra = rainfall [mm w.e.] or [kg m-3]
		re = effective grain radius [mm]
		gdn = grain dentricity
		gsp = grain sphericity*/

	// MAIN FUNCTION
	// specify constants
	IssmDouble dSnow = 150.0;    // density of snow [kg m-3]
	IssmDouble reNew = 0.05;    // new snow grain size [mm]
	IssmDouble gdnNew = 1.0;     // new snow dendricity 
	IssmDouble gspNew = 0.5;   // new snow sphericity 

	/*intermediary: */
	IssmDouble* mInit=NULL;
	bool        top=true;
	IssmDouble  mass=0.0;
	IssmDouble  massinit=0.0;
	IssmDouble  mass_diff=0.0;

	/*output: */
	IssmDouble* T=NULL;
	IssmDouble* dz=NULL;
	IssmDouble* d=NULL;
	IssmDouble* W=NULL;
	IssmDouble* a=NULL;
	IssmDouble* adiff=NULL;
	IssmDouble* re=NULL;
	IssmDouble* gdn=NULL;
	IssmDouble* gsp=NULL;
	IssmDouble  ra=0.0;
	int         m=0;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   accumulation module\n");

	/*Recover pointers: */
	T=*pT;
	dz=*pdz;
	d=*pd;
	W=*pW;
	a=*pa;
	adiff=*padiff;
	re=*pre;
	gdn=*pgdn;
	gsp=*pgsp;
	m=*pm;
	ra=*pRa;

	//Density of fresh snow [kg m-3]
	switch (dsnowIdx){
		case 0: // Default value defined above
			break;

		case 1: // Density of Antarctica snow
			dSnow = 350.0;
			//dSnow = 360.0; //FirnMICE Lundin et al., 2017
			break;

		case 2: // Density of Greenland snow, Fausto et al., 2018
			dSnow = 315.0;
			//From Vionnet et al., 2012 (Crocus)
			gdnNew = min(max(1.29 - 0.17*V,0.20),1.0);
			gspNew = min(max(0.08*V + 0.38,0.5),0.9);
			reNew=max(1e-1*(gdnNew/.99+(1.0-1.0*gdnNew/.99)*(gspNew/.99*3.0+(1.0-gspNew/.99)*4.0))/2.0,Gdntol);
			break;

		case 3: //Surface snow accumulation density from Kaspers et al., 2004, Antarctica
			//dSnow = alpha1 + beta1*T + delta1*C + epsilon1*W
			//     7.36x10-2  1.06x10-3  6.69x10-2  4.77x10-3 
			dSnow=(7.36e-2 + 1.06e-3*min(Tmean,CtoK-Ttol) + 6.69e-2*C/1000. + 4.77e-3*Vmean)*1000.;
			break;

		case 4: // Kuipers Munneke and others (2015), Greenland
			dSnow = 481.0 + 4.834*(Tmean-CtoK);
			break;
	}

	// determine initial mass
	mInit=xNew<IssmDouble>(m);
	for(int i=0;i<m;i++) mInit[i]= d[i] * dz[i];
	massinit=0.0; 
	for(int i=0;i<m;i++)massinit+=mInit[i];

	if (P > 0.0+Ptol){

		if (T_air <= CtoK+Ttol){ // if snow

			IssmDouble z_snow = P/dSnow;               // depth of snow
			IssmDouble dfall = gdnNew;
			IssmDouble sfall = gspNew;
			IssmDouble refall = reNew;

			// if snow depth is greater than specified min dz, new cell created
			if (z_snow > dzMin+Dtol){

				newcell(&T,T_air,top,m); //new cell T
				newcell(&dz,z_snow,top,m); //new cell dz
				newcell(&d,dSnow,top,m); //new cell d
				newcell(&W,0.0,top,m); //new cell W
				newcell(&a,aSnow,top,m); //new cell a
				newcell(&adiff,aSnow,top,m); //new cell a
				newcell(&re,refall,top,m); //new cell grain size
				newcell(&gdn,dfall,top,m); //new cell grain dendricity
				newcell(&gsp,sfall,top,m); //new cell grain sphericity
				m=m+1;
			}
			else { // if snow depth is less than specified minimum dz snow

				mass = mInit[0] + P;         // grid cell adjust mass

				dz[0] = dz[0] + P/dSnow;    // adjust grid cell depth      
				d[0] = mass / dz[0];    // adjust grid cell density

				// adjust variables as a linearly weighted function of mass
				// adjust temperature (assume P is same temp as air)
				T[0] = (T_air * P + T[0] * mInit[0])/mass;

				// adjust a, re, gdn & gsp
				if(aIdx>0)a[0] = (aSnow * P + a[0] * mInit[0])/mass;
				gdn[0] = dfall; 
				gsp[0] = sfall; 
				re[0] = max(1e-1*(gdn[0]/.99+(1.0-1.0*gdn[0]/.99)*(gsp[0]/.99*3.0+(1.0-gsp[0]/.99)*4.0))/2.0,Gdntol);
			}
		}
		else{ // if rain    

			/*rain is added by increasing the mass and temperature of the ice
			  of the top grid cell.  Temperatures are set artifically high to
			  account for the latent heat of fusion.  This is the same as
			  directly adding liquid water to the the snow pack surface but
			  makes the numerics easier.*/

			// grid cell adjust mass
			mass = mInit[0] + P;

			// adjust temperature
			// liquid: must account for latent heat of fusion
			T[0] = (P *(T_air + LF/CI) + T[0] * mInit[0]) / mass;

			// adjust grid cell density
			d[0] = mass / dz[0];

			// if d > the density of ice, d = dIce
			if (d[0] > dIce-Dtol){
				d[0] = dIce;           // adjust d
				dz[0] = mass / d[0];    // dz is adjusted to conserve mass
			}

			ra=P;
		}

		// check for conservation of mass
		mass=0; for(int i=0;i<m;i++)mass+=d[i]*dz[i]; 

		mass_diff = mass - massinit - P;

		#ifndef _HAVE_AD_  //avoid round operation. only check in forward mode.
		mass_diff = round(mass_diff * 100.0)/100.0;
		if (mass_diff > 0) _error_("mass not conserved in accumulation function");
		#endif

	}
	/*Free resources:*/
	if(mInit)xDelete<IssmDouble>(mInit);

	/*Assign output pointers:*/
	*pT=T;
	*pdz=dz;
	*pd=d;
	*pW=W;
	*pa=a;
	*padiff=adiff;
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;
	*pm=m;
	*pRa=ra;
} /*}}}*/
void melt(IssmDouble* pM, IssmDouble* pMs, IssmDouble* pR, IssmDouble* pF, IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble Ra, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, IssmDouble zY, IssmDouble dIce, int sid){ /*{{{*/

	//// MELT ROUTINE

	// Description:
	// computes the quantity of meltwater due to snow temperature in excess of
	// 0 deg C, determines pore water content and adjusts grid spacing

	/*intermediary:*/
	IssmDouble* m=NULL;
	IssmDouble* maxF=NULL;
	IssmDouble* dW=NULL;
	IssmDouble* exsW=NULL;
	IssmDouble* exsT=NULL;
	IssmDouble* surpT=NULL;
	IssmDouble* surpE=NULL;
	IssmDouble* flxDn=NULL;
	IssmDouble  ER=0.0;
	IssmDouble* EI=NULL;
	IssmDouble* EW=NULL;
	IssmDouble* M=NULL;
	int*        D=NULL;

	IssmDouble sumM=0.0;
	IssmDouble sumER=0.0;
	IssmDouble addE=0.0;
	IssmDouble mSum0=0.0;
	IssmDouble sumE0=0.0;
	IssmDouble mSum1=0.0;
	IssmDouble sumE1=0.0;
	IssmDouble dE=0.0;
	IssmDouble dm=0.0;
	IssmDouble X=0.0;
	IssmDouble Wi=0.0;

	int        D_size = 0;

	/*outputs:*/
	IssmDouble  Msurf = 0.0;
	IssmDouble  mAdd = 0.0;
	IssmDouble  surplusE = 0.0;
	IssmDouble  surplusT = 0.0;
	IssmDouble  dz_add = 0.0;
	IssmDouble  Rsum = 0.0;
	IssmDouble  Fsum = 0.0;
	IssmDouble* T=*pT;
	IssmDouble* d=*pd;
	IssmDouble* dz=*pdz;
	IssmDouble* W=*pW;
	IssmDouble* a=*pa;
	IssmDouble* adiff=*padiff;
	IssmDouble* re=*pre;
	IssmDouble* gdn=*pgdn;
	IssmDouble* gsp=*pgsp;
	int         n=*pn;
	IssmDouble* R=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   melt module\n");

	//// INITIALIZATION

	/*Allocations: */
	M=xNewZeroInit<IssmDouble>(n); 
	maxF=xNewZeroInit<IssmDouble>(n); 
	dW=xNewZeroInit<IssmDouble>(n); 

	// store initial mass [kg] and energy [J]
	m=xNew<IssmDouble>(n); for(int i=0;i<n;i++) m[i] = dz[i]* d[i];                    // grid cell mass [kg]
	EI=xNew<IssmDouble>(n); for(int i=0;i<n;i++)EI[i] = m[i] * T[i] * CI;               // initial enegy of snow/ice
	EW=xNew<IssmDouble>(n); for(int i=0;i<n;i++)EW[i]= W[i] * (LF + CtoK * CI);     // initial enegy of water

	mSum0 = cellsum(W,n) + cellsum(m,n);        // total mass [kg]
	sumE0 = cellsum(EI,n) + cellsum(EW,n);      // total energy [J]

	// initialize melt and runoff scalars
	Rsum = 0.0;       // runoff [kg]
	Fsum = 0.0;       // refreeze [kg]
	sumM = 0.0;       // total melt [kg]
	mAdd = 0.0;       // mass added/removed to/from base of model [kg]
	addE = 0.0;       // energy added/removed to/from base of model [J]
	dz_add=0.0;      // thickness of the layer added/removed to/from base of model [m]

	// calculate temperature excess above 0 deg C
	exsT=xNewZeroInit<IssmDouble>(n);
	for(int i=0;i<n;i++) exsT[i]= max(0.0, T[i] - CtoK);        // [K] to [degC]

	// new grid point center temperature, T [K]
	// for(int i=0;i<n;i++) T[i]-=exsT[i];
	for(int i=0;i<n;i++) T[i]=min(T[i],CtoK);

	// specify irreducible water content saturation [fraction]
	const IssmDouble Swi = 0.07;                     // assumed constant after Colbeck, 1974
	const IssmDouble dPHC = 830.0;                     //Pore closeoff density

	//// REFREEZE PORE WATER
	// check if any pore water
	if (cellsum(W,n) > 0.0+Wtol){
		if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("      pore water refreeze\n");
		// calculate maximum freeze amount, maxF [kg]
		for(int i=0;i<n;i++) maxF[i] = max(0.0, -((T[i] - CtoK) * m[i] * CI) / LF);

		// freeze pore water and change snow/ice properties
		for(int i=0;i<n;i++) dW[i] = min(maxF[i], W[i]);    // freeze mass [kg]   
		for(int i=0;i<n;i++) W[i] -= dW[i];                                            // pore water mass [kg]
		for(int i=0;i<n;i++) m[i] += dW[i];                                            // new mass [kg]
		for(int i=0;i<n;i++) d[i] = m[i] / dz[i];                                    // density [kg m-3]   
		for(int i=0;i<n;i++) if(m[i]>Wtol) T[i] = T[i] + (dW[i]*(LF+(CtoK - T[i])*CI)/(m[i]*CI));      // temperature [K]

		// if pore water froze in ice then adjust d and dz thickness
		for(int i=0;i<n;i++)if(d[i]> dIce-Dtol)d[i]=dIce;
		for(int i=0;i<n;i++) dz[i]= m[i]/d[i];

	}

	// squeeze water from snow pack
	exsW=xNew<IssmDouble>(n); 
	for(int i=0;i<n;i++){
		Wi= (dIce - d[i]) * Swi * (m[i] / d[i]);        // irreducible water content [kg]
		exsW[i] = max(0.0, W[i] - Wi);                  // water "squeezed" from snow [kg]
	}

	//// MELT, PERCOLATION AND REFREEZE

	// initialize refreeze, runoff, flxDn and dW vectors [kg]
	IssmDouble* F = xNewZeroInit<IssmDouble>(n);

	// Add previous refreeze to F and reset dW
	for(int i=0;i<n;i++){
		F[i]=F[i]+dW[i];
		dW[i] = 0.0;
	}

	// run melt algorithm if there is melt water or excess pore water
	if ((cellsum(exsT,n) > 0.0+Ttol) || (cellsum(exsW,n) > 0.0+Wtol)){
		// _printf_(""MELT OCCURS");
		// check to see if thermal energy exceeds energy to melt entire cell
		// if so redistribute temperature to lower cells (temperature surplus)
		// (maximum T of snow before entire grid cell melts is a constant
		// LF/CI = 159.1342)
		surpT=xNew<IssmDouble>(n); for(int i=0;i<n;i++)surpT[i] = max(0.0, exsT[i]- LF/CI);

		if (cellsum(surpT,n) > 0.0+Ttol ){
			// _printf_("T Surplus");
			// calculate surplus energy
			surpE=xNew<IssmDouble>(n); for(int i=0;i<n;i++)surpE[i] = surpT[i] * CI * m[i];

			int i = 0;
			while (cellsum(surpE,n) > 0.0+Ttol && i<n){

				if (i<n-1){
					// use surplus energy to increase the temperature of lower cell
					T[i+1] = surpE[i]/m[i+1]/CI + T[i+1];

					exsT[i+1] = max(0.0, T[i+1] - CtoK) + exsT[i+1];
					T[i+1] = min(CtoK, T[i+1]);

					surpT[i+1] = max(0.0, exsT[i+1] - LF/CI);
					surpE[i+1] = surpT[i+1] * CI * m[i+1];
				}
				else{
					surplusT=max(0.0, exsT[i] - LF/CI);
					surplusE=surpE[i];
					if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0){
						_printf0_(" WARNING: surplus energy at the base of GEMB column\n");
					}
				}

				// adjust current cell properties (again 159.1342 is the max T)
				exsT[i] = LF/CI;
				surpE[i] = 0.0;
				i = i + 1;
			}
		}

		// convert temperature excess to melt [kg]
		IssmDouble Mmax=0.0;
		for(int i=0;i<n;i++){
			Mmax=exsT[i] * d[i] * dz[i] * CI / LF;
			M[i] = min(Mmax, m[i]);  // melt
		}
		Msurf = M[0];
		sumM = max(0.0,cellsum(M,n)-Ra);  // total melt [kg] minus the liquid rain that had been added

		// calculate maximum refreeze amount, maxF [kg]
		for(int i=0;i<n;i++)maxF[i] = max(0.0, -((T[i] - CtoK) * d[i] * dz[i] * CI)/ LF);

		// initialize refreeze, runoff, flxDn and dW vectors [kg]
		IssmDouble* R = xNewZeroInit<IssmDouble>(n);

		flxDn=xNewZeroInit<IssmDouble>(n+1);

		// determine the deepest grid cell where melt/pore water is generated
		X = 0;
		for(int i=n-1;i>=0;i--){
			if(M[i]> 0.0+Wtol || exsW[i]> 0.0+Wtol){
				X=i;
				break;
			}
		}

		IssmDouble depthice=0.0;
		int Xi=0;
		//// meltwater percolation
		for(int i=0;i<n;i++){
			// calculate total melt water entering cell
			IssmDouble inM = M[i]+ flxDn[i];

			depthice=0.0;
			if (d[i] >= dPHC-Dtol){
				for(int l=i;(l<n && d[l]>=dPHC-Dtol);l++) depthice=depthice+dz[l];
			}

			// break loop if there is no meltwater and if depth is > mw_depth
			if (fabs(inM) < Wtol && i > X){
				break;
			}
			// if reaches impermeable ice layer all liquid water runs off (R)
			else if (d[i] >= dIce-Dtol || (d[i] >= dPHC-Dtol && depthice>0.1+Dtol)){  // dPHC = pore hole close off [kg m-3]
				// _printf_("ICE LAYER");
				// no water freezes in this cell
				// no water percolates to lower cell
				// cell ice temperature & density do not change

				m[i] = m[i] - M[i];                     // mass after melt
				Wi = (dIce-d[i]) * Swi * (m[i]/d[i]);    // irreducible water 
				dW[i] = max(min(inM, Wi - W[i]),-1*W[i]);            // change in pore water
				R[i] = max(0.0, inM - dW[i]);             // runoff
			}
			// check if no energy to refreeze meltwater
			else if (fabs(maxF[i]) < Dtol){
				// _printf_("REFREEZE == 0");
				// no water freezes in this cell
				// cell ice temperature & density do not change

				m[i] = m[i] - M[i];                     // mass after melt
				Wi = (dIce-d[i]) * Swi * (m[i]/d[i]);    // irreducible water 
				dW[i] = max(min(inM, Wi - W[i]),-1*W[i]);              // change in pore water
				flxDn[i+1] = max(0.0, inM - dW[i]);         // meltwater out
				R[i] = 0.0;
			}
			// some or all meltwater refreezes
			else{
				// change in density and temperature
				// _printf_("MELT REFREEZE");
				//-----------------------melt water-----------------------------
				m[i] = m[i] - M[i];
				IssmDouble dz_0 = m[i]/d[i];          
				IssmDouble dMax = (dIce - d[i])*dz_0;              // d max = dIce
				IssmDouble F1 = min(min(inM,dMax),maxF[i]);         // maximum refreeze               
				m[i] = m[i] + F1;                       // mass after refreeze
				d[i] = m[i]/dz_0;

				//-----------------------pore water-----------------------------
				Wi = (dIce-d[i])* Swi * dz_0;            // irreducible water 
				dW[i] = max(min(inM - F1, Wi-W[i]),-1*W[i]);         // change in pore water
				IssmDouble F2 = 0.0;                                 

				if (dW[i] < 0.0-Wtol){                         // excess pore water
					dMax = (dIce - d[i])*dz_0;          // maximum refreeze                                             
					IssmDouble maxF2 = min(dMax, maxF[i]-F1);      // maximum refreeze
					F2 = min(-1.0*dW[i], maxF2);            // pore water refreeze
					m[i] = m[i] + F2;                   // mass after refreeze
					d[i] = m[i]/dz_0;
				}

				F[i] = F[i] + F1 + F2;

				flxDn[i+1] = max(0.0,inM - F1 - dW[i]);     // meltwater out        
				if (m[i]>Wtol){
					T[i] = T[i] + ((F1+F2)*(LF+(CtoK - T[i])*CI)/(m[i]*CI));// change in temperature
				}

				// check if an ice layer forms 
				if (fabs(d[i] - dIce) < Dtol){
					// _printf_("ICE LAYER FORMS");
					// excess water runs off
					R[i] = flxDn[i+1];
					// no water percolates to lower cell
					flxDn[i+1] = 0.0;
				}
			}
			Xi=Xi+1;
		}

		//// GRID CELL SPACING AND MODEL DEPTH
		for(int i=0;i<n;i++)if (W[i] < 0.0-Wtol) _error_("negative pore water generated in melt equations");

		// delete all cells with zero mass
		// adjust pore water
		for(int i=0;i<n;i++)W[i] += dW[i];

		//calculate Rsum:
		Rsum=cellsum(R,n) + flxDn[Xi];

		// delete all cells with zero mass
		D_size=0; for(int i=0;i<n;i++)if(m[i]> (0.0+Wtol))D_size++; 
		D=xNew<int>(D_size);
		D_size=0; for(int i=0;i<n;i++)if(m[i]> (0.0+Wtol)){ D[D_size] = i; D_size++;}

		celldelete(&m,n,D,D_size);
		celldelete(&W,n,D,D_size);
		celldelete(&d,n,D,D_size);
		celldelete(&dz,n,D,D_size);
		celldelete(&T,n,D,D_size);
		celldelete(&a,n,D,D_size);
		celldelete(&adiff,n,D,D_size);
		celldelete(&re,n,D,D_size);
		celldelete(&gdn,n,D,D_size);
		celldelete(&gsp,n,D,D_size);
		celldelete(&EI,n,D,D_size);
		celldelete(&EW,n,D,D_size);
		n=D_size;
		xDelete<int>(D);

		// calculate new grid lengths
		for(int i=0;i<n;i++)dz[i] = m[i] / d[i];

		/*Free resources:*/
		xDelete<IssmDouble>(R);
	}

	//calculate Fsum:
	Fsum=cellsum(F,n);

	/*Free resources:*/
	xDelete<IssmDouble>(F);

	//Manage the layering to match the user defined requirements
	managelayers(&mAdd, &dz_add, &addE, &m, &EI, &EW, &T, &d, &dz, &W, &a, &adiff, &re, &gdn, &gsp, &n, dzMin, zMax, zMin, zTop, zY);

	//// CHECK FOR MASS AND ENERGY CONSERVATION

	// calculate final mass [kg] and energy [J]
	sumER = Rsum * (LF + CtoK * CI);
	for(int i=0;i<n;i++)EI[i] = m[i] * T[i] * CI;
	for(int i=0;i<n;i++)EW[i] = W[i] * (LF + CtoK * CI);

	mSum1 = cellsum(W,n) + cellsum(m,n) + Rsum;
	sumE1 = cellsum(EI,n) + cellsum(EW,n);

	/*checks: */
	for(int i=0;i<n;i++) if (W[i]<0.0-Wtol) _error_("negative pore water generated in melt equations\n");

	/*only in forward mode! avoid round in AD mode as it is not differentiable: */
	#ifndef _HAVE_AD_
	dm = round((mSum0 - mSum1 + mAdd)*100.0)/100.0;
	dE = round(sumE0 - sumE1 - sumER +  addE - surplusE);
	if (dm !=0  || dE !=0) _error_("mass or energy are not conserved in melt equations\n"
			<< "dm: " << dm << " dE: " << dE << "\n");
	#endif

	/*Free resources:*/
	if(m)xDelete<IssmDouble>(m);
	if(EI)xDelete<IssmDouble>(EI);
	if(EW)xDelete<IssmDouble>(EW);
	if(maxF)xDelete<IssmDouble>(maxF);
	if(dW)xDelete<IssmDouble>(dW);
	if(exsW)xDelete<IssmDouble>(exsW);
	if(exsT)xDelete<IssmDouble>(exsT);
	if(surpT)xDelete<IssmDouble>(surpT);
	if(surpE)xDelete<IssmDouble>(surpE);
	if(flxDn)xDelete<IssmDouble>(flxDn);
	if(D)xDelete<int>(D);
	if(M)xDelete<IssmDouble>(M);

	/*Assign output pointers:*/
	*pMs=Msurf;
	*pM=sumM;
	*pR=Rsum;
	*pF=Fsum;
	*pmAdd=mAdd;
	*pdz_add=dz_add;

	*pT=T;
	*pd=d;
	*pdz=dz;
	*pW=W;
	*pa=a;
	*padiff=adiff;
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;
	*pn=n;

} /*}}}*/ 
void managelayers(IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble* paddE, IssmDouble** pm, IssmDouble** pEI, IssmDouble** pEW, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, IssmDouble zY){ /*{{{*/

	/*intermediary:*/
	IssmDouble* Zcum=NULL;
	IssmDouble* dzMin2=NULL;
	IssmDouble* dzMax2=NULL;
	int*        D=NULL;

	IssmDouble zY2=zY;
	IssmDouble X=0.0;
	int X1=0;
	int X2=0;
	int D_size = 0;

	IssmDouble Ztot=0.0;
	IssmDouble T_bot=0.0;
	IssmDouble m_bot=0.0;
	IssmDouble dz_bot=0.0;
	IssmDouble d_bot=0.0;
	IssmDouble W_bot=0.0;
	IssmDouble a_bot=0.0;
	IssmDouble adiff_bot=0.0;
	IssmDouble re_bot=0.0;
	IssmDouble gdn_bot=0.0;
	IssmDouble gsp_bot=0.0;
	IssmDouble EI_bot=0.0;
	IssmDouble EW_bot=0.0;
	bool       top=false;

	/*outputs:*/
	IssmDouble  mAdd = 0.0;
	IssmDouble  addE = 0.0;
	IssmDouble  dz_add = 0.0;
	IssmDouble* T=*pT;
	IssmDouble* d=*pd;
	IssmDouble* dz=*pdz;
	IssmDouble* W=*pW;
	IssmDouble* a=*pa;
	IssmDouble* adiff=*padiff;
	IssmDouble* re=*pre;
	IssmDouble* gdn=*pgdn;
	IssmDouble* gsp=*pgsp;
	IssmDouble* m=*pm;
	IssmDouble* EI=*pEI;
	IssmDouble* EW=*pEW;
	int         n=*pn;

	//Merging of cells as they are burried under snow.
	Zcum=xNew<IssmDouble>(n);
	dzMin2=xNew<IssmDouble>(n);
	dzMax2=xNew<IssmDouble>(n);

	X=0;
	Zcum[0]=dz[0]; // Compute a cumulative depth vector
	for (int i=0;i<n;i++){
		if (i==0){
			dzMin2[i]=dzMin;
		}
		else{
			Zcum[i]=Zcum[i-1]+dz[i];
			if (Zcum[i]<=zTop+Dtol){
				dzMin2[i]=dzMin;
				X=i;
			}
			else{
				dzMin2[i]=zY2*dzMin2[i-1];
			}
		}
	}

	// Check to see if any cells are too small and need to be merged
	for (int i=0; i<n; i++){
		if ( (i<=X && dz[i]<dzMin-Dtol) || (i>X && dz[i]<dzMin2[i]-Dtol) ) {

			if (i==n-1){
				X2=i;
				//find closest cell to merge with
				for(int j=n-2;j>=0;j--){
					if(m[j]!=Delflag){
						X1=j;
						break;
					}
				}
			}
			else{
				X1=i+1;
				X2=i;
			}

			// adjust variables as a linearly weighted function of mass
			IssmDouble m_new = m[X2] + m[X1];
			T[X1] = (T[X2]*m[X2] + T[X1]*m[X1]) / m_new;
			a[X1] = (a[X2]*m[X2] + a[X1]*m[X1]) / m_new;
			adiff[X1] = (adiff[X2]*m[X2] + adiff[X1]*m[X1]) / m_new;
         //use grain properties from lower cell
			re[X1] = re[X2]; 
			gdn[X1] = gdn[X2]; 
			gsp[X1] = gsp[X2]; 

			// merge with underlying grid cell and delete old cell
			dz[X1] = dz[X2] + dz[X1];                 // combine cell depths
			d[X1] = m_new / dz[X1];                   // combine top densities
			W[X1] = W[X1] + W[X2];                     // combine liquid water
			m[X1] = m_new;                             // combine top masses

			// set cell to -99999 for deletion
			m[X2] = Delflag;
		}
	}

	// delete combined cells
	D_size=0; for(int i=0;i<n;i++)if(m[i]> Delflag+Wtol)D_size++; 
	D=xNew<int>(D_size); 
	D_size=0; for(int i=0;i<n;i++)if(m[i]> Delflag+Wtol){ D[D_size] = i; D_size++;}

	celldelete(&m,n,D,D_size);
	celldelete(&W,n,D,D_size);
	celldelete(&dz,n,D,D_size);
	celldelete(&d,n,D,D_size);
	celldelete(&T,n,D,D_size);
	celldelete(&a,n,D,D_size);
	celldelete(&adiff,n,D,D_size);
	celldelete(&re,n,D,D_size);
	celldelete(&gdn,n,D,D_size);
	celldelete(&gsp,n,D,D_size);
	celldelete(&EI,n,D,D_size);
	celldelete(&EW,n,D,D_size);
	n=D_size;
	xDelete<int>(D);

	// check if any of the cell depths are too large
	X=0;
	Zcum[0]=dz[0]; // Compute a cumulative depth vector
	for (int i=0;i<n;i++){
		if (i==0){
			dzMax2[i]=dzMin*2.0;
		}
		else{
			Zcum[i]=Zcum[i-1]+dz[i];
			if (Zcum[i]<=zTop+Dtol){
				dzMax2[i]=dzMin*2.0;
				X=i;
			}
			else{
				dzMax2[i]=max(zY2*dzMin2[i-1],dzMin*2.0);
			}
		}
	}

	for (int j=n-1;j>=0;j--){
		if ((j<X && dz[j] > dzMax2[j]+Dtol) || (dz[j] > dzMax2[j]*zY2+Dtol)){
			// _printf_("dz > dzMin * 2");
			// split in two
			cellsplit(&dz, n, j,.5);
			cellsplit(&W, n, j,.5);
			cellsplit(&m, n, j,.5);
			cellsplit(&T, n, j,1.0);
			cellsplit(&d, n, j,1.0);
			cellsplit(&a, n, j,1.0);
			cellsplit(&adiff, n, j,1.0);
			cellsplit(&EI, n, j,.5);
			cellsplit(&EW, n, j,.5);
			cellsplit(&re, n, j,1.0);
			cellsplit(&gdn, n, j,1.0);
			cellsplit(&gsp, n, j,1.0);
			n++;
		}
	}

	//// CORRECT FOR TOTAL MODEL DEPTH
	// WORKS FINE BUT HAS BEEN DISABLED FOR CONVIENCE OF MODEL OUTPUT
	// INTERPRETATION
	// calculate total model depth
	Ztot = cellsum(dz,n);

	if (Ztot < zMin-Dtol){
		// printf("Total depth < zMin %f \n", Ztot);
		// mass and energy to be added
		mAdd = m[n-1]+W[n-1];
		addE = T[n-1]*m[n-1]*CI + W[n-1]*(LF+CtoK*CI);

		// add a grid cell of the same size and temperature to the bottom
		dz_bot=dz[n-1];
		T_bot=T[n-1];
		W_bot=W[n-1];
		m_bot=m[n-1];
		d_bot=d[n-1];
		a_bot=a[n-1];
		adiff_bot=adiff[n-1];
		re_bot=re[n-1];
		gdn_bot=gdn[n-1];
		gsp_bot=gsp[n-1];
		EI_bot=EI[n-1];
		EW_bot=EW[n-1];

		dz_add=dz_bot;

		newcell(&dz,dz_bot,top,n);
		newcell(&T,T_bot,top,n);
		newcell(&W,W_bot,top,n);
		newcell(&m,m_bot,top,n);
		newcell(&d,d_bot,top,n);
		newcell(&a,a_bot,top,n);
		newcell(&adiff,adiff_bot,top,n);
		newcell(&re,re_bot,top,n);
		newcell(&gdn,gdn_bot,top,n);
		newcell(&gsp,gsp_bot,top,n);
		newcell(&EI,EI_bot,top,n);
		newcell(&EW,EW_bot,top,n);
		n=n+1;
	}
	else if (Ztot > zMax+Dtol){
		// printf("Total depth > zMax %f \n", Ztot);
		// mass and energy loss
		mAdd = -(m[n-1]+W[n-1]);
		addE = -(T[n-1]*m[n-1]*CI) - (W[n-1]*(LF+CtoK*CI));
		dz_add=-(dz[n-1]);

		// remove a grid cell from the bottom
		D_size=n-1;
		D=xNew<int>(D_size);

		for(int i=0;i<n-1;i++) D[i]=i;
		celldelete(&dz,n,D,D_size);
		celldelete(&T,n,D,D_size);
		celldelete(&W,n,D,D_size);
		celldelete(&m,n,D,D_size);
		celldelete(&d,n,D,D_size);
		celldelete(&a,n,D,D_size);
		celldelete(&adiff,n,D,D_size);
		celldelete(&re,n,D,D_size);
		celldelete(&gdn,n,D,D_size);
		celldelete(&gsp,n,D,D_size);
		celldelete(&EI,n,D,D_size);
		celldelete(&EW,n,D,D_size);
		n=D_size;
		xDelete<int>(D);
	}

	/*Free resources:*/
 	xDelete<IssmDouble>(Zcum);
	xDelete<IssmDouble>(dzMin2);
	xDelete<IssmDouble>(dzMax2);
	if(D)xDelete<int>(D);

	/*Assign output pointers:*/
	*pT=T;
	*pd=d;
	*pdz=dz;
	*pW=W;
	*pa=a;
	*padiff=adiff;
	*pre=re;
	*pgdn=gdn;
	*pgsp=gsp;
	*pn=n;
	*pm=m;
	*pEI=EI;
	*pEW=EW;

	*pmAdd=mAdd;
	*paddE=addE;
	*pdz_add=dz_add;

} /*}}}*/ 
void densification(IssmDouble** pd,IssmDouble** pdz, IssmDouble* T, IssmDouble* re, int denIdx, int aIdx, int swIdx, IssmDouble adThresh, IssmDouble C, IssmDouble dt, IssmDouble Tmean, IssmDouble dIce, int m, int sid){ /*{{{*/

	//// THIS NEEDS TO BE DOUBLE CHECKED AS THERE SEAMS TO BE LITTLE DENSIFICATION IN THE MODEL OUTOUT [MAYBE COMPACTION IS COMPENSATED FOR BY TRACES OF SNOW???]

	//// FUNCTION INFO

	// Author: Alex Gardner, University of Alberta
	// Date last modified: FEB, 2008 

	// Description: 
	//   computes the densification of snow/firn using the emperical model of
	//   Herron and Langway (1980) or the semi-emperical model of Anthern et al.
	//   (2010)

	// Inputs:
	//   denIdx = densification model to use:
	//       1 = emperical model of Herron and Langway (1980)
	//       2 = semi-imerical model of Anthern et al. (2010)
	//       3 = physical model from Appendix B of Anthern et al. (2010)
	//   d   = initial snow/firn density [kg m-3]
	//   T   = temperature [K]
	//   dz  = grid cell size [m]
	//   C   = average accumulation rate [kg m-2 yr-1]
	//   dt  = time lapsed [s]
	//   re  = effective grain radius [mm];
	//   Ta  = mean annual temperature                                          

	// Reference: 
	// Herron and Langway (1980), Anthern et al. (2010)

	//// FOR TESTING
	// denIdx = 2;
	// d = 800;
	// T = 270;
	// dz = 0.005;
	// C = 200;
	// dt = 60*60;
	// re = 0.7;
	// Tmean = 273.15-18;

	//// MAIN FUNCTION
	// specify constants
	dt      = dt / dts;  // convert from [s] to [d]
	// R     = 8.314        // gas constant [mol-1 K-1]
	// Ec    = 60           // activation energy for self-diffusion of water
	//                      // molecules through the ice tattice [kJ mol-1]
	// Eg    = 42.4         // activation energy for grain growth [kJ mol-1]

	/*intermediary: */
	IssmDouble c0=0.0;
	IssmDouble c1=0.0;
   IssmDouble c2=0.0;
	IssmDouble H=0.0;
	IssmDouble M0=0.0;
	IssmDouble M1=0.0;
   IssmDouble M2=0.0;
	IssmDouble c0arth=0.0;
	IssmDouble c1arth=0.0;

	/*output: */
	IssmDouble* dz=NULL;
	IssmDouble* d=NULL;

	if(VerboseSmb() && sid==0 && IssmComm::GetRank()==0)_printf0_("   densification module\n");

	/*Recover pointers: */
	dz=*pdz;
	d=*pd;

	// initial mass
	IssmDouble* mass_init = xNew<IssmDouble>(m);for(int i=0;i<m;i++) mass_init[i]=d[i] * dz[i];

	/*allocations and initialization of overburden pressure and factor H: */
	IssmDouble* cumdz = xNew<IssmDouble>(m-1);
	cumdz[0]=dz[0];
	for(int i=1;i<m-1;i++)cumdz[i]=cumdz[i-1]+dz[i];

	IssmDouble* obp = xNew<IssmDouble>(m);
	obp[0]=0.0;
	for(int i=1;i<m;i++)obp[i]=cumdz[i-1]*d[i-1];

	// calculate new snow/firn density for:
	//   snow with densities <= 550 [kg m-3]
	//   snow with densities > 550 [kg m-3]

	for(int i=0;i<m;i++){
		switch (denIdx){
			case 1: // Herron and Langway (1980)
				c0 = (11.0 * exp(-10160.0 / (T[i] * R))) * C/1000.0;
				c1 = (575.0 * exp(-21400.0 / (T[i]* R))) * pow(C/1000.0,.5);
				break;

			case 2: // Arthern et al. (2010) [semi-emperical]
				// common variable
				// NOTE: Ec=60000, Eg=42400 (i.e. should be in J not kJ)
				H = exp((-60000.0/(T[i] * R)) + (42400.0/(Tmean * R))) * (C * 9.81);
				c0 = 0.07 * H;
				c1 = 0.03 * H;
				break;

			case 3: // Arthern et al. (2010) [physical model eqn. B1]

				// common variable
				H = exp((-60000.0/(T[i] * R))) * obp[i] / pow(re[i]/1000.0,2.0);
				c0 = 9.2e-9 * H;
				c1 = 3.7e-9 * H;
				break;

			case 4: // Li and Zwally (2004)
				c0 = (C/dIce) * max(139.21 - 0.542*Tmean,1.0)*8.36*pow(max(CtoK - T[i],1.0),-2.061);
				c1 = c0;
				break;

			case 5: // Helsen et al. (2008)
				// common variable
				c0 = (C/dIce) * max(76.138 - 0.28965*Tmean,1.0)*8.36*pow(max(CtoK - T[i],1.0),-2.061);
				c1 = c0;
				break;

			case 6: // Ligtenberg and others (2011) [semi-emperical], Antarctica 
				// common variable
				// From literature: H = exp((-60000.0/(Tmean * R)) + (42400.0/(Tmean * R))) * (C * 9.81);
				H = exp((-60000.0/(T[i] * R)) + (42400.0/(Tmean * R))) * (C * 9.81);
				c0arth = 0.07 * H;
				c1arth = 0.03 * H;
				//ERA-5 old
				//M0 = max(2.3128 - (0.2480 * log(C)),0.25);
				//M1 = max(2.7950 - (0.3318 * log(C)),0.25);
				// ERA5 new aIdx=1, swIdx=0
				if (aIdx==1 && swIdx==0){
					if (fabs(adThresh - 820.0) < Dtol){
                  // ERA5 v4
                  M0 = max(1.5131 - (0.1317 * log(C)),0.25);
                  M1 = max(1.8819 - (0.2158 * log(C)),0.25);
					}
					else{
						// ERA5 new aIdx=1, swIdx=0
						//M0 = max(1.8785 - (0.1811 * log(C)),0.25);
						//M1 = max(2.0005 - (0.2346 * log(C)),0.25);
						// ERA5 new aIdx=1, swIdx=0, bare ice
						M0 = max(1.8422 - (0.1688 * log(C)),0.25);
						M1 = max(2.4979 - (0.3225 * log(C)),0.25);
					}
				}
				// ERA5 new aIdx=2, swIdx=1
				else if (aIdx<3 && swIdx>0){
					M0 = max(2.2191 - (0.2301 * log(C)),0.25);
					M1 = max(2.2917 - (0.2710 * log(C)),0.25);
				}
				// ERA5 new aIdx=2, swIdx=0
				//else if (aIdx==2){
				//}
				//From Ligtenberg
				//H = exp((-60000.0/(Tmean * R)) + (42400.0/(Tmean * R))) * (C * 9.81);
				//M0 = max(1.435 - (0.151 * log(C)),0.25);
				//M1 = max(2.366 - (0.293 * log(C)),0.25);
            //RACMO
            M0 = max(1.6383 - (0.1691 * log(C)),0.25);
            M1 = max(1.9991 - (0.2414 * log(C)),0.25);
				c0 = M0*c0arth;
				c1 = M1*c1arth;
				c2 = M2*c1arth;
				break;

			case 7: // Kuipers Munneke and others (2015) [semi-emperical], Greenland
				// common variable
				// From literature: H = exp((-60000.0/(T[i] * R)) + (42400.0/(T[i] * R))) * (C * 9.81);
				H = exp((-60000.0/(T[i] * R)) + (42400.0/(Tmean * R))) * (C * 9.81);
				c0arth = 0.07 * H;
				c1arth = 0.03 * H;
				// ERA5 old
				//M0 = max(1.8554 - (0.1316 * log(C)),0.25);
				//M1 = max(2.8901 - (0.3014 * log(C)),0.25);
				// ERA5 new aIdx=1, swIdx=0
				if (aIdx==1 && swIdx==0){
					if (fabs(adThresh - 820.0) < Dtol){
						// ERA5 v4
						M0 = max(1.3566 - (0.1350 * log(C)),0.25);
						M1 = max(1.8705 - (0.2290 * log(C)),0.25);
					}
					else{
						// ERA5 new aIdx=1, swIdx=0
						//M0 = max(1.4574 - (0.1123 * log(C)),0.25);
						//M1 = max(2.0238 - (0.2070 * log(C)),0.25);
						// ERA5 new aIdx=1, swIdx=0, bare ice
						M0 = max(1.4318 - (0.1055 * log(C)),0.25);
						M1 = max(2.0453 - (0.2137 * log(C)),0.25);
					}
				}
				// ERA5 new aIdx=2, swIdx=1
				else if (aIdx<3 && swIdx>0){
					M0 = max(1.7834 - (0.1409 * log(C)),0.25);
					M1 = max(1.9260 - (0.1527 * log(C)),0.25);
				}
				// ERA5 new aIdx=2, swIdx=0
				//else if (aIdx==2){
				//}
				// From Kuipers Munneke
				//M0 = max(1.042 - (0.0916 * log(C)),0.25);
				//M1 = max(1.734 - (0.2039 * log(C)),0.25);
				// RACMO
            M0 = max(1.2691 - (0.1184 * log(C)),0.25);
            M1 = max(1.9983 - (0.2511 * log(C)),0.25);
            c0 = M0*c0arth;
            c1 = M1*c1arth;
            c2 = M2*c1arth;
				break;
		}

      // new snow density
      if(d[i] <= 550.0+Dtol) d[i] = d[i] + (c0 * (dIce - d[i]) / 365.0 * dt);
      else if(d[i] <= 830.0+Dtol | fabs(c2)<Dtol) d[i] = d[i] + (c1 * (dIce - d[i]) / 365.0 * dt);
      else d[i] = d[i] + (c2 * (dIce - d[i]) / 365.0 * dt);

		// do not allow densities to exceed the density of ice
		if(d[i] > dIce-Ptol) d[i]=dIce;

		// calculate new grid cell length
		dz[i] = mass_init[i] / d[i];
	}

	/*Free resources:*/
	xDelete<IssmDouble>(mass_init);
	xDelete<IssmDouble>(cumdz);
	xDelete<IssmDouble>(obp);

	/*Assign output pointers:*/
	*pdz=dz;
	*pd=d;

} /*}}}*/
