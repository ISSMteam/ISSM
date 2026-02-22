/* file:  PddSurfaceMassBlanceSicopolis.cpp
   Calculating the surface mass balance using the adapted PDD routine from SICOPOLIS.
 */

#include "./elements.h"
#include "../Numerics/numerics.h"
#include "../Exceptions/exceptions.h"
#include <cmath>

IssmDouble PddSurfaceMassBalanceSicopolis(IssmDouble* monthlytemperatures, IssmDouble* monthlyprec,
				 IssmDouble* melt, IssmDouble* accu, IssmDouble* melt_star, IssmDouble* t_ampl, IssmDouble* p_ampl,
				 IssmDouble yts, IssmDouble s, IssmDouble desfac,
				 IssmDouble s0t, IssmDouble s0p, IssmDouble rlaps,
				 IssmDouble rho_water,IssmDouble rho_ice,IssmDouble pdd_fac_ice,IssmDouble pdd_fac_snow){

  int			imonth;				// month counter
  IssmDouble B;					// output: surface mass balance (m/a IE), melt+accumulation
  IssmDouble frac_solid, snowfall, rainfall, runoff; 
  IssmDouble saccu;				// yearly surface accumulation (m/a IE)
  IssmDouble smelt;				// yearly melt (m/a IE)
  IssmDouble smelt_star;		// yearly ...
  IssmDouble precip;				// total precipitation during 1 year
  IssmDouble sconv;				//rhow_rain/rhoi / 12 months
  IssmDouble st;					// elevation between altitude of the temp record and current altitude
  IssmDouble sp;					// elevation between altitude of the prec record and current altitude
  IssmDouble q;					// q is desert/elev. fact
  IssmDouble pdd;					// pdd factor (a * degC)
  IssmDouble tstar;				// monthly temp. after lapse rate correction (degC)
  IssmDouble precip_star;		// monthly precip after correction (m/a IE)
  IssmDouble Pmax = 0.6;
  IssmDouble inv_twelve=1./12.; 

  sconv=(rho_water/rho_ice);		//rhow_rain/rhoi

  pdd_fac_snow=pdd_fac_snow*(0.001*365)*sconv; // (mm WE)/(d*deg C) --> (m IE)/(a*deg C)
  pdd_fac_ice=pdd_fac_ice*(0.001*365)*sconv; // (mm WE)/(d*deg C) --> (m IE)/(a*deg C)

  /* initalize fields */
  precip=0.0;
  tstar=0.0;
  snowfall=0.0;
  pdd=0.0;
  /* seasonal loop */
  for(imonth=0;imonth<12;imonth++){

    /********* Surface temperature correction *******/    
    st=(s-s0t)/1000.;

    /******** Monhtly temperature correction *******/
    monthlytemperatures[imonth]=monthlytemperatures[imonth]-rlaps*st;//*max(st,1e-3);
    tstar=monthlytemperatures[imonth]+t_ampl[0];

    /********* Precipitation correction *************/
    /* Ref: Vizcaino et al 2010; DOI 10.1007/s00382-009-0591-y */
    if(s0p<2000.0)
      q=exp(desfac*(max(s,2000.0)-2000.0));
	 else
      q=exp(desfac*(max(s,2000.0)-s0p));

    precip_star=q*monthlyprec[imonth]*sconv*p_ampl[0]*yts; // convert precip from m/s -> m/a
    precip=precip+precip_star*inv_twelve;

    /********* compute PDD **************************/
    /* Ref: Calov & Greve 2005 Journal of Glaciology, Vol. 51, No. 172, 2005, Correspondence */
	 IssmDouble s_stat=5.0;
    IssmDouble inv_sqrt2pi =1.0/sqrt(2.0*PI);
    IssmDouble inv_s_stat  =1.0/s_stat;
    IssmDouble inv_sqrt2   =1.0/sqrt(2.0);

	 #if !defined(_HAVE_ADOLC_)
    pdd=pdd+(s_stat*inv_sqrt2pi*exp(-0.5*pow(tstar*inv_s_stat,2))
				 +0.5*tstar*erfc(-tstar*inv_s_stat*inv_sqrt2))*inv_twelve;
	 #else
	 _error_("Cannot differentiate erfc, talk to ADOLC folks (http://functions.wolfram.com/GammaBetaErf/Erfc/20/01/)");
	 #endif

	 /*Partition of precip in solid and liquid parts, Bales et al. (2009) */
	 IssmDouble temp_rain=7.2;		// Threshold monthly mean temperature for
											// precipitation = 101% rain, in deg C
	 IssmDouble temp_snow=-11.6;  // Threshold monthly mean temperature for
											// precipitation = 100% snow, in deg C

	 IssmDouble coeff1=5.4714e-01;	// Coefficients
	 IssmDouble coeff2=-9.1603e-02;	// of
	 IssmDouble coeff3=-3.314e-03;	// the
	 IssmDouble coeff4= 4.66e-04;		// fifth-order
	 IssmDouble coeff5=3.8e-05;		// polynomial
	 IssmDouble coeff6=6.0e-07;		// fit

	 if(tstar>=temp_rain)
	  frac_solid = 0.0;
	 else if(tstar<=temp_snow)
	  frac_solid = 1.0;
	 else{ 
		 frac_solid=coeff1+tstar*(coeff2
					 +tstar*(coeff3+tstar*(coeff4+tstar*(coeff5+tstar*coeff6))));
	 }

	 snowfall=snowfall+precip_star*frac_solid*inv_twelve;
  } 
  /* end of seasonal loop */ 

  rainfall=precip-snowfall;
  if(snowfall<0.0) snowfall=0.0;   // correction of
  if(rainfall<0.0) rainfall=0.0;   // negative values

  if(rainfall<=(Pmax*snowfall)){
	  if((rainfall+pdd_fac_snow*pdd)<=(Pmax*snowfall)) {
		  smelt_star = rainfall+pdd_fac_snow*pdd;
		  smelt      = 0.0;
		  runoff     = smelt;
	  }
	  else{
		  smelt_star = Pmax*snowfall;
		  smelt      = pdd_fac_ice*(pdd-(smelt_star-rainfall)/pdd_fac_snow);
		  runoff     = smelt;
	  }
  } 
  else{
	  smelt_star = Pmax*snowfall;
	  smelt      = pdd_fac_ice*pdd;
	  runoff     = smelt+rainfall-Pmax*snowfall;
  }

  saccu = precip;	

  /* asign output*/
  melt[0]=runoff/yts;
  accu[0]=saccu/yts;
  melt_star[0]=smelt_star/yts;
  B=(saccu-runoff)/yts;

  return B;
}
