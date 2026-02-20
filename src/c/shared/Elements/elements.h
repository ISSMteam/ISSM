/*!\file: elements.h
 * \brief prototypes for elements.h
 */ 

#ifndef _SHARED_ELEMENTS_H_
#define _SHARED_ELEMENTS_H_

#include "../Numerics/types.h"

IssmDouble Cuffey(IssmDouble temperature);
IssmDouble BuddJacka(IssmDouble temperature);
IssmDouble CuffeyTemperate(IssmDouble temperature, IssmDouble waterfraction, IssmDouble stressexp);
IssmDouble Paterson(IssmDouble temperature);
IssmDouble Arrhenius(IssmDouble temperature,IssmDouble depth,IssmDouble n);
IssmDouble NyeH2O(IssmDouble temperature);
IssmDouble NyeCO2(IssmDouble temperature);
IssmDouble LliboutryDuval(IssmDouble enthalpy, IssmDouble pressure, IssmDouble n, IssmDouble betaCC, IssmDouble referencetemperature, IssmDouble heatcapacity, IssmDouble latentheat);
// IssmDouble LliboutryDuval(IssmDouble temperature, IssmDouble waterfraction, IssmDouble depth,IssmDouble n);
IssmDouble EstarLambdaS(IssmDouble epseff, IssmDouble epsprime_norm);
void EstarOmega(IssmDouble* omega,IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble vmag,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz, IssmDouble* dvmag);
void EstarStrainrateQuantities(IssmDouble *pepsprime_norm, IssmDouble vx,IssmDouble vy,IssmDouble vz,IssmDouble vmag,IssmDouble* dvx,IssmDouble* dvy,IssmDouble* dvz,IssmDouble* dvmag);
IssmDouble PddSurfaceMassBalance(IssmDouble* monthlytemperatures,  IssmDouble* monthlyprec,
				 IssmDouble* pdds, IssmDouble* pds, IssmDouble* melt, IssmDouble* accu, IssmDouble signorm, 
				 IssmDouble yts, IssmDouble h, IssmDouble s, IssmDouble desfac,IssmDouble s0t,
				 IssmDouble s0p, IssmDouble rlaps, IssmDouble rlapslgm,
				 IssmDouble TdiffTime,IssmDouble sealevTime,IssmDouble pddsnowfac,IssmDouble pddicefac,
				 IssmDouble rho_water, IssmDouble rho_ice);
IssmDouble PddSurfaceMassBalanceSicopolis(IssmDouble* monthlytemperatures,  IssmDouble* monthlyprec,
				 IssmDouble* melt, IssmDouble* accu, IssmDouble* melt_star, IssmDouble* t_ampl, IssmDouble* p_ampl,
				 IssmDouble yts, IssmDouble s, IssmDouble desfac,IssmDouble s0t,
				 IssmDouble s0p, IssmDouble rlaps, IssmDouble rho_water, IssmDouble rho_ice, IssmDouble pdd_fac_ice, IssmDouble pdd_fac_snow);
void ComputeDelta18oTemperaturePrecipitation(IssmDouble Delta18oSurfacePresent, IssmDouble Delta18oSurfaceLgm, IssmDouble Delta18oSurfaceTime,
					     IssmDouble Delta18oPresent, IssmDouble Delta18oLgm, IssmDouble Delta18oTime,
					     IssmDouble* PrecipitationsPresentday,
					     IssmDouble* TemperaturesLgm, IssmDouble* TemperaturesPresentday, 
					     IssmDouble* monthlytemperaturesout, IssmDouble* monthlyprecout);
void ComputeMungsmTemperaturePrecipitation(IssmDouble TdiffTime, IssmDouble PfacTime,
					   IssmDouble* PrecipitationsLgm,IssmDouble* PrecipitationsPresentday,
					   IssmDouble* TemperaturesLgm, IssmDouble* TemperaturesPresentday, 
					   IssmDouble* monthlytemperaturesout, IssmDouble* monthlyprecout);
void ComputeD18OTemperaturePrecipitationFromPD(IssmDouble d018,IssmDouble dpermil,bool isTemperatureScaled,
			bool isPrecipScaled, IssmDouble f, IssmDouble* PrecipitationPresentday,IssmDouble* TemperaturePresentday,
			IssmDouble* PrecipitationReconstructed,IssmDouble* TemperatureReconstructed, IssmDouble* monthlytemperaturesout,
			IssmDouble* monthlyprecout);
IssmDouble DrainageFunctionWaterfraction(IssmDouble waterfraction, IssmDouble dt=0.);
IssmDouble StressIntensityIntegralWeight(IssmDouble depth, IssmDouble water_depth, IssmDouble thickness);

/*Enthalphy*/
void EnthalpyToThermal(IssmDouble* ptemperature,IssmDouble* pwaterfraction,IssmDouble enthalpy,IssmDouble pressure);

/*Print arrays*/
void printarray(IssmPDouble* array,int lines,int cols=1);
#if _HAVE_AD_  && !defined(_WRAPPERS_)
void printarray(IssmDouble* array,int lines,int cols=1);
#endif
void printarray_matlab(const char* filename,int* array,int lines,int cols=1);
void printarray(int* array,int lines,int cols=1);
void printarray(bool* array,int lines,int cols=1);
void printsparsity(IssmPDouble* array,int lines,int cols=1);
void printbinary(int n);
void InversionStatsHeader(int NJ);
void InversionStatsIter(int iter,double J, double Gnorm, double* Jlist, int N);
void InversionStatsFooter(int NJ);
#endif //ifndef _SHARED_ELEMENTS_H_
