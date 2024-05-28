/* file:  PddSurfaceMassBlance.cpp
   Calculating the surface mass balance using the positive degree day method.
   Updating the precipitation and temperature to the new elevation
 */

#include "../io/io.h"
#include "./elements.h"
#include "../Numerics/numerics.h"
#include <cmath>

IssmDouble PddSurfaceMassBalance(IssmDouble* monthlytemperatures, IssmDouble* monthlyprec,
         IssmDouble* pdds, IssmDouble* pds, IssmDouble* melt, IssmDouble* accu, 
         IssmDouble signorm, IssmDouble yts, IssmDouble h, IssmDouble s, IssmDouble desfac,
         IssmDouble s0t,IssmDouble s0p, IssmDouble rlaps,IssmDouble rlapslgm,
         IssmDouble TdiffTime,IssmDouble sealevTime, IssmDouble pddsnowfac,IssmDouble pddicefac,
         IssmDouble rho_water,IssmDouble rho_ice){

  // output:
  IssmDouble B;    // surface mass balance, melt+accumulation
  int    iqj,imonth;

  IssmDouble saccu;     // yearly surface accumulation
  IssmDouble smelt;     // yearly melt
  IssmDouble precrunoff;      // yearly runoff
  IssmDouble prect; // total precipitation during 1 year taking into account des. ef.
  IssmDouble water; //water=rain + snowmelt 
  IssmDouble runoff; //meltwater only, does not include rain 
  IssmDouble sconv; //rhow_rain/rhoi / 12 months

  //IssmDouble  sealev=0.;         // degrees per meter. 7.5 lev's 99 paper, 9 Marshall 99 paper
  //IssmDouble  Pfac=0.5,Tdiff=0.5;
  IssmDouble rtlaps;
  // IssmDouble lapser=6.5         // lapse rate
  // IssmDouble desfac = 0.3;      // desert elevation factor
  // IssmDouble s0p=0.;            // should be set to elevation from precip source
  // IssmDouble s0t=0.;         // should be set to elevation from temperature source
  IssmDouble st;             // elevation between altitude of the temp record and current altitude
  IssmDouble sp;             // elevation between altitude of the prec record and current altitude
  IssmDouble deselcut=1.0;

  // PDD and PD constants and variables
  IssmDouble siglim;          // sigma limit for the integration which is equal to 2.5 sigmanorm
  IssmDouble signormc = signorm - 0.5;     // sigma of the temperature distribution for cloudy day
  IssmDouble siglimc, siglim0, siglim0c;
  IssmDouble PDup, pddsig, PDCUT = 2.0; // PDcut: rain/snow cutoff temperature (C)
  IssmDouble DT = 0.02;
  IssmDouble pddt, pd; // pd: snow/precip fraction, precipitation falling as snow

  IssmDouble q, qmpt;   // q is desert/elev. fact, hnpfac is huybrect fact, and pd is normal dist.
  IssmDouble qm = 0.;   // snow part of the precipitation 
  IssmDouble qmt = 0.;  // precipitation without desertification effect adjustment
  IssmDouble qmp = 0.;  // desertification taken into account
  IssmDouble pdd = 0.;     
  IssmDouble frzndd = 0.;  

  IssmDouble tstar;          // monthly mean surface temp
  IssmDouble Tsum= 0.;       // average summer (JJA) temperature
  IssmDouble Tsurf = 0.;     // average annual temperature    

  IssmDouble deltm=1./12.;
  int    ismon[12]={11,0,1,2,3,4,5,6,7,8,9,10};

  IssmDouble snwm;  // snow that could have been melted in a year.
  IssmDouble snwmf; //  ablation factor for snow per positive degree day.
  IssmDouble smf;   //  ablation factor for ice per pdd (Braithwaite 1995 from tarasov 2002).

  IssmDouble dfrz=1.5, CovrLm=2009./3.35e+5, dCovrLm=dfrz*CovrLm; //m*J kg^-1 C^-1 /(J kg^-1)=m/C yr
  IssmDouble supice,supcap,diffndd;
  IssmDouble fsupT=0.5,  fsupndd=0.5;  // Tsurf mode factors for supice
  IssmDouble pddtj, hmx2;
  IssmDouble pddsnowfac0=4.3, pddicefac0=8.3;
  IssmDouble snowfac, icefac;

  sconv=(rho_water/rho_ice)/12.; //rhow_rain/rhoi / 12 months

  /*PDD constant*/
  siglim = 2.5*signorm; 
  siglimc = 2.5*signormc;
  siglim0 = siglim/DT + 0.5;
  siglim0c = siglimc/DT + 0.5;
  PDup = siglimc+PDCUT;

  // seasonal loop
  for (iqj = 0; iqj < 12; iqj++){
    imonth =  ismon[iqj];

    /*********compute lapse rate ****************/
    st=(s-s0t)/1000.;
    rtlaps=TdiffTime*rlapslgm + (1.-TdiffTime)*rlaps; // lapse rate

    /*********compute Surface temperature *******/
    monthlytemperatures[imonth]=monthlytemperatures[imonth] - rtlaps *max(st,sealevTime*0.001);
    tstar = monthlytemperatures[imonth];
    Tsurf = tstar*deltm+Tsurf;        

    /*********compute PD ****************/
    if (tstar < PDup){
      pd = 1.;
      if (tstar >= -siglimc){ pd = pds[reCast<int,IssmDouble>(tstar/DT + siglim0c)];}}
    else { 
      pd = 0.;}

    /******exp des/elev precip reduction*******/
    sp=(s-s0p)/1000.-deselcut; // deselev effect is wrt chng in topo
    if (sp>0.0){q = exp(-desfac*sp);}
    else {q = 1.0;}

    qmt= qmt + monthlyprec[imonth]*sconv;  //*sconv to convert in m of ice equivalent per month  
    qmpt= q*monthlyprec[imonth]*sconv;
    qmp= qmp + qmpt;
    qm= qm + qmpt*pd;

    /*********compute PDD************/
    // ndd(month)=-(tstar-pdd(month)) since ndd+pdd gives expectation of
    // gaussian=T_m, so ndd=-(Tsurf-pdd)
    if (iqj>5 && iqj<9){ Tsum=Tsum+tstar;} 

    if (tstar >= siglim) {pdd = pdd + tstar*deltm;}
    else if (tstar> -siglim){
      pddsig=pdds[reCast<int,IssmDouble>(tstar/DT + siglim0)];
      pdd = pdd + pddsig*deltm;
      frzndd = frzndd - (tstar-pddsig)*deltm;}
    else{frzndd = frzndd - tstar*deltm; }

    /*Assign output pointer*/
    *(monthlytemperatures+imonth) = monthlytemperatures[imonth];
    *(monthlyprec+imonth) = monthlyprec[imonth];      
  } // end of seasonal loop 
  //******************************************************************

  saccu = qm;
  prect = qmp;     // total precipitation during 1 year taking into account des. ef.
  Tsum=Tsum/3;

  snowfac=pddsnowfac0;
  icefac=pddicefac0;
  if (pddsnowfac>0) {
    if (pddsnowfac<1.65) {
      _printf0_("WARNING: Pdd snow factor input, " << pddsnowfac << ", results in a negative value. It will be ignored. \n");
    }
    else{
    snowfac=pddsnowfac;
    }
  }
  if (pddicefac>0) {
    if (pddicefac>17.22) {
      _printf0_("WARNING: Pdd ice factor input, " << pddicefac << ", results in a negative value. It will be ignored. \n");
    }
    else{
      icefac=pddicefac;
    }
  }

  /***** determine PDD factors *****/
  if(Tsum< -1.) {
    snwmf=(2.65+snowfac-pddsnowfac0)*0.001;   //  ablation factor for snow per positive degree day.*0.001 to go from mm to m/ppd
    smf=17.22*0.001;    //  ablation factor for ice per pdd (Braithwaite 1995 from tarasov 2002)
  } 
  else if(Tsum< 10){
    snwmf = (0.15*(Tsum+1) + (2.65+snowfac-pddsnowfac0))*0.001;
    smf = (((17.22-icefac)/(pow(11,3)))*pow((10.-Tsum),3) + pddicefac0)*0.001;
    //JC,smf = (((icefac-pddicefac0)/(Tsum+1))*pow((10.-Tsum),3) + pddicefac0)*0.001;
  }
  else{
    snwmf=snowfac*0.001;
    smf=icefac*0.001;
  }
  //snwmf=0.95*snwmf;
  //smf=0.95*smf;

  /*****  compute PDD ablation and refreezing *****/
  pddt = pdd *365;
  snwm = snwmf*pddt;       // snow that could have been melted in a year
  hmx2 = min(h,dfrz);   // refreeze active layer max depth: dfrz

  if(snwm < saccu) {
    water=prect-saccu + snwm; //water=rain + snowmelt
    //     l 2.2= capillary factor
    //     Should refreezing be controlled by frzndd or by mean annual Tsurf?
    //     dCovrLm concept is of warming of active layer (thickness =d@=1-
    //     >2m)
    //     problem with water seepage into ice: should be sealed after
    //     refreezing
    //     so everything needs to be predicated on 1 year scale, except for
    //     thermal
    //     conductivity through ice
    //     also, need to account that melt season has low accum, so what's
    //     going to
    //     hold the meltwater around for refreezing? And melt-time will have
    //     low seasonal frzndd

    //      Superimposed ice :  Pfeffer et al. 1991, Tarasov 2002

    supice= min(hmx2*CovrLm*frzndd+2.2*(saccu-snwm), water); // superimposed ice
    supcap=min(2.2*(saccu-snwm),water);
    runoff=snwm - supice;  //meltwater only, does not include rain
  }
  else {  //all snow melted
    supice= min(hmx2*CovrLm*frzndd, prect );
    runoff= saccu + smf*(pddt-saccu/snwmf) - supice;
    supcap=0;
  }
  //     pdd melting doesn't cool Tsurf, so ndd refreeze shouldn't warm it
  //     except pdd melt heat source is atmosphere, while refreeze is
  //     ground/ice stored interim
  //     assume pdd=ndd, then melt should equal refreeze and Tsurf should=0
  //     assume ndd=2*pdd, then all supice is refrozen, but Tsurf should be
  //     <0
  //     assume ndd>pdd, little melt => little supice 
  //     bottom line: compare for Tsurf<0 : supice and no supice case,
  //     expect Tsurf difference
  //     except some of cooling flux comes from atmosphere//
  //     1 dm supice should not raise Tsurf by 1/dCovrLm = 16.675C
  //     does supice make sense when H< 0.1m? then d=thermoactive ice layer ////
  //     < 0.1 

  //     make more sense to just use residual pdd-ndd except that pdd
  //     residual not clear yet
  //     frzndd should not be used up by refreezing in snow, so stick in
  //     supcap.
  diffndd=0;
  if (frzndd>0) {
    diffndd=fsupndd*min((supice-supcap)/dCovrLm ,frzndd);
    frzndd=frzndd-diffndd;
  }
  if(runoff<0){
    saccu= saccu -runoff;
    smelt = 0;
    precrunoff=prect-saccu;
    //here assume pdd residual is 0, => 
    Tsurf= max(Tsurf,-frzndd);
  }
  else {
    smelt = runoff;
    precrunoff=prect-max(0.,supice)-saccu;}
  //here really need pdd balance, try 0.5 fudge factor?
  //at least runoff>0 => it's fairly warm, so Tsurf is !<<0,
  //yet from site plots, can be ice free with Tsurf=-5.5C
  if(Tsurf<0) {
    Tsurf= min(Tsurf+fsupT*diffndd , 0.);}

  melt[0]=smelt/yts;
  accu[0]=saccu/yts;
  B = saccu - smelt;
  B = B/yts;
  pddtj=pddt;

  return B;
}
