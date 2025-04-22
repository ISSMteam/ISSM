/*!\file:  SurfaceMassBalancex.h
 * \brief header file for SMB
 */

#ifndef _SurfaceMassBalancex_H
#define _SurfaceMassBalancex_H

#include "../../classes/classes.h"

/* local prototypes: */
void SurfaceMassBalancex(FemModel* femmodel);
void SmbForcingx(FemModel* femmodel);
void SmbGradientsx(FemModel* femmodel);
void SmbGradientsElax(FemModel* femmodel);
void Smbarmax(FemModel* femmodel);
void Delta18oParameterizationx(FemModel* femmodel);
void MungsmtpParameterizationx(FemModel* femmodel);
void Delta18opdParameterizationx(FemModel* femmodel);
void PositiveDegreeDayx(FemModel* femmodel);
void PositiveDegreeDaySicopolisx(FemModel* femmodel);
void PositiveDegreeDayGCMx(FemModel* femmodel);
void SmbHenningx(FemModel* femmodel);
void SmbComponentsx(FemModel* femmodel);
void SmbMeltComponentsx(FemModel* femmodel);
void SmbGradientsComponentsx(FemModel* femmodel);
void SmbDebrisEvattx(FemModel* femmodel);
/* SEMIC: */
void SmbSemicx(FemModel* femmodel, int ismethod);
/*GEMB: */
void       Gembx(FemModel* femmodel);
void       GembgridInitialize(IssmDouble** pdz, int* psize, IssmDouble z_top, IssmDouble dz_top, IssmDouble z_max, IssmDouble beta);
IssmDouble Marbouty(IssmDouble T, IssmDouble d, IssmDouble dT);
IssmDouble gardnerAlb(IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, IssmDouble dPHC, int m);
void grainGrowth(IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, IssmDouble* T,IssmDouble* dz,IssmDouble* d, IssmDouble* W,IssmDouble smb_dt,int m,int aIdx, int sid);
void albedo(IssmDouble** a, IssmDouble** adiff, int aIdx, IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble cldFrac, IssmDouble aIce, IssmDouble aSnow, IssmDouble aValue, IssmDouble adThresh, IssmDouble* T, IssmDouble* W, IssmDouble P, IssmDouble EC, IssmDouble Msurf, IssmDouble clabSnow, IssmDouble clabIce, IssmDouble SZA, IssmDouble COT, IssmDouble t0wet, IssmDouble t0dry, IssmDouble K, IssmDouble dt, IssmDouble dIce, int m, int sid);
void shortwave(IssmDouble** pswf, int swIdx, int aIdx, IssmDouble dsw, IssmDouble dswdiff, IssmDouble as, IssmDouble asdiff, IssmDouble* d, IssmDouble* dz, IssmDouble* re, IssmDouble dIce, int m, int sid);
void thermo(IssmDouble* pshf, IssmDouble* plhf, IssmDouble* pEC, IssmDouble** pT, IssmDouble* pulwrf, IssmDouble* re, IssmDouble* dz, IssmDouble* d, IssmDouble* swf, IssmDouble dlw, IssmDouble Ta, IssmDouble V, IssmDouble eAir, IssmDouble pAir, int tcIdx, int eIdx, IssmDouble teValue, IssmDouble dulwrfValue, IssmDouble teThresh, IssmDouble Ws, IssmDouble dt0, IssmDouble dzMin, int m, IssmDouble Vz, IssmDouble Tz, IssmDouble thermo_scaling, IssmDouble dIce, int sid, bool isconstrainsurfaceT, bool isdeltaLWup);
void accumulation(IssmDouble** pT, IssmDouble** pdz, IssmDouble** pd, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, IssmDouble* pRa, int* pm, int aIdx, int dsnowIdx, IssmDouble Tmean, IssmDouble Ta, IssmDouble P, IssmDouble dzMin, IssmDouble aSnow, IssmDouble C, IssmDouble V, IssmDouble Vmean, IssmDouble dIce, int sid);
void melt(IssmDouble* pM, IssmDouble* pMs, IssmDouble* pR, IssmDouble* pF, IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble Ra, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, IssmDouble zY, IssmDouble dIce, int sid);
void managelayers(IssmDouble* pmAdd, IssmDouble* pdz_add, IssmDouble* paddE, IssmDouble** pm, IssmDouble** pEI, IssmDouble** pEW, IssmDouble** pT, IssmDouble** pd, IssmDouble** pdz, IssmDouble** pW, IssmDouble** pa, IssmDouble** padiff, IssmDouble** pre, IssmDouble** pgdn, IssmDouble** pgsp, int* pn, IssmDouble dzMin, IssmDouble zMax, IssmDouble zMin, IssmDouble zTop, IssmDouble zY);
void densification(IssmDouble** pd,IssmDouble** pdz, IssmDouble* T, IssmDouble* re, int denIdx, int aIdx, int swIdx, IssmDouble adThresh, IssmDouble C, IssmDouble dt, IssmDouble Tmean, IssmDouble dIce, int m, int sid);
#endif  /* _SurfaceMassBalancex_H*/
