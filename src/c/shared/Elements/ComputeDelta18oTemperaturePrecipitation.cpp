/* file:  ComputeTemperaturePrecipitation.cpp
   Compute the temperature and precipitation at time t from 
   the data at present day and delta18O
 */

#include "./elements.h"
#include "../Numerics/numerics.h"
#include <cmath>

void ComputeDelta18oTemperaturePrecipitation(IssmDouble Delta18oSurfacePresent, IssmDouble Delta18oSurfaceLgm, IssmDouble Delta18oSurfaceTime,
				     IssmDouble Delta18oPresent, IssmDouble Delta18oLgm, IssmDouble Delta18oTime, 
				     IssmDouble* PrecipitationsPresentday,
				     IssmDouble* TemperaturesLgm, IssmDouble* TemperaturesPresentday, 
				     IssmDouble* monthlytemperaturesout, IssmDouble* monthlyprecout){

  IssmDouble monthlytemperaturestmp[12],monthlyprectmp[12];
  IssmDouble delta18oLapseRate=-6.2*pow(10.,-3);
  IssmDouble glacialindex; // used to vary present day temperature

  glacialindex = (Delta18oTime-Delta18oPresent-delta18oLapseRate*(Delta18oSurfaceTime-Delta18oSurfacePresent))
    /(Delta18oLgm-Delta18oPresent-delta18oLapseRate*(Delta18oSurfaceLgm-Delta18oSurfacePresent)); // Tarasov 2004 paper

  for (int imonth = 0; imonth<12; imonth++){
    monthlytemperaturestmp[imonth] = glacialindex*TemperaturesLgm[imonth] + (1.-glacialindex)*TemperaturesPresentday[imonth];
    monthlyprectmp[imonth] = PrecipitationsPresentday[imonth];

    /*Assign output pointer*/
    *(monthlytemperaturesout+imonth) = monthlytemperaturestmp[imonth];
    *(monthlyprecout+imonth) = monthlyprectmp[imonth];
  }
}
