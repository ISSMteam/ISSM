/* file:  ComputeMungsmTemperaturePrecipitation.cpp
   Compute the temperature and precipitation at time t from 
   the data at present day and lgm.
   The interpolation is done from some factors extracted from the MUNGSM
 */

#include "./elements.h"
#include "../Numerics/numerics.h"
#include <cmath>

void ComputeMungsmTemperaturePrecipitation(IssmDouble TdiffTime, IssmDouble PfacTime,
					   IssmDouble* PrecipitationsLgm, IssmDouble* PrecipitationsPresentday,
					   IssmDouble* TemperaturesLgm, IssmDouble* TemperaturesPresentday,
					   IssmDouble* monthlytemperaturesout, IssmDouble* monthlyprecout){ 

  IssmDouble monthlytemperaturestmp[12],monthlyprectmp[12];
  IssmDouble tdiffh;  

  for (int imonth = 0; imonth<12; imonth++){
    tdiffh = TdiffTime*( TemperaturesLgm[imonth] - TemperaturesPresentday[imonth] );
    monthlytemperaturestmp[imonth] = tdiffh + TemperaturesPresentday[imonth] ;

    monthlyprectmp[imonth] =min(1.5, PrecipitationsPresentday[imonth] * pow(PrecipitationsLgm[imonth],PfacTime));   // [m/yr]

    /*Assign output pointer*/
    *(monthlytemperaturesout+imonth) = monthlytemperaturestmp[imonth];
    *(monthlyprecout+imonth) = monthlyprectmp[imonth];
  }
  // printf(" tempera %f\n",monthlytemperaturestmp[1]);
}
