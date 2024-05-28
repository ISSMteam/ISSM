/* file:  ComputeTemperaturePrecipitationfrom018.cpp
 Scale present day monthly precipitation and temperature fields
 along the NGRIP oxygen isotope record.
 */

#include "./elements.h"
#include "../Numerics/numerics.h"
#include <cmath>

void ComputeD18OTemperaturePrecipitationFromPD(IssmDouble d018,IssmDouble dpermil,bool isTemperatureScaled,
			bool isPrecipScaled, IssmDouble f, IssmDouble* PrecipitationPresentday,IssmDouble* TemperaturePresentday,
			IssmDouble* PrecipitationReconstructed,IssmDouble* TemperatureReconstructed, IssmDouble* monthlytemperaturesout, 
			IssmDouble* monthlyprecout){

  IssmDouble monthlytemperaturestmp[12],monthlyprectmp[12];
  IssmDouble deltaTemp;

  /* Constants */
  // dpermil = 2.4;/*degrees C per mil*/

  /*Create Delta Temp to be applied to monthly temps and used in precip scaling*/
  deltaTemp = dpermil * (d018+34.83);   

  for(int imonth = 0; imonth<12; imonth++){

	 if(isTemperatureScaled)monthlytemperaturestmp[imonth] = TemperaturePresentday[imonth] + deltaTemp;
	 else{
		 monthlytemperaturestmp[imonth] = TemperatureReconstructed[imonth];
		 deltaTemp=TemperatureReconstructed[imonth]-TemperaturePresentday[imonth];
	 }

	 if (isPrecipScaled)monthlyprectmp[imonth] = PrecipitationPresentday[imonth]*exp((f/dpermil)*deltaTemp);
	 else monthlyprectmp[imonth] = PrecipitationReconstructed[imonth];

    /*Assign output pointer*/
    *(monthlytemperaturesout+imonth) = monthlytemperaturestmp[imonth];
    *(monthlyprecout+imonth) = monthlyprectmp[imonth];
  }
}
