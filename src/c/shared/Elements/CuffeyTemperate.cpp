/* \file CuffeyTemperate.cpp
 * \brief figure out B of ice for a certain temperature and waterfraction
 *	  INPUT function B=Cuffey(temperature, waterfraction)
 *    where rigidity (in s^(1/3)Pa) is the flow law parameter in the flow law sigma=B*e(1/3) (Cuffey, p75). 
 */

#include <math.h>
#include "./elements.h"
#include "../Numerics/numerics.h"

IssmDouble CuffeyTemperate(IssmDouble temperature, IssmDouble waterfraction, IssmDouble stressexp){

	return Cuffey(temperature)*pow(1+181.25*max(0., min(0.01, waterfraction)), -1./stressexp); 

}
