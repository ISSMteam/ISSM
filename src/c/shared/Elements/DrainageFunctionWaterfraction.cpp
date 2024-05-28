/*!\file DrainageFunctionWaterfraction.cpp
 * \brief: drain excess water fraction
 */

#include <math.h>
#include "../Numerics/types.h"
#include "../Exceptions/exceptions.h"

IssmDouble DrainageFunctionWaterfraction(IssmDouble waterfraction, IssmDouble dt=0.){
	/* DrainageFunctionWaterfraction returns how much of the waterfraction is drained per year */
	_assert_(waterfraction>=0.);
	_assert_(dt>=0.);

	IssmDouble w0=0.01, w1=0.02, w2=0.03;
	IssmDouble yts=365.*24.*60.*60.;
	IssmDouble Dret, D0=0., D1=0.005/yts, D2=0.05/yts;

	/*get drainage function value*/
	if((w0==w1)||(w1==w2)||(w0==w2))
		_error_("Error: equal ordinates in DrainageFunctionWaterfraction -> division by zero. Abort");

	if(waterfraction<=w0)
		Dret=D0;
	else if((waterfraction>w0) && (waterfraction<=w1))
		Dret=(D1-D0)/(w1-w0)*(waterfraction-w0)+D0;
	else if((waterfraction>w1) && (waterfraction<=w2))
		Dret=(D2-D1)/(w2-w1)*(waterfraction-w1)+D1;
	else 
		Dret=D2;

	/*drain only up to w0*/
	if(dt==0.){
		if(waterfraction>w0) 
			return waterfraction-w0;
		else
			return Dret;
	}
	else{
		if((waterfraction>w0) && (waterfraction-dt*Dret<w0))
			return (waterfraction-w0)/dt;
		else
			return Dret;
	}
}
