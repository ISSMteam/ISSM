#ifndef _OPPOSITEANGLE_H_
#define _OPPOSITEANGLE_H_

#include "../Numerics/constants.h"

/*Return the opposite angle modulo 2 Pi*/
namespace bamg {
	inline float  OppositeAngle(float  a){return a<0 ? PI+a:a-PI;}
	inline double OppositeAngle(double a){return a<0 ? PI+a:a-PI;}
}

#endif
