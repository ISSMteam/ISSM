#ifndef _BAMGDET_H_
#define _BAMGDET_H_

#include "./include.h"

namespace bamg {

	long long inline det(const I2 &a,const I2 & b,const I2 &c){
		long long bax = b.x - a.x ,bay = b.y - a.y; 
		long long cax = c.x - a.x ,cay = c.y - a.y; 
		return  bax*cay - bay*cax;
	}

}
#endif
