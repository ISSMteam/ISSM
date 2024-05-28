/*!\file:  cross.cpp
 * \brief cross product for 2 vectors
 */ 

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./types.h"

void cross(IssmDouble* result,IssmDouble* vector1,IssmDouble* vector2){

	/*result,vector1 and vector2 are all assumed to be of size 3: */

	result[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
	result[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];
	result[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];

}
