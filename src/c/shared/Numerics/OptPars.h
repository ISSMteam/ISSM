/*!\file:  OptPars.h
 * \brief place holder for optimization parameters
 */ 

#ifndef _OPTPARS_H_
#define _OPTPARS_H_

#include "./types.h"

struct OptPars{

	IssmDouble  xmin;
	IssmDouble  xmax;
	IssmDouble *cm_jump;
	int* maxiter;
	int  nsteps;
	int  nsize;

};

#endif
