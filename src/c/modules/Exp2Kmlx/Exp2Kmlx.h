/*!\file:  Exp2Kmlx.h
 * \brief header file for exp to kml conversion routines.
 */ 

#ifndef _EXP2KMLX_H
#define _EXP2KMLX_H

#include <float.h>    /*  DBL_MAX  */
#include "../../classes/classes.h"

/* local prototypes: */
int Exp2Kmlx(char* filexp,char* filkml,int sgn,bool holes);
int Exp2Kmlx(char* filexp,char* filkml,int sgn,double cm,double sp,bool holes);

#endif  /* _EXP2KMLX_H */
