/*!\file:  Trianglex.h
 * \brief header file for Trianglex module
 */ 

#ifndef _TRIANGLEX_H_
#define _TRIANGLEX_H_

#include <string.h>
#include "../../classes/classes.h"

/* local prototypes: */
void Trianglex(int** pindex,IssmPDouble** px,IssmPDouble** py,int** psegments,int** psegmentmarkerlist,int* pnels,int* pnods, int* pnseg,Contours* domain,Contours* rifts,double area);
#endif  /* _TRIANGLEX_H */
