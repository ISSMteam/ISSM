/*!\file:  GaussPoints.h
 * \brief
 */ 

#include "../../shared/Numerics/types.h"
#ifndef _GAUSSPOINTS_H
#define _GAUSSPOINTS_H

#define MAX_LINE_GAUS_PTS    7
void GaussLegendreLinear(IssmPDouble** pxgaus, IssmPDouble** pxwgt, int ngaus);
#define MAX_TRIA_SYM_ORD    20
void GaussLegendreTria(int* pngaus, IssmDouble** pl1, IssmDouble** pl2, IssmDouble** pl3, IssmDouble** pwgt, int iord);
#define MAX_TETRA_SYM_ORD    6
void GaussLegendreTetra(int* pngaus, IssmDouble** pl1, IssmDouble** pl2, IssmDouble** pl3, IssmDouble** pl4, IssmDouble** pwgt, int iord);
#define MAX_GAUS_ITER   30
void GaussRecur(IssmPDouble* zero, IssmPDouble* weight, int n, IssmPDouble* alpha, IssmPDouble* beta);

#endif
