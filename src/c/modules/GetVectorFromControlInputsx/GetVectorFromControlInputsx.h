/*!\file:  GetVectorFromControlInputsx.h
 */ 

#ifndef _GETVECTORFROMCONTROLINPUTSXX_H
#define _GETVECTORFROMCONTROLINPUTSXX_H

#include "../../classes/classes.h"

/* local prototypes: */
void	GetVectorFromControlInputsx( Vector<IssmDouble>** pvector, Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads, Materials* materials,  Parameters* parameters,const char* data="value");
void	GetVectorFromControlInputsx( IssmDouble** pvector,int* pN, Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads, Materials* materials,  Parameters* parameters,const char* data="value");

void	GetPassiveVectorFromControlInputsx(double** pvector,int* pN, Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads, Materials* materials,  Parameters* parameters,const char* data="value");

#endif  /* _GETVECTORFROMCONTROLINPUTSXX_H */
