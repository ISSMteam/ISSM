/*!\file:  OutputDefinitionsResponsexx.h
*/ 

#ifndef _OUTPUTDEFINITIONSRESPONSEXX_H
#define _OUTPUTDEFINITIONSRESPONSEXX_H

#include "../../classes/classes.h"

/* local prototypes: */
int OutputDefinitionsResponsex(IssmDouble* presponse, FemModel* femmodel,const char* output_string);
int OutputDefinitionsResponsex(IssmDouble* presponse, FemModel* femmodel,int output_enum);

#endif  /* _OUTPUTDEFINITIONSRESPONSEXX_H */
