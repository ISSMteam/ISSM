
#ifndef _ISSM_ENUM_
#define _ISSM_ENUM_

#include "./EnumDefinitions.h"
const char* EnumToStringx(int enum_in);
void        EnumToStringx(char** string,int enum_in);
int         StringToEnumx(const char* string_in,bool notfounderror=true);
bool        IsInputEnum(int enum_in);
bool        IsParamEnum(int enum_in);

#endif
