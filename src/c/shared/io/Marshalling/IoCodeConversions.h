#ifndef _IOCODECONVERSION_H_
#define _IOCODECONVERSION_H_

void FieldAndEnumFromCode(int* out_enum,char** pfield,const char* string_in);

int IoCodeToEnumSMB(int enum_in);
int IoCodeToEnumBasal(int enum_in);
int IoCodeToEnumCalving(int enum_in);
int IoCodeToEnumFrontalforcings(int enum_in);
int IoCodeToEnumHydrology(int enum_in);
int IoCodeToEnumMaterials(int enum_in);
int IoCodeToEnumNature(int enum_in);
int IoCodeToEnumTimestepping(int enum_in);
int IoCodeToEnumAmr(int enum_in);
int IoCodeToEnumGrd(int enum_in);

int IoCodeToEnumVertexEquation(int enum_in);
int IoCodeToEnumElementEquation(int enum_in);

int IoRiftfillToEnum(int enum_in);

#endif	
