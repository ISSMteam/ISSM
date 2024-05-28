/*\file matlabio.h
 *\brief: I/O for ISSM in matlab mode
 */

#ifndef _MATLAB_IO_H_
#define _MATLAB_IO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

#include "../include/matlabincludes.h"
#include "../../../c/bamg/bamgobjects.h"
#include "../../../c/classes/classes.h"
#include "../../../c/toolkits/toolkits.h"
#include "../../../c/shared/shared.h"

void WriteData(mxArray** pdataref,IssmDenseMat<double>* matrix);
void WriteData(mxArray** pdataref,double* matrix, int M,int N);
void WriteData(mxArray** pdataref,int*    matrix, int M,int N);
void WriteData(mxArray** pdataref,IssmSeqVec<double>* vector);
void WriteData(mxArray** pdataref,double* vector, int M);
void WriteData(mxArray** pdataref,short* vector, int M);
void WriteData(mxArray** pdataref,int* vector, int M);
void WriteData(mxArray** pdataref,int integer);
void WriteData(mxArray** pdataref,bool boolean);
void WriteData(mxArray** pdataref,double scalar);
void WriteData(mxArray** pdataref,const char* string);
void WriteData(mxArray** pdataref);
void WriteData(mxArray** pdataref,BamgGeom* bamggeom);
void WriteData(mxArray** pdataref,BamgMesh* bamgmesh);
void WriteData(mxArray** pdataref,RiftStruct* riftstruct);
void WriteData(mxArray** pdataref,Contours* contours);

void FetchData(double** pmatrix,int* pM,int *pN,const mxArray* dataref);
void FetchData(int** pmatrix,int* pM,int *pN,const mxArray* dataref);
void FetchData(bool** pmatrix,int* pM,int *pN,const mxArray* dataref);
void FetchData(int** pvector,int* pM,const mxArray* dataref);
void FetchData(float** pvector,int* pM,const mxArray* dataref);
void FetchData(double** pvector,int* pM,const mxArray* dataref);
void FetchData(bool** pvector,int* pM,const mxArray* dataref);
void FetchData(char** pstring,const mxArray* dataref);
void FetchData(double* pscalar,const mxArray* dataref);
void FetchData(int* pinteger,const mxArray* dataref);
void FetchData(bool* pbool,const mxArray* dataref);
void FetchData(BamgGeom** bamggeom,const mxArray* dataref);
void FetchData(BamgMesh** bamgmesh,const mxArray* dataref);
void FetchData(BamgOpts** bamgopts,const mxArray* dataref);
void FetchData(Options** poptions,int istart, int nrhs,const mxArray** pdataref);
void FetchData(Contours** pcontours,const mxArray* dataref);
void FetchChacoData(int* pnvtxs,int** padjacency,int** pstart,float** pewgts,const mxArray* A_IN, const mxArray* EWGTS_IN);

Option* OptionParse(char* name, const mxArray* prhs[]);
GenericOption<double>*  OptionDoubleParse(char* name, const mxArray* prhs[]);
GenericOption<double*>* OptionDoubleArrayParse(char* name, const mxArray* prhs[]);
GenericOption<char*>*   OptionCharParse(char* name, const mxArray* prhs[]);

mxArray* mxGetAssignedField(const mxArray* pmxa_array,int number, const char* field);
void SetStructureField(mxArray* dataref,const char* fieldname,int fieldrows,int fieldcols,double* fieldpointer);
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int fieldrows,int fieldcols,double* fieldpointer);
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int fieldrows,int fieldcols,int*    fieldpointer);
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,int field);
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,double field);
void SetStructureFieldi(mxArray* dataref,int i,const char* fieldname,const char* string);
int CheckNumMatlabArguments(int nlhs,int NLHS, int nrhs,int NRHS, const char* THISFUNCTION, void (*function)( void ));

/*Matlab to double* routines: */
int MatlabMatrixToDoubleMatrix(double** pmatrix,int* pmatrix_rows,int* pmatrix_cols,const mxArray* mxmatrix);

/*Print*/
void ApiPrintf(const char* string);
#endif	/* _IO_H_ */
