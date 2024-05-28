/*\file pythonio.h
 *\brief: I/O for ISSM in python mode
 */

#ifndef _PYTHON_IO_H_
#define _PYTHON_IO_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif 

#include "../include/pythonincludes.h"
#include "../../c/bamg/bamgobjects.h"
#include "../../c/classes/classes.h"
#include "../../c/shared/shared.h"

void WriteData(PyObject* py_tuple,int index, double* matrix, int M,int N);
void WriteData(PyObject* py_tuple,int index, int* matrix, int M,int N);
void WriteData(PyObject* py_tuple,int index, bool* matrix, int M,int N);
void WriteData(PyObject* py_tuple,int index, int integer);
void WriteData(PyObject* py_tuple,int index, double* vector, int M);
void WriteData(PyObject* py_tuple,int index, short* vector, int M);
void WriteData(PyObject* py_tuple,int index, int* vector, int M);
void WriteData(PyObject* py_tuple,int index, char* string);
void WriteData(PyObject* py_tuple,int index);
void WriteData(PyObject* py_tuple,int index, IssmDenseMat<double>* matrix);
void WriteData(PyObject* py_tuple,int index, IssmSeqVec<double>* vector);
void WriteData(PyObject* py_tuple,int index, BamgGeom* bamggeom);
void WriteData(PyObject* py_tuple,int index, BamgMesh* bamgmesh);
void WriteData(PyObject* py_tuple,int index, RiftStruct* riftstruct);

void FetchData(double** pmatrix,int* pM,int *pN,PyObject* py_array);
void FetchData(int** pmatrix,int* pM,int *pN,PyObject* py_matrix);
void FetchData(bool** pmatrix,int* pM,int *pN,PyObject* py_matrix);
void FetchData(double** pvector,int* pM,PyObject* py_ref);
void FetchData(float** pvector,int* pM,PyObject* dataref);
void FetchData(int** pvector,int* pM,PyObject* py_ref);
void FetchData(bool** pvector,int* pM,PyObject* py_ref);
void FetchData(char** pstring,PyObject* py_unicode);
void FetchData(double* pscalar,PyObject* py_float);
void FetchData(short* pscalar,PyObject* py_float);
void FetchData(int* pscalar,PyObject* py_long);
void FetchData(bool* pbool,PyObject* py_boolean);
void FetchData(BamgGeom** bamggeom,PyObject* py_dict);
void FetchData(BamgMesh** bamgmesh,PyObject* py_dict);
void FetchData(BamgOpts** bamgopts,PyObject* py_dict);
void FetchData(Options** poptions,int istart, int nrhs,PyObject* py_tuple);
void FetchData(Contours** pcontours,PyObject* py_list);
void FetchChacoData(int* pnvtxs,int** padjacency,int** pstart,float** pewgts,PyObject* A_IN,PyObject* EWGTS_IN);

void pyGetJc(double* a, int nvtxs, int* Jc);
void pyGetIr(double* a, int nvtxs, int nzmax, int* Ir);

int CheckNumPythonArguments(PyObject* inputs,int NRHS, void (*function)( void ));

/*Utils*/
PyObject* PyArrayFromCopiedData(int dims[2],double* data);
PyObject* PyArrayFromCopiedData(int dimi,int dimj,double* data);
PyObject* PyArrayFromCopiedData(int dimi,int dimj,int* data);
PyObject* PyArrayFromCopiedData(int dimi,int dimj,bool* data);

/*Print*/
void ApiPrintf(const char* string);

#endif	/* _IO_H_ */
