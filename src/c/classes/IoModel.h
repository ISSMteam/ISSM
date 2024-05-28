/* \file IoModel.h
 * \brief  Header file defining the IoModel structure that will help in processing the input data coming 
 * into ISSM, from Matlab, or through a binary file opened for reading.
 * \sa IoModel.cpp
 */

#ifndef _IOMODEL_H
#define _IOMODEL_H

#include "../shared/Enum/Enum.h"
#include <vector>

class Parameters;
class Elements;
class Inputs;
class Param;
class Options;

class IoConstant { /*holds single IssmDouble, int, bool and char from input*/
	public:
		Param* constant; 
		bool   isindependent;
		char*  name;

		~IoConstant();
		IoConstant();
		IoConstant(bool value,const char* name_in);
		IoConstant(int value,const char* name_in);
		IoConstant(IssmDouble value,const char* name_in);
		IoConstant(char* value,const char* name_in);
		IoConstant(char** value,int numstrings,const char* name_in);
};

class IoData { /*holds temporary data (array), memory intensive*/
	public:
		int         code;
		IssmDouble* data;
		bool        isindependent;
		int         layout;
		int         M,N;
		char*       name;

		~IoData();
		IoData();
		IoData(IssmDouble* matrix,int code,int layout_in,int M,int N,const char* name_in);
};

class IoModel {

	private: 
		std::vector<IoConstant*> constants; //this dataset holds all IssmDouble, int, bool and char from input
		std::vector<IoData*>     data;      //this dataset holds temporary data, memory intensive

		/*for AD mode: to keep track of our independent variables we fetch:*/
		//bool    *independents;
		//DataSet *independent_objects;

	public:
		/*pointer to input file*/
		FILE *fid;

		/*Solution*/
		int   solution_enum;

		/*Partitioning*/
		bool *my_elements;
		bool *my_faces;
		bool *my_vfaces;
		bool *my_edges;
		bool *my_vedges;
		bool *my_hedges;
		bool *my_vertices;
		int  *my_vertices_lids;
		int  *epart;

		/*Mesh properties and connectivity tables*/
		int  domaindim;
		int  domaintype;
		int *elements;
		int *edges;
		int *verticaledges;
		int *horizontaledges;
		int *elementtoedgeconnectivity;
		int *elementtoverticaledgeconnectivity;
		int *elementtohorizontaledgeconnectivity;
		int *elementtofaceconnectivity;
		int *elementtoverticalfaceconnectivity;
		int *faces;
		int *verticalfaces;
		int  facescols;
		int  meshelementtype;
		int *numbernodetoelementconnectivity;
		int  numberofedges;
		int  numberofverticaledges;
		int  numberofhorizontaledges;
		int  numberofelements;
		int  numberoffaces;
		int  numberofverticalfaces;
		int  numberofvertices;
		int *singlenodetoelementconnectivity;

		/*Methods*/
		~IoModel();
		IoModel();
		IoModel(FILE* iomodel_handle,int solution_enum_in,bool trace,IssmPDouble* X);

		/*NEW*/
		void        AddConstant(IoConstant* constant_in);
		void        AddConstantIndependent(IoConstant* constant_in);
		void        AddData(IoData* data_in);
		void        AddDataIndependent(IoData* data_in);
		void        FetchIndependentConstant(int* pXcount,IssmPDouble* X,const char* name);
		void        FetchIndependentData(int* pXcount,IssmPDouble* X,const char* name);
		void        FillIndependents(IssmDouble* xp);
		void        FindConstant(bool* pvalue,const char* constant_name);
		void        FindConstant(int* pvalue,const char* constant_name);
		void        FindConstant(IssmDouble* pvalue,const char* constant_name);
		void        FindConstant(char **pvalue,const char* constant_name);
		void        FindConstant(char ***pvalue,int* psize,const char* constant_name);
		int         NumIndependents();

		/*Input/Output*/
		void        CheckFile(void);
		Param      *CopyConstantObject(const char* constant_name,int param_enum);
		void        ConstantToInput(Inputs* inputs,Elements* elements,IssmDouble value, int vector_enum,int type);
		IssmDouble *Data(const char* data_name);
		void        DeclareIndependents(bool trace,IssmPDouble* X);
		void        DeleteData(int num,...);
		void        DeleteData(IssmDouble* vector,const char* data_name);
		void        DeleteData(char*** pstringarray, int numstrings,const char* data_name);
		void        FetchConstants(void);
		void        FetchData(bool* pboolean,const char* data_name);
		void        FetchData(int* pinteger,const char* data_name);
		void        FetchData(IssmDouble* pscalar,const char* data_name);
		void        FetchData(IssmDouble** pscalar, const char* data_name);	
		void        FetchData(char** pstring,const char* data_name);
		void        FetchData(char*** pstrings,int* pnumstrings,const char* data_name);
		void        FetchData(int** pmatrix,int* pM,int* pN,const char* data_name);
		void        FetchData(bool**  pboolmatrix,int* pM,int* pN,const char* data_name);
		void        FetchData(IssmDouble**  pscalarmatrix,int* pM,int* pN,const char* data_name);
#if _HAVE_AD_  && !defined(_WRAPPERS_)
		void        FetchData(IssmPDouble**  pscalarmatrix,int* pM,int* pN,const char* data_name);
		void        FetchData(IssmPDouble** pscalar,const char* data_name);
#endif
		void        FetchData(IssmDouble*** pmatrixarray,int** pmdims,int** pndims, int* pnumrecords,const char* data_name);
		void        FetchData(IssmDouble** pmatrix,int* pM,int* pN, int layer_number,const char* data_name);
		void        FetchData(Options *options,const char* data_name);
		void        FetchData(int num,...);
		void        FetchDataToInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum);
		void        FetchDataToInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum,IssmDouble default_value);
		void        FetchDataToDatasetInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum);
		void        FetchIndependent(const char* dependent_name);
		void        FetchMultipleData(char***   pstringarray,int* pnumstrings,const char* data_name);
		void        FetchMultipleData(IssmDouble*** pmatrixarray,int** pmdims,int** pndims, int* pnumrecords,const char* data_name);
		void        FetchMultipleData(int*** pmatrices,int** pmdims,int** pndims, int* pnumrecords,const char* data_name);
		void        FetchMultipleData(int** pvector, int* pnum_instances,const char* data_name);
		void        FetchMultipleData(IssmDouble** pvector, int* pM,const char* data_name);
		fpos_t*     SetFilePointersToData(int** pcodes,int** pvector_types, int* pnum_instances, const char* data_name);
		FILE*       SetFilePointerToData(int* pcode,int* pvector_type, const char* data_name);
		void        StartTrace(bool trace);
};

#endif  /* _IOMODEL_H */
