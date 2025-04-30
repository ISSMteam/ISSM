#ifndef _CONTAINER_PARAMETERS_H_
#define  _CONTAINER_PARAMETERS_H_
#include <stdio.h>

/*forward declarations */
class Param;
class DataSet;
class MarshallHandle;
template <class doublematrix> class Matrix;
template <class doubletype> class Vector;
#include "../../shared/shared.h"

#define NUMPARAMS ParametersENDEnum - ParametersSTARTEnum -1

/*!\brief Declaration of Parameters class.  
 *
 * Declaration of Parameters class.  Parameters are a static array of Parameter objects.
 */ 
class Parameters{

	private:
		Param* params[NUMPARAMS];
		int    EnumToIndex(int enum_in);

	public:

		/*constructors, destructors*/ 
		Parameters();
		~Parameters();

		/*numerics*/
		void  AddObject(Param* newparam);
		Parameters* Copy(void);
		void  DeepEcho();
		void  Echo();
		void  Delete(int enum_type);
		bool  Exist(int enum_type);
		void  Marshall(MarshallHandle* marshallhandle);

		void  FindParam(bool* pinteger,int enum_type);
		void  FindParam(int* pinteger,int enum_type);
		void  FindParam(IssmDouble* pscalar, int enum_type);
		void  FindParam(IssmDouble* pscalar, int enum_type, IssmDouble time);
		void  FindParam(IssmDouble* pscalar, int enum_type, IssmDouble time, int timestepping, IssmDouble dt);
		void  FindParam(IssmDouble* pscalar, int row, IssmDouble time, int enum_type);
		void  FindParam(IssmDouble* pscalar, int row, int column, IssmDouble time, int enum_type);
		void  FindParam(IssmDouble* pscalar, int row, IssmDouble time, int timestepping, IssmDouble dt, int enum_type);
		void  FindParam(char** pstring,int enum_type);
		void  FindParam(char*** pstringarray,int* pM,int enum_type);
		void  FindParam(int** pintarray,int* pM,int enum_type);
		void  FindParam(int** pintarray,int* pM,int* PN,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble time,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int* pN,IssmDouble starttime,IssmDouble endtime,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int enum_type);
		void  FindParam(IssmDouble** pIssmDoublearray,int* pM,int* pN,int enum_type);
		void  FindParam(IssmDouble*** parray,int* pM, int** pmdims_array,int** pndims_array,int enum_type);
		void  FindParam(Vector<IssmDouble>** pvec,int enum_type);
		void  FindParam(Matrix<IssmDouble>** pmat,int enum_type);
		void  FindParam(FILE** pfid,int enum_type);
		void  FindParam(DataSet** pdataset, int enum_type);
		void  FindParamAndMakePassive(IssmPDouble* pscalar, int enum_type);
		void  FindParamAndMakePassive(IssmPDouble** pvec,int* pM,int enum_type);
		void  FindControlParam(IssmDouble** pvec,int* pM, int param_enum, const char* data);
		void  FindControlParamAndMakePassive(IssmPDouble** pvec,int* pM, int param_enum, const char* data);
		IssmDouble FindParam(int enum_type);

		void  SetParam(bool boolean,int enum_type);
		void  SetParam(int integer,int enum_type);
		void  SetParam(IssmDouble scalar, int enum_type);
		void  SetParam(char* string,int enum_type);
		void  SetParam(char** stringarray,int M,int enum_type);
		void  SetParam(IssmDouble* IssmDoublearray,int M,int enum_type);
		void  SetParam(IssmDouble* IssmDoublearray,int M,int N,int enum_type);
		void  SetParam(IssmDouble* IssmDoublearray, int enum_type);
		void  SetParam(int* intarray,int M,int enum_type);
		void  SetParam(int* intarray,int M,int N,int enum_type);
		void  SetParam(Vector<IssmDouble>* vec,int enum_type);
		void  SetParam(Matrix<IssmDouble>* mat,int enum_type);
		void  SetParam(FILE* fid,int enum_type);
		void  SetParam(DataSet* dataset,int enum_type);
		void  SetControlFromVector(IssmDouble* array, int enum_type, int M, int N, int offset);
		void  SetGradientFromVector(IssmDouble* array, int enum_type, int M, int N, int offset);
		void  GetVectorFromControl(Vector<IssmDouble>* vector,int control_enum,int control_index,int N,const char* data,int offset);
		Param* FindParamObject(int enum_type);

};

/*Methods relating to parameters: */
char *OptionsFromAnalysis(char** ptoolkit,Parameters *parameters,int analysis_type);
void  ToolkitsOptionsFromAnalysis(Parameters* parameters,int analysis_type);

#endif //ifndef _PARAMETERS_H_
