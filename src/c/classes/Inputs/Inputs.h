#ifndef _CONTAINER_INPUTS2_H_
#define _CONTAINER_INPUTS2_H_

/*forward declarations */
class Input;
class SegInput;
class TriaInput;
class PentaInput;
class TransientInput;
class ElementInput;
class DatasetInput;
class ArrayInput;
class IntArrayInput;
class ControlInput;
class Parameters;
#include "../../shared/shared.h"

#define NUMINPUTS InputsENDEnum - InputsSTARTEnum -1

/*!\brief Declaration of Inputs class.
 *
 * Declaration of Inputs class.  Inputs are a static array of Input objects.
 */
class Inputs{

	private:
		/*Private fields*/
		Input* inputs[NUMINPUTS];
		int     numberofelements_local;
		int     numberofvertices_local;

		/*Private functions*/
		int     EnumToIndex(int enum_in);

	public:

		/*constructors, destructors*/
		Inputs();
		Inputs(int nbe,int nbv);
		~Inputs();

		/*numerics*/
		void     AddInput(Input* in_input);
		void     ChangeEnum(int enumtype,int new_enumtype);
		void     Configure(Parameters* parameters);
		Inputs* Copy(void);
		int      DeleteInput(int enum_type);
		void     DuplicateInput(int original_enum,int new_enum);
		void     ZAXPY(IssmDouble alpha, int xenum, int yenum, int zenum);
		void     AXPY(IssmDouble alpha, int xenum, int yenum);
		void     Shift(int inputenum, IssmDouble alpha);
		void     AverageAndReplace(int inputenum);
		void     DeepEcho(void);
		void     DeepEcho(int enum_in);
		void     Echo(void);
		void     Echo(int enum_in);
		bool     Exist(int enum_type);
		void     GetInputsInterpolations(int* pnuminputs,int** pinterpolations,int** penum);
		void             GetArray(int enum_in,int row,IssmDouble** pvalues,int* pN);
		void             GetIntArray(int enum_in,int row,int** pvalues,int* pN);
		void             GetArrayPtr(int enum_in,int row,IssmDouble** pvalues,int* pN);
		void             GetIntArrayPtr(int enum_in,int row,int** pvalues,int* pN);
		SegInput*       GetSegInput(int enum_type);
		TriaInput*      GetTriaInput(int enum_type);
		TriaInput*      GetTriaInput(int enum_type,IssmDouble time);
		TriaInput*      GetTriaInput(int enum_in,IssmDouble start_time,IssmDouble end_time,int averaging_method);
		PentaInput*     GetPentaInput(int enum_type);
		PentaInput*     GetPentaInput(int enum_type,IssmDouble time);
		PentaInput*     GetPentaInput(int enum_in,IssmDouble start_time,IssmDouble end_time,int averaging_method);
		TransientInput* GetTransientInput(int enum_type);
		ElementInput*   GetControlInputData(int enum_type,const char* data);
		DatasetInput*   GetDatasetInput(int enum_type);
		ControlInput*   GetControlInput(int enum_type);
		void  Marshall(MarshallHandle* marshallhandle);
		int   GetInputObjectEnum(int enum_type);
		void  GetInputValue(bool* pvalue,int enum_in,int index);
		void  GetInputValue(int*  pvalue,int enum_in,int index);
		void  GetInputValue(IssmDouble*  pvalue,int enum_in,int index);
		bool  IsFileInputUpdate(IssmDouble time);
		void  ResultInterpolation(int* pinterpolation,int*nodesperelement,int* parray_size, int output_enum);
		void  SetInput(int enum_in,int index,bool value);
		void  SetInput(int enum_in,int index,int value);
		void  SetDoubleInput(int enum_in,int index,IssmDouble value);
		void  SetTransientInput(int enum_in,IssmDouble* times,int numtimes);
		void  SetControlInput(int enum_in,int layout,int interpolation,int control_id);
		void  SetTransientControlInput(int enum_in,int control_id,IssmDouble* times,int numtimes);
		TransientInput* SetDatasetTransientInput(int enum_in,int id,IssmDouble* times,int numtimes);
		void  SetArrayInput(int enum_in,int row,IssmDouble* layers,int numlayers);
		void  SetIntArrayInput(int enum_in,int row,int* layers,int numlayers);
		void  SetTriaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetTriaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values,int n);
		void  SetTriaDatasetInput(int enum_in,int id,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetTriaInput(int enum_in,int interpolation,int row,IssmDouble values);
		void  SetTriaInput(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetTriaInput(int enum_in,int interpolation,int row,int numindices,IssmDouble* values);
		void  SetPentaControlInputGradient(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetPentaDatasetInput(int enum_in,int id,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetPentaInput(int enum_in,int interpolation,int row,IssmDouble values);
		void  SetPentaInput(int enum_in,int interpolation,int numindices,int* indices,IssmDouble* values);
		void  SetPentaInput(int enum_in,int interpolation,int row,int numindices,IssmDouble* values);
};

#endif //ifndef _INPUTS_H_
