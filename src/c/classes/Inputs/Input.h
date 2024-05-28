/*!\file:  Input.h
 * \brief abstract class for Input object
 */

#ifndef _INPUT_H_
#define _INPUT_H_

/*Headers:*/
#include "../../shared/shared.h"
#include "../../datastructures/Object.h"
class Gauss;
class Parameters;
class SegInput;
class TriaInput;
class PentaInput;
template <class doubletype> class Vector;

class Input: public Object{

	private:
		int enum_type;
	public:

		/*Non virtual functions*/
		int  InstanceEnum(){return this->enum_type;};
		void ChangeEnum(int newenumtype){this->enum_type=newenumtype;};

		/*Virtual functions*/
		virtual ~Input(){};
		virtual void Configure(Parameters* parameters){return;};
		virtual Input* copy()=0;
		//virtual void GetInputAllTimeAverages(IssmDouble** pvalues,IssmDouble** ptimes, int* pnumtimes){_error_("Not implemented yet");};
		virtual void  GetInputAverage(IssmDouble* pvalue){_error_("Not implemented yet");};
		virtual IssmDouble GetInputMax(void){_error_("Not implemented yet");};
		virtual IssmDouble GetInputMaxAbs(void){_error_("Not implemented yet");};
		virtual IssmDouble GetInputMin(void){_error_("Not implemented yet");};
		virtual void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss){_error_("Not implemented yet");};
		virtual void GetInputValue(IssmDouble* pvalue,Gauss* gauss){int* temp = xNew<int>(3); _error_("Not implemented yet for");};
		virtual int  GetInputInterpolationType(){_error_("Not implemented yet");};
		virtual SegInput*   GetSegInput(){ int* temp = xNew<int>(3); this->Echo(); _error_("Not implemented yet");};
		virtual TriaInput*  GetTriaInput(){ int* temp = xNew<int>(3); this->Echo(); _error_("Not implemented yet");};
		virtual PentaInput* GetPentaInput(){int* temp = xNew<int>(3); this->Echo(); _error_("Not implemented yet");};
		//virtual void GetInputUpToCurrentTimeAverages(IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){_error_("Not implemented yet");};

		virtual void   AXPY(Input* xinput,IssmDouble scalar){_error_("Not implemented yet");};
		virtual void   Shift(IssmDouble scalar){_error_("Not implemented yet");};
		virtual void   PointWiseMult(Input* xinput){_error_("Not implemented yet");};
		virtual void   Pow(IssmDouble scale_factor){_error_("Not implemented yet");};
		virtual void   Scale(IssmDouble scale_factor){_error_("Not implemented yet");};
		virtual void   AverageAndReplace(void){_error_("Not implemented yet");};

		virtual int  GetResultArraySize(void){_error_("Not implemented yet");};
		virtual int  GetResultInterpolation(void){_error_("Not implemented yet");};
		virtual int  GetResultNumberOfNodes(void){_error_("Not implemented yet");};
		//virtual void ResultToMatrix(IssmDouble* values,int ncols,int sid){_error_("not supported yet");};
		//virtual void ResultToPatch(IssmDouble* values,int nodesperelement,int sid){_error_("not supported yet");};
};
#endif
