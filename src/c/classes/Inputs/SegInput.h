#ifndef _SEGINPUT2_H_
#define _SEGINPUT2_H_

/*Headers:*/
#include "./ElementInput.h"
#include "../Elements/SegRef.h"

class SegInput: public ElementInput, public SegRef{

	public:
		/*SegInput constructors, destructors: {{{*/
		SegInput();
		SegInput(int nbe_in,int nbv_in,int interp_in);
		~SegInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*SegInput management: {{{*/
		void SetInput(int interp_in,int row,IssmDouble value_in);
		void SetInput(int interp_in,int numinds,int* rows,IssmDouble* values_in);
		void SetInput(int interp_in,int row,int numinds,IssmDouble* values_in);
		int  GetInterpolation();
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss);
		void GetInputAverage(IssmDouble* pvalue);
		IssmDouble GetInputMin();
		IssmDouble GetInputMax();
		IssmDouble GetInputMaxAbs();
		SegInput* GetSegInput(){return this;};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void Scale(IssmDouble scalar);
		void Pow(IssmDouble scalar);
		void AXPY(Input* xinput,IssmDouble scalar);
		void PointWiseMult(Input* xinput);
		void Serve(int numindices,int* indices);
		void Serve(int row,int numindices);
		int  GetResultArraySize(void);
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		/*}}}*/
		void Reset(int interp_in);

};
#endif  /* _SEGINPUT_H */
