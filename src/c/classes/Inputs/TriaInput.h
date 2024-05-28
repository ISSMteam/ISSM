#ifndef _TRIAINPUT2_H_
#define _TRIAINPUT2_H_

/*Headers:*/
#include "./ElementInput.h"
#include "../Elements/TriaRef.h"

class TriaInput: public ElementInput, public TriaRef{

	private:
		int isserved_collapsed;
		int collapsed_ids[2];
	public:
		/*TriaInput constructors, destructors: {{{*/
		TriaInput();
		TriaInput(int nbe_in,int nbv_in,int interp_in);
		~TriaInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*TriaInput management: {{{*/
		void SetInput(int interp_in,int row,IssmDouble value_in);
		void SetInput(int interp_in,int numinds,int* rows,IssmDouble* values_in);
		void SetInput(int interp_in,int row,int numinds,IssmDouble* values_in);
		int  GetInterpolation();
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss);
		void GetInputAverage(IssmDouble* pvalue);
		IssmDouble GetInputMin();
		IssmDouble GetInputMax();
		IssmDouble GetInputMaxAbs();
		TriaInput* GetTriaInput(){return this;};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void Scale(IssmDouble scalar);
		void Pow(IssmDouble scalar);
		void AXPY(Input* xinput,IssmDouble scalar);
		void Shift(IssmDouble scalar);
		void AverageAndReplace(void);
		void PointWiseMult(Input* xinput);
		void Serve(int numindices,int* indices);
		void Serve(int row,int numindices);
		void ServeCollapsed(int row,int id0,int in1);
		void SetServeCollapsed(bool);
		int  GetResultArraySize(void);
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		/*}}}*/
		void Reset(int interp_in);

};
#endif  /* _TRIAINPUT_H */
