#ifndef _PENTAINPUT2_H_
#define _PENTAINPUT2_H_

/*Headers:*/
#include "./ElementInput.h"
#include "../Elements/PentaRef.h"

class PentaInput: public ElementInput, public PentaRef{

	private:
		int isserved_collapsed;
	public:
		/*PentaInput constructors, destructors: {{{*/
		PentaInput();
		PentaInput(int nbe_in,int nbv_in,int interp_in);
		~PentaInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*PentaInput management: {{{*/
		void SetInput(int interp_in,int row,IssmDouble value_in);
		void SetInput(int interp_in,int numinds,int* rows,IssmDouble* values_in);
		void SetInput(int interp_in,int row,int numinds,IssmDouble* values_in);
		int  GetInterpolation();
		void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss);
		void GetInputAverage(IssmDouble* pvalue);
		IssmDouble GetInputMin();
		IssmDouble GetInputMax();
		IssmDouble GetInputMaxAbs();
		PentaInput* GetPentaInput(){return this;};
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss);
		void Scale(IssmDouble scalar);
		void Pow(IssmDouble scalar);
		void AXPY(Input* xinput,IssmDouble scalar);
		void PointWiseMult(Input* xinput);
		void Serve(int numindices,int* indices);
		void Serve(int row,int numindices);
		void ServeCollapsed(int row,int state);
		void SetServeCollapsed(int);
		int  GetResultArraySize(void);
		int  GetResultInterpolation(void);
		int  GetResultNumberOfNodes(void);
		/*}}}*/
		void Reset(int interp_in);

};
#endif  /* _TRIAINPUT_H */
