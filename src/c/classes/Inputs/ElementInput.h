#ifndef _ELEMENTINPUT2_H_
#define _ELEMENTINPUT2_H_

/*Headers:*/
#include "./Input.h"

class ElementInput: public Input{

	protected:
		int         numberofelements_local;
		int         numberofvertices_local;
		int         interpolation;
		int         M,N;
		bool        isserved;
		IssmDouble* values;

	public:
		IssmDouble* element_values;

		/*ElementInput constructors, destructors*/ 
		ElementInput();
		~ElementInput();

		int  GetInputInterpolationType();

		/*Object virtual functions definitions:*/
		virtual Input *copy()=0;
		virtual void    DeepEcho()=0;
		virtual void    Echo()=0;
		virtual int     Id()=0;
		virtual void    Marshall(MarshallHandle* marshallhandle)=0;
		virtual int     ObjectEnum()=0;
		/*Other*/
		virtual void SetInput(int interp_in,int row,IssmDouble value_in)=0;
		virtual void SetInput(int interp_in,int numinds,int* rows,IssmDouble* values_in)=0;
		virtual void SetInput(int interp_in,int row,int numinds,IssmDouble* values_in)=0;
		virtual int  GetInterpolation()=0;
		virtual void GetInputDerivativeValue(IssmDouble* derivativevalues, IssmDouble* xyz_list, Gauss* gauss)=0;
		virtual void GetInputValue(IssmDouble* pvalue,Gauss* gauss)=0;
		virtual void Serve(int numindices,int* indices)=0;
		virtual void Serve(int row,int numindices)=0;
		virtual int  GetResultArraySize(void)=0;
		virtual int  GetResultInterpolation(void)=0;
		virtual int  GetResultNumberOfNodes(void)=0;
};
#endif  /* _ELEMENTINPUT_H */
