/*! \file ControlInput.h 
 *  \brief: header file for triavertexinput object
 */

#ifndef _CONTROLINPUT2_H_
#define _CONTROLINPUT2_H_

/*Headers:*/
#include "./Input.h"
class Gauss;
class ElementInput;
class TransientInput;

class ControlInput: public Input{

	public:
		int    control_id;
		int    enum_type;
		int    layout_enum;
		Input *maxvalues;
		Input *minvalues;
		Input *values;

		/*ControlInput constructors, destructors:*/
		ControlInput();
		ControlInput(int nbe, int nbv,int input_layout_enum,int interp,int id);
		ControlInput(int enum_in,int nbe, int nbv,int id,IssmDouble* times, int numtimes);
		~ControlInput();

		/*Object virtual functions definitions:*/
		Input* copy();
		void   Configure(Parameters* params);
		void   DeepEcho();
		void   Echo();
		void   Marshall(MarshallHandle* marshallhandle);
		int    ObjectEnum(){return ControlInputEnum;}

		void SetInput(Input* in_input){_error_("not impelemented");};
		void SetInput(Input* in_input,int timeoffset){_error_("not impelemented");};
		ElementInput*   GetInput(const char* data);
		TransientInput* GetTransientInput(const char* data);
		void        SetControl(int interp,int numindices,int* indices,IssmDouble* values_in,IssmDouble* values_min,IssmDouble* values_max);
		TriaInput*  GetTriaInput();
		PentaInput* GetPentaInput();
		void AverageAndReplace(void);
};
#endif  /* _CONTROLINPUT_H */
