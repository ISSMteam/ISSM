/*! \file TransientInput.h
 *  \brief: header file for transientinput object
 */

#ifndef _TRANSIENTINPUT2_H_
#define _TRANSIENTINPUT2_H_

/*Headers:*/
#include "./Input.h"
class Gauss;
class Parameters;

class TransientInput: public Input{

	private:
		int     numberofelements_local;
		int     numberofvertices_local;

	public:
		int          enum_type;
		int          numtimesteps;
		Input      **inputs;
		IssmDouble  *timesteps;
		Parameters  *parameters;      //to find current time.

		IssmDouble   current_step;
		Input       *current_input;

		/*TransientInput constructors, destructors: {{{*/
		TransientInput();
		TransientInput(int in_enum_type,int nbe,int nbv,IssmDouble* times,int N);
		~TransientInput();
		void AddTriaTimeInput(IssmDouble time,int numindices,int* indices,IssmDouble* values_in,int interp_in);
		void AddPentaTimeInput(IssmDouble time,int numindices,int* indices,IssmDouble* values_in,int interp_in);
		void AddTriaTimeInput(int step,int numindices,int* indices,IssmDouble* values_in,int interp_in);
		void AddPentaTimeInput(int step,int numindices,int* indices,IssmDouble* values_in,int interp_in);
		/*}}}*/
		/*Object virtual functions definitions:{{{*/
		Input*  copy();
		void    Configure(Parameters* params);
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*TransientInput management:*/
		void         GetAllTimes(IssmDouble** ptimesteps,int* pnumtimesteps);
		TriaInput*  GetTriaInput();
		TriaInput*  GetTriaInput(IssmDouble time);
		TriaInput*  GetTriaInput(IssmDouble start_time,IssmDouble end_time,int averaging_method);
		TriaInput*  GetTriaInput(int offset);
		PentaInput* GetPentaInput();
		PentaInput* GetPentaInput(IssmDouble time);
		PentaInput* GetPentaInput(int offset);
		PentaInput* GetPentaInput(IssmDouble start_time,IssmDouble end_time,int averaging_method);
		Input*      GetTimeInput(IssmDouble time){_error_("This should not happen!");};
		IssmDouble   GetTimeByOffset(int offset);
		int          GetTimeInputOffset(IssmDouble time);
		void         SetCurrentTimeInput(IssmDouble time);
		void         SetAverageAsCurrentTimeInput(IssmDouble start_time,IssmDouble end_time,int averaging_method);
		/*numerics:*/

};
#endif  /* _TRANSIENTINPUT_H */
