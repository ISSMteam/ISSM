/*! \file TransientFileInput.h
 *  \brief: header file for transientinput object
 */

#ifndef _TRANSIENTFILEINPUT_H_
#define _TRANSIENTFILEINPUT_H_

/*Headers:*/
#include "./Input.h"
class Gauss;
class Parameters;

class TransientFileInput: public Input{

	private:
		int         numberofelements_local;
		int         numberofvertices_local;
		char       *filename;
		IssmPDouble loading_period;

	public:
		int          enum_type;
		int          numtimesteps;
		Input      **inputs;
		IssmDouble  *timesteps;
		Parameters  *parameters;      //to find current time.

		IssmDouble   current_step;
		Input       *current_input;

		/*TransientFileInput constructors, destructors*/
		TransientFileInput();
		TransientFileInput(int in_enum_type,int nbe,int nbv,char* in_filename,IssmPDouble period_in);
		~TransientFileInput();

		/*Object virtual functions definitions:*/
		Input*  copy();
		void    Configure(Parameters* params);
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();

		/*TransientFileInput management:*/
		bool IsFileInputUpdate(IssmDouble time);
};
#endif  /* _TRANSIENTFILEINPUT_H */
