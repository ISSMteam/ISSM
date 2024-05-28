/*! \file DatasetInput.h 
 *  \brief: header file for datasetinput object
 */

#ifndef _DATASETINPUT2_H_
#define _DATASETINPUT2_H_

/*Headers:*/
#include "./Input.h"
class TriaInput;
class PentaInput;
class TransientInput;

class DatasetInput: public Input{

	private:
		int             numids;
		Input        **inputs;
		int            *ids;
		int             numberofelements_local;
		int             numberofvertices_local;

	public:
		int GetNumIds() const {return this->numids;};
		/*DatasetInput constructors, destructors: {{{*/
		DatasetInput();
		DatasetInput(int nbe, int nbv);
		~DatasetInput();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Input* copy();
		void    Configure(Parameters* params);
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		void    SetTriaInput(int interp_in,int numinds,int* rows,IssmDouble* values_in);
		/*}}}*/
		IssmDouble GetInputMin();
		void SetTriaInput(int id,int interp_in,int numinds,int* rows,IssmDouble* values_in);
		void SetPentaInput(int id,int interp_in,int numinds,int* rows,IssmDouble* values_in);
		TransientInput* SetTransientInput(int id,IssmDouble* times,int numtimes);
		PentaInput* GetPentaInputByOffset(int i);
		TriaInput*  GetTriaInput(void);
		TriaInput*  GetTriaInputByOffset(int i);
		TransientInput* GetTransientInputByOffset(int offset);
		void GetInputValue(IssmDouble* pvalue,Gauss* gauss,int index);
};
#endif  /* _DATASETINPUT2_H */
