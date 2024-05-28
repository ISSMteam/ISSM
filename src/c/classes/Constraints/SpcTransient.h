/*!\file SpcTransient.h
 * \brief: header file for spc object
 */

#ifndef _SPCTRANSIENT_H_
#define _SPCTRANSIENT_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
/*}}}*/

class SpcTransient: public Constraint{

	private: 
		int         id;              /* id, to track it             */
		int         nodeid;          /*node id                      */
		int         dof;             /*component                    */
		IssmDouble *values;          /*different values in time     */
		IssmDouble *times;           /*different time steps         */
		int         nsteps;          /*number of time steps         */
		int         analysis_type;
		bool        penalty;         /*Is this a penalty constraint */

	public:

		/*SpcTransient constructors, destructors:{{{*/
		SpcTransient();
		SpcTransient(int id_in,int nodeid, int dof,int nsteps, IssmDouble* times, IssmDouble* values,int analysis_type);
		~SpcTransient();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void    DeepEcho();
		void    Echo();
		int     Id(); 
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();
		/*}}}*/
		/*Constraint virtual functions definitions: {{{*/
		void   ActivatePenaltyMethod(void);
		void   ConstrainNode(Nodes* nodes,Parameters* parameters);
		void   PenaltyDofAndValue(int* dof,IssmDouble* value,Nodes* nodes,Parameters* parameters);
		void InputUpdateFromVectorDakota(IssmDouble* vector,Nodes* nodes,int name,int type){_error_("not implemented yet");};
		/*}}}*/
		/*SpcTransient management:{{{ */
		int        GetDof();
		int        GetNodeId();
		IssmDouble GetValue();
		/*}}}*/

};

#endif  /* _SPCTRANSIENT_H_ */
