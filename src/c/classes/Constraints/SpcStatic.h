/*!\file SpcStatic.h
 * \brief: header file for spc object
 */

#ifndef _SPCStatic_H_
#define _SPCStatic_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
/*}}}*/

class SpcStatic: public Constraint{

	private: 
		int        id;
		int        nodeid;
		int        dof;
		IssmDouble value;
		int        analysis_type;
		bool       penalty;

	public:

		/*SpcStatic constructors, destructors:{{{*/
		SpcStatic();
		SpcStatic(int id,int nodeid, int dof,IssmDouble value,int analysis_type);
		~SpcStatic();
		/*}}}*/
		/*Object virtual functions definitions:{{{ */
		Object* copy();
		void  DeepEcho();
		void  Echo();
		int   Id(); 
		void  Marshall(MarshallHandle* marshallhandle);
		int   ObjectEnum();
		/*}}}*/
		/*Constraint virtual functions definitions: {{{*/
		void ActivatePenaltyMethod(void);
		void ConstrainNode(Nodes* nodes,Parameters* parameters);
		void PenaltyDofAndValue(int* dof,IssmDouble* value,Nodes* nodes,Parameters* parameters);
		void InputUpdateFromVectorDakota(IssmDouble* vector,Nodes* nodes,int name,int type);
		/*}}}*/
		/*SpcStatic management:{{{ */
		int    GetDof();
		int    GetNodeId();
		IssmDouble GetValue();
		/*}}}*/
		void UpdateSpcThicknessAD(IssmDouble* vector,Nodes* nodes);

};

#endif  /* _SPCStatic_H_*/
