/*!\file SpcDynamic.h
 * \brief: header file for spc object
 */

#ifndef _SPCDynamic_H_
#define _SPCDynamic_H_

/*Headers:*/
/*{{{*/
#include "../../datastructures/datastructures.h"
/*}}}*/

class SpcDynamic: public Constraint{

	private: 
		int        id;              /*! id, to track it */
		int        nodeid;          /*!node id          */
		int        dof;             /*!component        */
		IssmDouble value;           /*value             */
		bool       isset;
		int        analysis_type;
		bool       penalty;         /*Is this a penalty constraint */

	public:

		/*SpcDynamic constructors, destructors*/
		SpcDynamic();
		SpcDynamic(int id_in,int nodeid, int dof,IssmDouble value_in,int analysis_type);
		~SpcDynamic();

		/*Object virtual functions definitions*/
		Object *copy();
		void    DeepEcho();
		void    Echo();
		int     Id();
		void    Marshall(MarshallHandle* marshallhandle);
		int     ObjectEnum();

		/*Constraint virtual functions definitions*/
		void ActivatePenaltyMethod(void);
		void ConstrainNode(Nodes* nodes,Parameters* parameters);
		void PenaltyDofAndValue(int* dof,IssmDouble* value,Nodes* nodes,Parameters* parameters){_error_("not implemented yet");};
		void InputUpdateFromVectorDakota(IssmDouble* vector,Nodes* nodes,int name,int type){_error_("not implemented yet");};

		/*SpcDynamic management*/
		int        GetDof();
		int        GetNodeId();
		IssmDouble GetValue();
		void       SetDynamicConstraint(Nodes  *nodes,IssmDouble *yg_serial);

};

#endif  /* _SPCStatic_H_*/
