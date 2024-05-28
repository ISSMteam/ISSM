/*!\file:  Constraint.h
 * \brief abstract class for Constraint object
 * This class is a place holder for constraints
 * It is derived from Object, so DataSets can contain them.
 */ 

#ifndef _CONSTRAINT_H_
#define _CONSTRAINT_H_

/*Headers:*/
/*{{{*/
class Nodes;
#include "../../datastructures/datastructures.h"
#include "../../toolkits/toolkits.h"
/*}}}*/

class Constraint: public Object{

	public: 

		virtual      ~Constraint(){};
		virtual void ActivatePenaltyMethod(void)=0;
		virtual void ConstrainNode(Nodes* nodes,Parameters* parameters)=0;
		virtual void PenaltyDofAndValue(int* dof,IssmDouble* value,Nodes* nodes,Parameters* parameters)=0;
		virtual void InputUpdateFromVectorDakota(IssmDouble* vector,Nodes* nodes,int name,int type) = 0;

};
#endif
