/*!\file CreateNodalConstraintsx
 * \brief: establish degrees of freedom for all nodes, and return partitioning vector. Do only once.
 */

#include "./CreateNodalConstraintsx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void CreateNodalConstraintsx( Vector<IssmDouble>** pys, Nodes* nodes){

	bool  oldalloc  = false;

	/*output: */
	Vector<IssmDouble>* ys=NULL;
	int ssize;
	int slocalsize;

	if(VerboseModule()) _printf0_("   Create nodal constraints\n");

	/*figure out how many dofs we have: */
	ssize=nodes->NumberOfDofs(SsetEnum);
	slocalsize = nodes->NumberOfDofsLocal(SsetEnum);

	/*allocate:*/
	if(oldalloc)
	 ys=new Vector<IssmDouble>(ssize);
	else
	 ys=new Vector<IssmDouble>(slocalsize,ssize);

	/*go through all nodes, and for the ones corresponding to this configuration_type, fill the
	 * constraints vector with the constraint values: */
	for(Object* & object: nodes->objects){
		Node* node=xDynamicCast<Node*>(object);
		node->CreateNodalConstraints(ys);
	}

	/*Assemble: */
	ys->Assemble();

	/*Assign output pointers: */
	*pys=ys;
}
