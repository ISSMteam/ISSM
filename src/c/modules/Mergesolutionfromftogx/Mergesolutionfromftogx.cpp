/*!\file Mergesolutionfromftogx
 * \brief merge solution back from f set into g set
 */

#include "../../shared/io/io.h"
#include "./Mergesolutionfromftogx.h"

void	Mergesolutionfromftogx( Vector<IssmDouble>** pug, Vector<IssmDouble>* uf, Vector<IssmDouble>* ys, Nodes* nodes, Parameters* parameters, bool flag_ys0){

	/*Display message*/
	if(VerboseModule()) _printf0_("   Merging solution vector from fset to gset\n");

	/*first, get gsize, fsize and ssize: */
	int gsize=nodes->NumberOfDofs(GsetEnum);
	int gsize_local=nodes->NumberOfDofsLocal(GsetEnum);
	int fsize=nodes->NumberOfDofs(FsetEnum);
	int ssize=nodes->NumberOfDofs(SsetEnum);

	/*serialize uf and ys: those two vectors will be indexed by the nodes, who are the only ones
	 *that know which values should be plugged into ug and where: */
	if(ssize) if(flag_ys0) ys->Set(0.0);

	/*Get local vectors ys and uf*/
	int        *indices_ys = NULL;
	IssmDouble *local_ys   = NULL;
	ys->GetLocalVector(&local_ys,&indices_ys);
	int        *indices_uf = NULL;
	IssmDouble *local_uf   = NULL;
	uf->GetLocalVector(&local_uf,&indices_uf);

	/*initialize ug: */
	Vector<IssmDouble>* ug=new Vector<IssmDouble>(gsize_local,gsize);

	/*Let nodes figure it out*/
	for(Object* & object: nodes->objects){
		Node* node=xDynamicCast<Node*>(object);
		node->VecMerge(ug,local_uf,indices_uf,local_ys,indices_ys);
	}

	/*Assemble vector: */
	ug->Assemble();

	/*Cleanup and assign output pointer*/
	xDelete<int>(indices_uf);
	xDelete<int>(indices_ys);
	xDelete<IssmDouble>(local_uf);
	xDelete<IssmDouble>(local_ys);
	*pug=ug;
}
