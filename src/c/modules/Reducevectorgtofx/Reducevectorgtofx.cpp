/*!\file Reducevectorgtofx
 * \brief reduce petsc vector from g set to s set (free dofs), using the nodeset partitioning
 * vectors.
 */

#include "./Reducevectorgtofx.h"

void Reducevectorgtofx(Vector<IssmDouble>** puf, Vector<IssmDouble>* ug, Nodes* nodes,Parameters* parameters){

	if(VerboseModule()) _printf0_("   Reduce vector from g to f set\n");

	/*first figure out fsize: */
	int fsize=nodes->NumberOfDofs(FsetEnum);
	int flocalsize = nodes->NumberOfDofsLocal(FsetEnum);

	/*If fsize is 0, return NULL vector*/
	if(fsize==0){
		Vector<IssmDouble>* uf=new Vector<IssmDouble>(0);
		*puf=uf;
		//*puf=NULL;
		return;
	}

	/*Get local vectors ug*/
	int        *indices_ug = NULL;
	IssmDouble *local_ug   = NULL;
	ug->GetLocalVector(&local_ug,&indices_ug);

	/*Allocate output*/
	Vector<IssmDouble>* uf=new Vector<IssmDouble>(flocalsize,fsize);

	/*Let nodes figure it out*/
	for(Object* & object : nodes->objects){
		Node* node=(Node*)object;
		node->VecReduce(uf,local_ug,indices_ug);
	}

	/*Assemble vector: */
	uf->Assemble();

	/*Cleanup and assing output pointer*/
	xDelete<int>(indices_ug);
	xDelete<IssmDouble>(local_ug);
	*puf=uf;
}
