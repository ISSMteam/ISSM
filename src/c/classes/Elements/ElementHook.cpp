/*!\file ElementHook.c
 * \brief: implementation of the ElementHook object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
/*}}}*/

/*Object constructors and destructor*/
ElementHook::ElementHook(){/*{{{*/
	numanalyses=UNDEF;
	this->hnodes     = NULL;
	this->hvertices  = NULL;
	this->hmaterial  = NULL;
	this->hneighbors = NULL;
}
/*}}}*/
ElementHook::~ElementHook(){/*{{{*/

	if(this->hnodes){
		for(int i=0;i<this->numanalyses;i++){
			if(this->hnodes[i]) delete this->hnodes[i]; 
		}
		delete [] this->hnodes;
	}
	delete hvertices;
	delete hmaterial;
	delete hneighbors;
}
/*}}}*/
ElementHook::ElementHook(int in_numanalyses,int element_id,int numvertices,IoModel* iomodel){/*{{{*/

	/*retrieve material_id*/
	int material_id;
	material_id = element_id;

	/*retrieve vertices ids*/
	int* vertex_ids = xNew<int>(numvertices);
	for(int i=0;i<numvertices;i++){ 
		vertex_ids[i]=reCast<int>(iomodel->elements[(element_id-1)*numvertices+i]);
	}

	this->numanalyses = in_numanalyses;
	this->hnodes      = new Hook*[in_numanalyses];
	this->hvertices   = new Hook(&vertex_ids[0],numvertices);
	this->hmaterial   = new Hook(&material_id,1);
	this->hneighbors  = NULL;

	/*Initialize hnodes as NULL*/
	for(int i=0;i<this->numanalyses;i++){
		this->hnodes[i]=NULL;
	}

	/*Clean up*/
	xDelete<int>(vertex_ids);

}
/*}}}*/
void ElementHook::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int i;
	bool* hnodesi_null=NULL; /*intermediary needed*/
	bool  hnodes_null=true; /*this could be NULL on empty constructor*/
	bool  hneighbors_null=true; /*don't deal with hneighbors, unless explicitely asked to*/

	_assert_(this);

	/*preliminary, before marshall starts: */
	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){
		if(this->hneighbors)hneighbors_null=false;
		if(this->hnodes){
			hnodes_null=false;
			hnodesi_null=xNew<bool>(numanalyses);
			for(i=0;i<numanalyses;i++){
				if(this->hnodes[i])hnodesi_null[i]=false;
				else hnodesi_null[i]=true;
			}
		}
	}

	/*ok, marshall operations: */
	int object_enum = ElementHookEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->numanalyses);
	marshallhandle->call(hneighbors_null);
	marshallhandle->call(hnodes_null);
	marshallhandle->call(hnodesi_null,numanalyses);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){

		if (!hnodes_null)this->hnodes = new Hook*[numanalyses];
		else this->hnodes=NULL;
		this->hvertices   = new Hook();
		this->hmaterial   = new Hook();
		if(!hneighbors_null)this->hneighbors  = new Hook();
		else this->hneighbors=NULL;

		/*Initialize hnodes: */
		if (this->hnodes){
			for(int i=0;i<this->numanalyses;i++){
				if(!hnodesi_null[i])this->hnodes[i]=new Hook();
				else this->hnodes[i]=NULL;
			}
		}
	}

	if (this->hnodes){ 
		for (i=0;i<numanalyses;i++) if(this->hnodes[i])this->hnodes[i]->Marshall(marshallhandle);
	}
	this->hvertices->Marshall(marshallhandle);
	this->hmaterial->Marshall(marshallhandle);
	if(this->hneighbors)this->hneighbors->Marshall(marshallhandle);

	/*Free resources: */
	if(hnodesi_null) xDelete<bool>(hnodesi_null);

}
/*}}}*/

void ElementHook::DeepEcho(){/*{{{*/

	_printf_(" ElementHook DeepEcho:\n");
	_printf_("  numanalyses : "<< this->numanalyses <<"\n");

	_printf_("  hnodes:\n");
	if(hnodes){
		for(int i=0;i<this->numanalyses;i++) {
			if(hnodes[i]) hnodes[i]->DeepEcho();
			else _printf_("  hnodes["<< i << "] = NULL\n"); 
		}
	}
	else _printf_("  hnodes = NULL\n");

	_printf_("  hvertices:\n");
	if(hvertices) hvertices->DeepEcho();
	else _printf_("  hvertices = NULL\n");

	_printf_("  hmaterial:\n");
	if(hmaterial) hmaterial->DeepEcho();
   else _printf_("  hmaterial = NULL\n");

	_printf_("  hneighbors:\n");
	if(hneighbors) hneighbors->DeepEcho();
   else _printf_("  hneighbors = NULL\n");

	return;
}
/*}}}*/
void ElementHook::Echo(){/*{{{*/

	_printf_(" ElementHook Echo:\n");
	_printf_("  numanalyses : "<< this->numanalyses <<"\n");

	_printf_("  hnodes:\n");
	if(hnodes){
		for(int i=0;i<this->numanalyses;i++) {
			if(hnodes[i]) hnodes[i]->Echo();
		}
	}
	else _printf_("  hnodes = NULL\n");

	_printf_("  hvertices:\n");
	if(hvertices) hvertices->Echo();
	else _printf_("  hvertices = NULL\n");

	_printf_("  hmaterial:\n");
	if(hmaterial) hmaterial->Echo();
   else _printf_("  hmaterial = NULL\n");

	_printf_("  hneighbors:\n");
	if(hneighbors) hneighbors->Echo();
   else _printf_("  hneighbors = NULL\n");

	return;
}
/*}}}*/
void ElementHook::InitHookNeighbors(int* element_ids){/*{{{*/
	this->hneighbors=new Hook(element_ids,2);
}
/*}}}*/
void ElementHook::SetHookNodes(int* node_ids,int numnodes,int analysis_counter){/*{{{*/
	if(this->hnodes) this->hnodes[analysis_counter]= new Hook(node_ids,numnodes);
}
/*}}}*/
void ElementHook::SpawnSegHook(ElementHook* triahook,int index1,int index2){/*{{{*/

	triahook->numanalyses=this->numanalyses;

	int indices[2];
	indices[0]=index1;
	indices[1]=index2;

	/*Spawn nodes hook*/
	triahook->hnodes=new Hook*[this->numanalyses];
	for(int i=0;i<this->numanalyses;i++){
		/*Do not do anything if Hook is empty*/
		if (!this->hnodes[i] || this->hnodes[i]->GetNum()==0){
			triahook->hnodes[i]=NULL;
		}
		else{
			triahook->hnodes[i]=this->hnodes[i]->Spawn(indices,2);
		}
	}

	/*do not spawn hmaterial. material will be taken care of by Tria*/
	triahook->hmaterial=NULL;
	triahook->hvertices=(Hook*)this->hvertices->Spawn(indices,2);
}
/*}}}*/
void ElementHook::SpawnTriaHook(ElementHook* triahook,int index1,int index2,int index3){/*{{{*/

	/*Create arrow of indices depending on location (0=base 1=surface)*/
	int indices[3];
	indices[0] = index1;
	indices[1] = index2;
	indices[2] = index3;

	triahook->numanalyses=this->numanalyses;

	/*Spawn nodes hook*/
	triahook->hnodes=new Hook*[this->numanalyses];
	for(int i=0;i<this->numanalyses;i++){
		/*Do not do anything if Hook is empty*/
		if (!this->hnodes[i] || this->hnodes[i]->GetNum()==0){
			triahook->hnodes[i]=NULL;
		}
		else{
			triahook->hnodes[i]=this->hnodes[i]->Spawn(indices,3);
		}
	}

	/*do not spawn hmaterial. material will be taken care of by Penta*/
	triahook->hmaterial=NULL;
	triahook->hvertices=(Hook*)this->hvertices->Spawn(indices,3);
}
/*}}}*/
