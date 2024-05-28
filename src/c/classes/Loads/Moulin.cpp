/*!\file Moulin.c
 * \brief: implementation of the Moulin object
 */

/*Headers*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
#include "../../analyses/analyses.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES   1

/*Moulin constructors and destructor*/
Moulin::Moulin(){/*{{{*/
	this->parameters=NULL;
	this->hnode=NULL;
	this->hvertex=NULL;
	this->node=NULL;
	this->helement=NULL;
	this->element=NULL;
}
/*}}}*/
Moulin::Moulin(int id, int index, IoModel* iomodel){ /*{{{*/

	int pengrid_node_id;
	int pengrid_element_id;

	/*Some checks if debugging activated*/
	_assert_(iomodel->singlenodetoelementconnectivity);
	_assert_(index>=0 && index<iomodel->numberofvertices);
	_assert_(id);

	/*id: */
	this->id=id;

	/*hooks: */
	pengrid_node_id=index+1;
	pengrid_element_id=iomodel->singlenodetoelementconnectivity[index];
	_assert_(pengrid_element_id);

	this->hnode=new Hook(&pengrid_node_id,1);
	this->hvertex=new Hook(&pengrid_node_id,1);
	this->helement=new Hook(&pengrid_element_id,1);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->node=NULL;
	this->vertex=NULL;
	this->element=NULL;
}
/*}}}*/
Moulin::~Moulin(){/*{{{*/
	delete hnode;
	delete hvertex;
	delete helement;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Moulin::copy() {/*{{{*/

	Moulin* pengrid=NULL;

	pengrid=new Moulin();

	/*copy fields: */
	pengrid->id=this->id;

	/*point parameters: */
	pengrid->parameters=this->parameters;

	/*now deal with hooks and objects: */
	pengrid->hnode=(Hook*)this->hnode->copy();
	pengrid->hvertex=(Hook*)this->hvertex->copy();
	pengrid->helement=(Hook*)this->helement->copy();

	/*corresponding fields*/
	pengrid->node  =(Node*)pengrid->hnode->delivers();
	pengrid->vertex=(Vertex*)pengrid->hvertex->delivers();
	pengrid->element=(Element*)pengrid->helement->delivers();

	return pengrid;
}
/*}}}*/
void    Moulin::DeepEcho(void){/*{{{*/

	_printf_("Moulin:\n");
	_printf_("   id: " << id << "\n");
	hnode->DeepEcho();
	hvertex->DeepEcho();
	helement->DeepEcho();
	_printf_("   parameters\n");
	parameters->DeepEcho();
}
/*}}}*/
void    Moulin::Echo(void){/*{{{*/

	_printf_("Moulin:\n");
	_printf_("   id: " << id << "\n");
	hnode->Echo();
	hvertex->Echo();
	helement->Echo();
	_printf_("   parameters\n");
	parameters->Echo();
	//this->DeepEcho();
}
/*}}}*/
int     Moulin::Id(void){ return id; }/*{{{*/
/*}}}*/
void    Moulin::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	int object_enum = MoulinEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->id);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->hnode    = new Hook();
		this->hvertex  = new Hook();
		this->helement = new Hook();
	}

	this->hnode->Marshall(marshallhandle);
	this->hvertex->Marshall(marshallhandle);
	this->helement->Marshall(marshallhandle);

	/*corresponding fields*/
	node   =(Node*)this->hnode->delivers();
	vertex =(Vertex*)this->hvertex->delivers();
	element=(Element*)this->helement->delivers();
}
/*}}}*/
int     Moulin::ObjectEnum(void){/*{{{*/

	return MoulinEnum;
}
/*}}}*/

/*Load virtual functions definitions:*/
void  Moulin::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this load, ie links the objects in the hooks to their respective
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnode->configure(nodesin);
	hvertex->configure(verticesin);
	helement->configure(elementsin);

	/*Get corresponding fields*/
	node=(Node*)hnode->delivers();
	vertex=(Vertex*)hvertex->delivers();
	element=(Element*)helement->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
}
/*}}}*/
void  Moulin::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*recover some parameters*/
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case HydrologyGlaDSAnalysisEnum:
			Ke = this->CreateKMatrixHydrologyGlaDS();
			break;
		case HydrologyShaktiAnalysisEnum:
			/*do nothing: */
			return;
		case HydrologyDCInefficientAnalysisEnum:
			/*do nothing: */
			return;
		case HydrologyDCEfficientAnalysisEnum:
			/*do nothing: */
			return;
		default:
			_error_("Don't know why we should be here");

	}
	/*Add to global matrix*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}

}
/*}}}*/
void  Moulin::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case HydrologyGlaDSAnalysisEnum:
			pe = this->CreatePVectorHydrologyGlaDS();
			break;
		case HydrologyShaktiAnalysisEnum:
			pe = this->CreatePVectorHydrologyShakti();
			break;
		case HydrologyDCInefficientAnalysisEnum:
			pe = this->CreatePVectorHydrologyDCInefficient();
			break;
		case HydrologyDCEfficientAnalysisEnum:
			pe = this->CreatePVectorHydrologyDCEfficient();
			break;
		default:
			_error_("Don't know why we should be here");
			/*No loads applied, do nothing: */
			return;
	}
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}

}
/*}}}*/
void  Moulin::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(node);

	lidlist[0]=node->Lid();
}
/*}}}*/
void  Moulin::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(node);

	sidlist[0]=node->Sid();
}
/*}}}*/
int   Moulin::GetNumberOfNodes(void){/*{{{*/

	return NUMVERTICES;
}
/*}}}*/
bool  Moulin::IsPenalty(void){/*{{{*/
	return false;
}
/*}}}*/
void  Moulin::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*Don't do anything for now*/

}
/*}}}*/
void  Moulin::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*Don't do anything for now*/
}
/*}}}*/
void  Moulin::ResetHooks(){/*{{{*/

	this->node=NULL;
	this->element=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnode->reset();
	this->hvertex->reset();
	this->helement->reset();

}
/*}}}*/
void  Moulin::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Moulin::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	if(!flags[this->node->Lid()]){

		/*flag current node so that no other element processes it*/
		flags[this->node->Lid()]=true;

		flagsindices[flagsindices_counter[0]]=this->node->Lid();
		flagsindices_counter[0]++;

		/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
		switch(set2_enum){
			case FsetEnum:
				if(node->FSize()){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			case GsetEnum:
				if(node->gsize){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			case SsetEnum:
				if(node->SSize()){
					if(this->node->IsClone())
					 o_nz += 1;
					else
					 d_nz += 1;
				}
				break;
			default: _error_("not supported");
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/

/*Update virtual functions definitions:*/
void  Moulin::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/
	/*Nothing*/
}
/*}}}*/
void  Moulin::InputUpdateFromConstant(int constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromConstant(bool constant, int name){/*{{{*/

	/*Don't do anything for now*/
}
/*}}}*/
void  Moulin::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Moulin::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/

ElementMatrix* Moulin::CreateKMatrixHydrologyGlaDS(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;

	/*Initialize Element matrix*/
	ElementMatrix* Ke=new ElementMatrix(&node,1,this->parameters);

	/*Get all inputs and parameters*/
	IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	IssmDouble Am        = 0.; //For now...

	/*Load vector*/
	if(dt>0){
		Ke->values[0] = +Am/(rho_water*g)/dt;
	}

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyGlaDS(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;

	/*Initialize Element vector*/
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	/*Get all inputs and parameters*/
	IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	IssmDouble Am        = 0.; //For now...

	/*Get hydraulic potential*/
	IssmDouble phi_old,moulin_load;
	element->GetInputValue(&phi_old,node,HydraulicPotentialOldEnum);
	element->GetInputValue(&moulin_load,node,HydrologyMoulinInputEnum);

	pe->values[0] = moulin_load;
	if(dt>0.){
		pe->values[0] += Am/(rho_water*g) * phi_old/dt;
	}

	/*Clean up and return*/
	return pe;
}
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyShakti(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;

	IssmDouble moulin_load;

	/*Initialize Element matrix*/
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	this->element->GetInputValue(&moulin_load,node,HydrologyMoulinInputEnum);
	pe->values[0]=moulin_load;

	/*Clean up and return*/
	return pe;
}
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyDCInefficient(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/
	if(node->IsClone()) return NULL;
	bool isefficientlayer, active_element;
	IssmDouble moulin_load,dt;
	IssmDouble epl_active;
	/*Initialize Element matrix*/
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	this->element->GetInputValue(&moulin_load,node,HydrologydcBasalMoulinInputEnum);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	parameters->FindParam(&isefficientlayer,HydrologydcIsefficientlayerEnum);

	//Test version input in EPL when active
	if(isefficientlayer){
		this->element->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);
		if(!active_element){
			/* this->element->GetInputValue(&epl_active,node,HydrologydcMaskEplactiveNodeEnum); */
			/* if(reCast<bool>(epl_active))pe->values[0]=0.0; */
			/* else { */
			pe->values[0]=moulin_load*dt;
			/* 	if (moulin_load>0)_printf_("MoulinInput in Sed is "<<pe->values[0]<<"\n"); */
			/* 	if (moulin_load>0)pe->Echo(); */
			/* } */
			//if (node->Sid()==4)_printf_("MoulinInput in Sed is "<<moulin_load*dt<<"\n");
		}
		else pe->values[0]=0.0;
	}
	else pe->values[0]=moulin_load*dt;

	//Test only input in sed
	/* pe->values[0]=moulin_load*dt; */

	/*Clean up and return*/
	return pe;
}
/*}}}*/
ElementVector* Moulin::CreatePVectorHydrologyDCEfficient(void){/*{{{*/

	/*If this node is not the master node (belongs to another partition of the
	 * mesh), don't add the moulin input a second time*/

	if(node->IsClone()) return NULL;
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	//Test Input in epl if active
	/* IssmDouble epl_active; */
	/* this->element->GetInputValue(&epl_active,node,HydrologydcMaskEplactiveNodeEnum); */
	/* //if(node->Sid()==4)_printf_("Activity is "<<epl_active<<" \n"); */
	/* if(reCast<bool>(epl_active)){ */
	/* 	IssmDouble moulin_load,dt; */
	/* 	this->element->GetInputValue(&moulin_load,node,HydrologydcBasalMoulinInputEnum); */
	/* 	parameters->FindParam(&dt,TimesteppingTimeStepEnum); */
	/* 	pe->values[0]=moulin_load*dt; */
	/* 	if (moulin_load>0)_printf_("MoulinInput in Epl is "<<pe->values[1]<<"\n"); */

	/* } */
	/* 	else{ */
	/* 		pe->values[0]=0.0; */
	/* } */
	// Test element only test
	bool active_element;
	this->element->GetInputValue(&active_element,HydrologydcMaskEplactiveEltEnum);
	if(active_element){
		IssmDouble moulin_load,dt;
		this->element->GetInputValue(&moulin_load,node,HydrologydcBasalMoulinInputEnum);
		parameters->FindParam(&dt,TimesteppingTimeStepEnum);
		pe->values[0]=moulin_load*dt;
	}
	else pe->values[0]=0.0;


	//Test only input is sed
	/* pe->values[0]=0.0; */

	//Clean up and return
	return pe;
}
/*}}}*/
