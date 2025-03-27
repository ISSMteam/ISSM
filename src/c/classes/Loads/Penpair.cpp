/*!\file Penpair.c
 * \brief: implementation of the Penpair object
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
/*}}}*/

/*Element macros*/
#define NUMVERTICES 2

/*Penpair constructors and destructor*/
Penpair::Penpair(){/*{{{*/

	this->hnodes=NULL;
	this->nodes=NULL;
	this->parameters=NULL;
	return;
}
/*}}}*/
Penpair::Penpair(int penpair_id, int* penpair_node_ids){/*{{{*/

	this->id=penpair_id;
	this->hnodes=new Hook(penpair_node_ids,2);
	this->parameters=NULL;
	this->nodes=NULL;

	return;
}
/*}}}*/
Penpair::~Penpair(){/*{{{*/
	delete hnodes;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Penpair::copy() {/*{{{*/

	Penpair* penpair=NULL;

	penpair=new Penpair();

	/*copy fields: */
	penpair->id=this->id;

	/*now deal with hooks and objects: */
	penpair->hnodes=(Hook*)this->hnodes->copy();
	penpair->nodes =(Node**)penpair->hnodes->deliverp();

	/*point parameters: */
	penpair->parameters=this->parameters;

	return penpair;

}
/*}}}*/
void    Penpair::DeepEcho(void){/*{{{*/

	_printf_("Penpair:\n");
	_printf_("   id: " << id << "\n");
	hnodes->DeepEcho();

	return;
}		
/*}}}*/
void    Penpair::Echo(void){/*{{{*/

	_printf_("Penpair:\n");
	_printf_("   id: " << id << "\n");
	hnodes->Echo();

	return;
}
/*}}}*/
int     Penpair::Id(void){ return id; }/*{{{*/
/*}}}*/
void    Penpair::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	int object_enum = PenpairEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->id);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->hnodes = new Hook();
	}
	this->hnodes->Marshall(marshallhandle);

	/*corresponding fields*/
	nodes = (Node**)this->hnodes->deliverp();
}
/*}}}*/
int     Penpair::ObjectEnum(void){/*{{{*/

	return PenpairEnum;
}
/*}}}*/

/*Load virtual functions definitions:*/
void  Penpair::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnodes->configure((DataSet*)nodesin);

	/*Initialize hooked fields*/
	this->nodes  =(Node**)hnodes->deliverp();

	/*point parameters to real dataset: */
	this->parameters=parametersin;

}
/*}}}*/
void  Penpair::CreateJacobianMatrix(Matrix<IssmDouble>* Jff){/*{{{*/
	this->CreateKMatrix(Jff,NULL);
}
/*}}}*/
void  Penpair::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/
	/*If you code this piece, don't forget that a penalty will be inactive if it is dealing with clone nodes*/
	/*No loads applied, do nothing: */
	return;

}
/*}}}*/
void  Penpair::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	/*No loads applied, do nothing: */
	return;

}
/*}}}*/
void  Penpair::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	for(int i=0;i<NUMVERTICES;i++) lidlist[i]=nodes[i]->Lid();
}
/*}}}*/
void  Penpair::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	for(int i=0;i<NUMVERTICES;i++) sidlist[i]=nodes[i]->Sid();
}
/*}}}*/
int   Penpair::GetNumberOfNodes(void){/*{{{*/

	return NUMVERTICES;
}
/*}}}*/
bool  Penpair::IsPenalty(void){/*{{{*/
	return true;
}
/*}}}*/
void  Penpair::PenaltyCreateJacobianMatrix(Matrix<IssmDouble>* Jff,IssmDouble kmax){/*{{{*/
	this->PenaltyCreateKMatrix(Jff,NULL,kmax);
}
/*}}}*/
void  Penpair::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*Retrieve parameters: */
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case StressbalanceAnalysisEnum:
			Ke=PenaltyCreateKMatrixStressbalanceHoriz(kmax);
			break;
		case MasstransportAnalysisEnum:
		case FreeSurfaceBaseAnalysisEnum:
		case FreeSurfaceTopAnalysisEnum:
			Ke=PenaltyCreateKMatrixMasstransport(kmax);
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global Vector*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
}
/*}}}*/
void  Penpair::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/
	/*No loads applied, do nothing: */
	return;
}
/*}}}*/
void  Penpair::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnodes->reset();

}
/*}}}*/
void  Penpair::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Penpair::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<NUMVERTICES;i++){

		if(!flags[this->nodes[i]->Lid()]){

			/*flag current node so that no other element processes it*/
			flags[this->nodes[i]->Lid()]=true;

			flagsindices[flagsindices_counter[0]]=this->nodes[i]->Lid();
         flagsindices_counter[0]++;
			int counter=0;

			/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
			switch(set2_enum){
				case FsetEnum:
					if(nodes[i]->FSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case GsetEnum:
					if(nodes[i]->gsize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case SsetEnum:
					if(nodes[i]->SSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				default: _error_("not supported");
			}
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/

/*Update virtual functions definitions:*/
void  Penpair::InputUpdateFromConstant(IssmDouble constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Penpair::InputUpdateFromConstant(int constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Penpair::InputUpdateFromConstant(bool constant, int name){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/
void  Penpair::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/
	/*Nothing updated yet*/
}
/*}}}*/

/*Penpair management:*/
ElementMatrix* Penpair::PenaltyCreateKMatrixMasstransport(IssmDouble kmax){/*{{{*/

	/*Initialize Element vector and return if necessary*/
	const int      numdof = NUMVERTICES*1;
	ElementMatrix* Ke     = new ElementMatrix(nodes,NUMVERTICES,this->parameters);

	/*recover parameters: */
	IssmDouble penalty_factor = parameters->FindParam(MasstransportPenaltyFactorEnum);

	/*Create elementary matrix: add penalty to*/
	Ke->values[0*numdof+0]=+kmax*pow(10.,penalty_factor);
	Ke->values[0*numdof+1]=-kmax*pow(10.,penalty_factor);
	Ke->values[1*numdof+0]=-kmax*pow(10.,penalty_factor);
	Ke->values[1*numdof+1]=+kmax*pow(10.,penalty_factor);

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementMatrix* Penpair::PenaltyCreateKMatrixStressbalanceFS(IssmDouble kmax){/*{{{*/

	int        numdof,numdof2,N;
	IssmDouble penalty_offset;

	/*Initialize Element vector and return if necessary*/
	ElementMatrix* Ke=new ElementMatrix(nodes,NUMVERTICES,this->parameters,FSvelocityEnum);

	/*recover parameters: */
	parameters->FindParam(&penalty_offset,StressbalancePenaltyFactorEnum);

	/*Get number of dof for these two nodes*/
	numdof =this->nodes[0]->GetNumberOfDofs(FSApproximationEnum,GsetEnum);
	numdof2=this->nodes[1]->GetNumberOfDofs(FSApproximationEnum,GsetEnum);
	N=NUMVERTICES*numdof;

	/*Add penalty to Element matrix*/
	for(int i=0;i<numdof;i++){
		Ke->values[         i*N+i       ]=+kmax*pow(10.,penalty_offset);
		Ke->values[         i*N+numdof+i]=-kmax*pow(10.,penalty_offset);
		Ke->values[(numdof+i)*N+i       ]=-kmax*pow(10.,penalty_offset);
		Ke->values[(numdof+i)*N+numdof+i]=+kmax*pow(10.,penalty_offset);
	}

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementMatrix* Penpair::PenaltyCreateKMatrixStressbalanceHoriz(IssmDouble kmax){/*{{{*/

	int    approximation0=nodes[0]->GetApproximation();
	int    approximation1=nodes[1]->GetApproximation();

	switch(approximation0){
		case SSAApproximationEnum:
			switch(approximation1){
				case SSAApproximationEnum: return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				case HOApproximationEnum:   return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				default: _error_("Approximation "<<EnumToStringx(approximation1)<<" not supported yet");
			}
		case HOApproximationEnum:
			switch(approximation1){
				case SSAApproximationEnum: return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				case HOApproximationEnum:   return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				default: _error_("Approximation "<<EnumToStringx(approximation1)<<" not supported yet");
			}
		case MOLHOApproximationEnum:
			switch(approximation1){
				case MOLHOApproximationEnum:   return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				default: _error_("Approximation "<<EnumToStringx(approximation1)<<" not supported yet");
			}
		case L1L2ApproximationEnum:
			switch(approximation1){
				case L1L2ApproximationEnum:   return PenaltyCreateKMatrixStressbalanceSSAHO(kmax); 
				default: _error_("Approximation "<<EnumToStringx(approximation1)<<" not supported yet");
			}
		case FSvelocityEnum:
			switch(approximation1){
				case FSvelocityEnum: return PenaltyCreateKMatrixStressbalanceFS(kmax); 
				case NoneApproximationEnum: return   PenaltyCreateKMatrixStressbalanceFS(kmax); 
				default: _error_("Approximation "<<EnumToStringx(approximation1)<<" not supported yet");
			}
		case NoneApproximationEnum:
			switch(approximation1){
				case FSvelocityEnum: return PenaltyCreateKMatrixStressbalanceFS(kmax); 
				case NoneApproximationEnum: return   PenaltyCreateKMatrixStressbalanceFS(kmax); 
			}
		default: _error_("Approximation "<<EnumToStringx(approximation0)<<" not supported yet");
	}
}
/*}}}*/
ElementMatrix* Penpair::PenaltyCreateKMatrixStressbalanceSSAHO(IssmDouble kmax){/*{{{*/

	int        numdof,numdof2,N;
	IssmDouble penalty_offset;

	/*Initialize Element vector and return if necessary*/
	ElementMatrix* Ke=new ElementMatrix(nodes,NUMVERTICES,this->parameters);

	/*recover parameters: */
	parameters->FindParam(&penalty_offset,StressbalancePenaltyFactorEnum);

	/*Get number of dof for these two nodes*/
	numdof =this->nodes[0]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
	numdof2=this->nodes[1]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
	_assert_(numdof==numdof2);
	N=NUMVERTICES*numdof;

	/*Add penalty to Element matrix*/
	for(int i=0;i<numdof;i++){
		Ke->values[         i*N+i       ]=+kmax*pow(10.,penalty_offset);
		Ke->values[         i*N+numdof+i]=-kmax*pow(10.,penalty_offset);
		Ke->values[(numdof+i)*N+i       ]=-kmax*pow(10.,penalty_offset);
		Ke->values[(numdof+i)*N+numdof+i]=+kmax*pow(10.,penalty_offset);
	}

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
