/*!\file Pengrid.c
 * \brief: implementation of the Pengrid object
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

/*Pengrid constructors and destructor*/
Pengrid::Pengrid(){/*{{{*/
	this->parameters=NULL;
	this->hnode=NULL;
	this->node=NULL;
	this->helement=NULL;
	this->element=NULL;

	/*not active, not zigzagging: */
	active=0;
	zigzag_counter=0;

}
/*}}}*/
Pengrid::Pengrid(int id, int index, IoModel* iomodel){/*{{{*/

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
	this->helement=new Hook(&pengrid_element_id,1);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;
	this->node=NULL;
	this->element=NULL;

	//let's not forget internals
	this->active=0;
	this->zigzag_counter=0;

}
/*}}}*/
Pengrid::~Pengrid(){/*{{{*/
	delete hnode;
	delete helement;
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Pengrid::copy() {/*{{{*/

	Pengrid* pengrid=NULL;

	pengrid=new Pengrid();

	/*copy fields: */
	pengrid->id=this->id;

	/*point parameters: */
	pengrid->parameters=this->parameters;

	/*now deal with hooks and objects: */
	pengrid->hnode=(Hook*)this->hnode->copy();
	pengrid->helement=(Hook*)this->helement->copy();

	/*corresponding fields*/
	pengrid->node  =(Node*)pengrid->hnode->delivers();
	pengrid->element=(Element*)pengrid->helement->delivers();

	//let's not forget internals
	pengrid->active=this->active=0;
	pengrid->zigzag_counter=this->zigzag_counter=0;

	return pengrid;

}
/*}}}*/
void    Pengrid::DeepEcho(void){/*{{{*/

	_printf_("Pengrid:\n");
	_printf_("   id: " << id << "\n");
	hnode->DeepEcho();
	helement->DeepEcho();
	_printf_("   active " << this->active << "\n");
	_printf_("   zigzag_counter " << this->zigzag_counter << "\n");
	_printf_("   parameters\n");
	parameters->DeepEcho();
}
/*}}}*/
void    Pengrid::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int     Pengrid::Id(void){ return id; }/*{{{*/
/*}}}*/
void    Pengrid::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	int object_enum = PengridEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->id);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->hnode    = new Hook();
		this->helement = new Hook();
	}

	this->hnode->Marshall(marshallhandle);
	this->helement->Marshall(marshallhandle);

	/*corresponding fields*/
	node   =(Node*)this->hnode->delivers();
	element=(Element*)this->helement->delivers();

	marshallhandle->call(this->active);
	marshallhandle->call(this->zigzag_counter);
}/*}}}*/
int     Pengrid::ObjectEnum(void){/*{{{*/

	return PengridEnum;
}
/*}}}*/

/*Load virtual functions definitions:*/
void  Pengrid::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this load, ie links the objects in the hooks to their respective
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnode->configure(nodesin);
	helement->configure(elementsin);

	/*Get corresponding fields*/
	node=(Node*)hnode->delivers();
	element=(Element*)helement->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
}
/*}}}*/
void  Pengrid::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*No loads applied, do nothing: */
	return;

}
/*}}}*/
void  Pengrid::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/
	/*No loads applied, do nothing, originaly used for moulin input: */
	return;

}
/*}}}*/
void  Pengrid::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(node);

	lidlist[0]=node->Lid();
}
/*}}}*/
void  Pengrid::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(node);

	sidlist[0]=node->Sid();
}
/*}}}*/
int   Pengrid::GetNumberOfNodes(void){/*{{{*/

	return NUMVERTICES;
}
/*}}}*/
bool  Pengrid::IsPenalty(void){/*{{{*/
	return true;
}
/*}}}*/
void  Pengrid::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*Retrieve parameters: */
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case ThermalAnalysisEnum:
			Ke=PenaltyCreateKMatrixThermal(kmax);
			break;
		case MeltingAnalysisEnum:
			Ke=PenaltyCreateKMatrixMelting(kmax);
			break;
		case HydrologyDCInefficientAnalysisEnum:
			Ke=PenaltyCreateKMatrixHydrologyDCInefficient(kmax);
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global matrix*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}
}
/*}}}*/
void  Pengrid::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*Retrieve parameters: */
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case ThermalAnalysisEnum:
			pe=PenaltyCreatePVectorThermal(kmax);
			break;
		case MeltingAnalysisEnum:
			pe=PenaltyCreatePVectorMelting(kmax);
			break;
		case StressbalanceAnalysisEnum: case AdjointHorizAnalysisEnum:
			break;
		case HydrologyDCInefficientAnalysisEnum:
			pe=PenaltyCreatePVectorHydrologyDCInefficient(kmax);
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global Vector*/
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}
}
/*}}}*/
void  Pengrid::ResetHooks(){/*{{{*/

	this->node=NULL;
	this->element=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnode->reset();
	this->helement->reset();

}
/*}}}*/
void  Pengrid::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Pengrid::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum){/*{{{*/

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

/*Pengrid management:*/
void           Pengrid::ConstraintActivate(int* punstable){/*{{{*/

	int analysis_type;

	/*Retrieve parameters: */
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
	case StressbalanceAnalysisEnum:
		/*No penalty to check*/
		return;
	case ThermalAnalysisEnum:
		ConstraintActivateThermal(punstable);
		break;
	case MeltingAnalysisEnum:
		/*No penalty to check*/
		return;
	case HydrologyDCInefficientAnalysisEnum:
		ConstraintActivateHydrologyDCInefficient(punstable);
		break;
	default:
		_error_("analysis: " << EnumToStringx(analysis_type) << " not supported yet");
	}
}
/*}}}*/
void           Pengrid::ConstraintActivateHydrologyDCInefficient(int* punstable){/*{{{*/

	//   The penalty is stable if it doesn't change during two consecutive iterations.
	int        unstable=0;
	int        new_active;
	int        penalty_lock;
	IssmDouble pressure;
	IssmDouble h;
	IssmDouble h_max;
	HydrologyDCInefficientAnalysis* inefanalysis = NULL;

	/*check that pengrid is not a clone (penalty to be added only once)*/
	if(node->IsClone()){
		unstable=0;
		*punstable=unstable;
		return;
	}
	if(!element->IsOnBase()){
		unstable=0;
		active=0;
		*punstable=unstable;
		return;
	}

	/*Get sediment water head h*/
	inefanalysis = new HydrologyDCInefficientAnalysis();
	element->GetInputValue(&h,node,SedimentHeadSubstepEnum);
	inefanalysis->GetHydrologyDCInefficientHmax(&h_max,element,node);
	parameters->FindParam(&penalty_lock,HydrologydcPenaltyLockEnum);

	if (h>h_max){
	 new_active=1;
	}
	else{
	 new_active=0;
	}

	if(this->active==new_active){
		unstable=0;
	}
	else{
		unstable=1;
		if(penalty_lock)zigzag_counter++;
	}

	/*If penalty keeps zigzagging more than penalty_lock times: */
	if(penalty_lock){
		if(zigzag_counter>penalty_lock){
			unstable=0;
			active=1;
		}
	}
	/*Set penalty flag*/
	this->active=new_active;

	/*Assign output pointers:*/
	delete inefanalysis;
	*punstable=unstable;
}
/*}}}*/
void           Pengrid::ConstraintActivateThermal(int* punstable){/*{{{*/

	//   The penalty is stable if it doesn't change during to successive iterations.
	IssmDouble pressure;
	IssmDouble temperature;
	IssmDouble t_pmp;
	int        new_active;
	int        unstable=0;
	int        penalty_lock;

	/*recover pointers: */
	Penta* penta=(Penta*)element;

	/*check that pengrid is not a clone (penalty to be added only once)*/
	if (node->IsClone()){
		unstable=0;
		*punstable=unstable;
		return;
	}

	//First recover pressure and temperature values, using the element: */
	penta->GetInputValue(&pressure,node,PressureEnum);
	penta->GetInputValue(&temperature,node,TemperaturePicardEnum);

	//Recover our data:
	parameters->FindParam(&penalty_lock,ThermalPenaltyLockEnum);

	//Compute pressure melting point
	t_pmp=element->TMeltingPoint(pressure);

	//Figure out if temperature is over melting_point, in which case, this penalty needs to be activated.

	if (temperature>t_pmp){
		new_active=1;
	}
	else{
		new_active=0;
	}

	//Figure out stability of this penalty
	if (active==new_active){
		unstable=0;
	}
	else{
		unstable=1;
		if(penalty_lock)zigzag_counter++;
	}

	/*If penalty keeps zigzagging more than 5 times: */
	if(penalty_lock){
		if(zigzag_counter>penalty_lock){
			unstable=0;
			active=1;
		}
	}

	//Set penalty flag
	active=new_active;

	//*Assign output pointers:*/
	*punstable=unstable;
}
/*}}}*/
ElementMatrix* Pengrid::PenaltyCreateKMatrixHydrologyDCInefficient(IssmDouble kmax){/*{{{*/
	IssmDouble    penalty_factor;

	/*Retrieve parameters*/
	parameters->FindParam(&penalty_factor,HydrologydcPenaltyFactorEnum);

	/*Initialize Element matrix and return if necessary*/
	if(!this->active) return NULL;
	ElementMatrix* Ke=new ElementMatrix(&node,NUMVERTICES,this->parameters);

	Ke->values[0]=kmax*pow(10.,penalty_factor);

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementMatrix* Pengrid::PenaltyCreateKMatrixMelting(IssmDouble kmax){/*{{{*/

	IssmDouble pressure,temperature,t_pmp;
	IssmDouble penalty_factor;

	Penta* penta=(Penta*)element;

	/*check that pengrid is not a clone (penalty to be added only once)*/
	if (node->IsClone()) return NULL;
	ElementMatrix* Ke=new ElementMatrix(&node,1,this->parameters);

	/*Retrieve all parameters*/
	penta->GetInputValue(&pressure,node,PressureEnum);
	penta->GetInputValue(&temperature,node,TemperatureEnum);
	parameters->FindParam(&penalty_factor,ThermalPenaltyFactorEnum);

	/*Compute pressure melting point*/
	t_pmp=parameters->FindParam(MaterialsMeltingpointEnum)-parameters->FindParam(MaterialsBetaEnum)*pressure;

	/*Add penalty load*/
	if (temperature<t_pmp){ //If T<Tpmp, there must be no melting. Therefore, melting should be  constrained to 0 when T<Tpmp, instead of using spcs, use penalties
		Ke->values[0]=kmax*pow(10.,penalty_factor);
	}

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementMatrix* Pengrid::PenaltyCreateKMatrixThermal(IssmDouble kmax){/*{{{*/

	IssmDouble    penalty_factor;

	/*Initialize Element matrix and return if necessary*/
	if(!this->active) return NULL;
	ElementMatrix* Ke=new ElementMatrix(&node,NUMVERTICES,this->parameters);

	/*recover parameters: */
	parameters->FindParam(&penalty_factor,ThermalPenaltyFactorEnum);

	Ke->values[0]=kmax*pow(10.,penalty_factor);

	/*Clean up and return*/
	return Ke;
}
/*}}}*/
ElementVector* Pengrid::PenaltyCreatePVectorHydrologyDCInefficient(IssmDouble kmax){/*{{{*/

	IssmDouble h_max;
	IssmDouble penalty_factor;
	HydrologyDCInefficientAnalysis* inefanalysis = NULL;

	/*Initialize Element matrix and return if necessary*/
	if(!this->active) return NULL;
	ElementVector* pe=new ElementVector(&node,1,this->parameters);
	inefanalysis = new HydrologyDCInefficientAnalysis();

	/*Retrieve parameters*/
	parameters->FindParam(&penalty_factor,HydrologydcPenaltyFactorEnum);

	/*Get h_max and compute penalty*/
	inefanalysis->GetHydrologyDCInefficientHmax(&h_max,element,node);

	pe->values[0]=kmax*pow(10.,penalty_factor)*h_max;

	/*Clean up and return*/
	delete inefanalysis;
	return pe;
}
/*}}}*/
ElementVector* Pengrid::PenaltyCreatePVectorMelting(IssmDouble kmax){/*{{{*/

	IssmDouble pressure;
	IssmDouble temperature;
	IssmDouble melting_offset;
	IssmDouble t_pmp;
	IssmDouble dt,penalty_factor;

	/*recover pointers: */
	Penta* penta=(Penta*)element;

	/*check that pengrid is not a clone (penalty to be added only once)*/
	if (node->IsClone()) return NULL;
	ElementVector* pe=new ElementVector(&node,NUMVERTICES,this->parameters);

	/*Retrieve all parameters*/
	penta->GetInputValue(&pressure,node,PressureEnum);
	penta->GetInputValue(&temperature,node,TemperatureEnum);
	parameters->FindParam(&melting_offset,MeltingOffsetEnum);
	parameters->FindParam(&dt,TimesteppingTimeStepEnum);
	parameters->FindParam(&penalty_factor,ThermalPenaltyFactorEnum);

	/*Compute pressure melting point*/
	t_pmp=parameters->FindParam(MaterialsMeltingpointEnum)-parameters->FindParam(MaterialsBetaEnum)*pressure;

	/*Add penalty load
	  This time, the penalty must have the same value as the one used for the thermal computation
	  so that the corresponding melting can be computed correctly
	  In the thermal computation, we used kmax=melting_offset, and the same penalty_factor*/
	if (temperature<t_pmp){ //%no melting
		pe->values[0]=0;
	}
	else{
		if (reCast<bool>(dt)) pe->values[0]=melting_offset*pow(10.,penalty_factor)*(temperature-t_pmp)/dt;
		else    pe->values[0]=melting_offset*pow(10.,penalty_factor)*(temperature-t_pmp);
	}

	/*Clean up and return*/
	return pe;
}
/*}}}*/
ElementVector* Pengrid::PenaltyCreatePVectorThermal(IssmDouble kmax){/*{{{*/

	IssmDouble pressure;
	IssmDouble t_pmp;
	IssmDouble penalty_factor;

	Penta* penta=(Penta*)element;

	/*Initialize Element matrix and return if necessary*/
	if(!this->active) return NULL;
	ElementVector* pe=new ElementVector(&node,1,this->parameters);

	/*Retrieve all parameters*/
	penta->GetInputValue(&pressure,node,PressureEnum);
	parameters->FindParam(&penalty_factor,ThermalPenaltyFactorEnum);

	/*Compute pressure melting point*/
	t_pmp=parameters->FindParam(MaterialsMeltingpointEnum)-parameters->FindParam(MaterialsBetaEnum)*pressure;

	pe->values[0]=kmax*pow(10.,penalty_factor)*t_pmp;

	/*Clean up and return*/
	return pe;
}
/*}}}*/
void           Pengrid::ResetConstraint(void){/*{{{*/
	active         = 0;
	zigzag_counter = 0;
}
/*}}}*/
void           Pengrid::ResetZigzagCounter(){/*{{{*/

	zigzag_counter=0;
}
/*}}}*/
