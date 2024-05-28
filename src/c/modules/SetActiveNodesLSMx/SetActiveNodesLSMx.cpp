/*!\file GetMaskOfIceVerticesLSMx 
 * \brief: Return a mask for all the vertices determining whether the node should be active or not. 
 */

#include "./SetActiveNodesLSMx.h"

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../modules.h"

void SetActiveNodesLSMx(FemModel* femmodel,bool ishydrology,bool isdebris){/*{{{*/
	/* activate/deactivate nodes for levelset method according to IceMaskNodeActivation */

	/*Determine which node activation mask to pull from*/
	int nodeactivationmask = IceMaskNodeActivationEnum;
	if(ishydrology) nodeactivationmask = HydrologyMaskNodeActivationEnum;
	if(isdebris) nodeactivationmask = DebrisMaskNodeActivationEnum;


	for(Object* & object : femmodel->elements->objects){
		Element    *element  = xDynamicCast<Element*>(object);
		int         numnodes = element->GetNumberOfNodes();
		IssmDouble *mask     = xNew<IssmDouble>(numnodes);

		/*include switch for elements with multiple different sets of nodes*/
		switch(element->GetElementType()){
			case P1P1GLSEnum: case P1P1Enum:/* added to allow P1-P1 GLS */
			case MINIEnum:case MINIcondensedEnum:
			case TaylorHoodEnum:case XTaylorHoodEnum:case LATaylorHoodEnum:
			case CrouzeixRaviartEnum:case LACrouzeixRaviartEnum:case OneLayerP4zEnum:{
				Input* input=element->GetInput(nodeactivationmask);
				if(!input) _error_("Input " << EnumToStringx(nodeactivationmask) << " not found in element");

				/* Start looping on the number of vertices: */
				Gauss* gauss=element->NewGauss();
				for(int iv=0;iv<element->NumberofNodesVelocity();iv++){
					gauss->GaussNode(element->VelocityInterpolation(),iv);
					input->GetInputValue(&mask[iv],gauss);
				}
				for(int iv=0;iv<element->NumberofNodesPressure();iv++){
					gauss->GaussNode(element->PressureInterpolation(),iv);
					input->GetInputValue(&mask[element->NumberofNodesVelocity()+iv],gauss);
				}
				delete gauss;
				break;
			}
			default:
				element->GetInputListOnNodes(&mask[0],nodeactivationmask);
				break;
		}

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(mask[in]==1.) node->Activate();
			else             node->Deactivate();
		}
		xDelete<IssmDouble>(mask);
	}
}/*}}}*/

void GetMaskOfIceVerticesLSMx0(FemModel* femmodel,bool ishydrology,bool isdebris){/*{{{*/

	/*Determine which node activation to construct*/
	int nodeactivationmask = IceMaskNodeActivationEnum;
	if(ishydrology) nodeactivationmask = HydrologyMaskNodeActivationEnum;
	if(isdebris) nodeactivationmask = DebrisMaskNodeActivationEnum;

	/*Initialize vector with number of vertices*/
	int numvertices=femmodel->vertices->NumberOfVertices();
	if(numvertices==0)  return;

	int numvert_local = femmodel->vertices->NumberOfVerticesLocal();
	Vector<IssmDouble>* vec_mask_ice=new Vector<IssmDouble>(numvert_local,numvertices);

	/*Fill vector with values: */
	if(ishydrology){
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			if(element->IsIceInElement() && element->IsGrounded()){
				int nbv = element->GetNumberOfVertices();
				for(int iv=0;iv<nbv;iv++){
					vec_mask_ice->SetValue(element->vertices[iv]->Pid(),1.,INS_VAL);
				}
			}
		}
	}else if(isdebris){
		for(Object* & object : femmodel->elements->objects){
                        Element* element=xDynamicCast<Element*>(object);
                        if(element->IsIceInElement() && !element->IsAllMinThicknessInElement()){
                                int nbv = element->GetNumberOfVertices();
                                for(int iv=0;iv<nbv;iv++){
                                        vec_mask_ice->SetValue(element->vertices[iv]->Pid(),1.,INS_VAL);
                                }
                        }
                }
	}else{
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			if(element->IsIceInElement()){
				int nbv = element->GetNumberOfVertices();
				for(int iv=0;iv<nbv;iv++){
					vec_mask_ice->SetValue(element->vertices[iv]->Pid(),1.,INS_VAL);
				}
			}
		}
	}

	/*Assemble vector and serialize */
	vec_mask_ice->Assemble();
	InputUpdateFromVectorx(femmodel,vec_mask_ice,nodeactivationmask,VertexPIdEnum);
	delete vec_mask_ice;
}/*}}}*/
void GetMaskOfIceVerticesLSMx(FemModel* femmodel,bool ishydrology,bool isdebris){/*{{{*/

	/*Set configuration to levelset*/
	if(ishydrology){
		/*We may not be running with ismovingfront so we can't assume LevelsetAnalysis is active*/
		int hydrology_model;
		femmodel->parameters->FindParam(&hydrology_model,HydrologyModelEnum);
		if(hydrology_model==HydrologyshaktiEnum){
			femmodel->SetCurrentConfiguration(HydrologyShaktiAnalysisEnum);
		}
		else if(hydrology_model==HydrologyGlaDSEnum){
			femmodel->SetCurrentConfiguration(HydrologyGlaDSAnalysisEnum);
		}
		else{
			_error_("hydrology model not supported yet");
		}
	}else if(isdebris){
		femmodel->SetCurrentConfiguration(DebrisAnalysisEnum);
	}else{
		femmodel->SetCurrentConfiguration(LevelsetAnalysisEnum);
	}

	/*Determine which node activation to construct*/
	int nodeactivationmask = IceMaskNodeActivationEnum;
	if(ishydrology) nodeactivationmask = HydrologyMaskNodeActivationEnum;
	if(isdebris) nodeactivationmask = DebrisMaskNodeActivationEnum;

	/*Create vector on gset*/
	int gsize              = femmodel->nodes->NumberOfDofs(GsetEnum);
	int glocalsize_masters = femmodel->nodes->NumberOfDofsLocal(GsetEnum);
	if(gsize==0)  return;
	Vector<IssmDouble>* vec_mask_ice=new Vector<IssmDouble>(glocalsize_masters,gsize);

	/*Fill vector with values: */
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		if(ishydrology){
			if(element->IsIceInElement() && element->IsGrounded()){
				int numnodes = element->GetNumberOfNodes();
				int  gsize_local=GetNumberOfDofs(element->nodes,numnodes,GsetEnum,NoneEnum);
				int* glist_local=GetGlobalDofList(element->nodes,numnodes,GsetEnum,NoneEnum);
				IssmDouble* ones = xNew<IssmDouble>(gsize_local);
				for(int n=0;n<gsize_local;n++) ones[n] = 1.;
				vec_mask_ice->SetValues(gsize_local,glist_local,ones,INS_VAL);
				xDelete<IssmDouble>(ones);
				xDelete<int>(glist_local);
			}
		}else if(isdebris){
			if(element->IsIceInElement() && !element->IsAllMinThicknessInElement()){
                                int numnodes = element->GetNumberOfNodes();
                                int  gsize_local=GetNumberOfDofs(element->nodes,numnodes,GsetEnum,NoneEnum);
                                int* glist_local=GetGlobalDofList(element->nodes,numnodes,GsetEnum,NoneEnum);
                                IssmDouble* ones = xNew<IssmDouble>(gsize_local);
                                for(int n=0;n<gsize_local;n++) ones[n] = 1.;
                                vec_mask_ice->SetValues(gsize_local,glist_local,ones,INS_VAL);
                                xDelete<IssmDouble>(ones);
                                xDelete<int>(glist_local);
			}
		}else{
			if(element->IsIceInElement()){
				int numnodes = element->GetNumberOfNodes();
				int  gsize_local=GetNumberOfDofs(element->nodes,numnodes,GsetEnum,NoneEnum);
				int* glist_local=GetGlobalDofList(element->nodes,numnodes,GsetEnum,NoneEnum);
				IssmDouble* ones = xNew<IssmDouble>(gsize_local);
				for(int n=0;n<gsize_local;n++) ones[n] = 1.;
				vec_mask_ice->SetValues(gsize_local,glist_local,ones,INS_VAL);
				xDelete<IssmDouble>(ones);
				xDelete<int>(glist_local);
			}
		}
	}

	/*Assemble vector and serialize */
	vec_mask_ice->Assemble();

	/*Get local vector with masters and slaves*/
	IssmDouble *local_ug = NULL;
	femmodel->GetLocalVectorWithClonesGset(&local_ug,vec_mask_ice);
	delete vec_mask_ice;

	/*Now update inputs (analysis specific)*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputUpdateFromSolutionOneDof(local_ug,nodeactivationmask);
	}

	/*cleanup and return*/
	xDelete<IssmDouble>(local_ug);
}/*}}}*/
