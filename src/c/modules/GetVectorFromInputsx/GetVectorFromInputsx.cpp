/*!\file GetVectorFromInputsx
 * \brief retrieve vector from inputs in elements
 */

#include "./GetVectorFromInputsx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void GetVectorFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type){ /*{{{*/

	Vector<IssmDouble>* vector=NULL;

	switch(type){
		case ElementSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->elements->NumberOfElements());
			break;
		case VertexPIdEnum: case VertexSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->vertices->NumberOfVertices());
			break;
		case NodesEnum:case NodeSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->nodes->NumberOfNodes());
			break;
		default:
			_error_("vector type: " << EnumToStringx(type) << " not supported yet!");
	}
	/*Look up in elements*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->GetVectorFromInputs(vector,name,type);
	}

	vector->Assemble();

	/*Assign output pointers:*/
	*pvector=vector;
} /*}}}*/
void GetVectoronBaseFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type){ /*{{{*/

	int domaintype;
	Vector<IssmDouble>* vector=NULL;

	switch(type){
		case ElementSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->elements->NumberOfElements());
			break;
		case VertexPIdEnum: case VertexSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->vertices->NumberOfVertices());
			break;
		case NodesEnum:case NodeSIdEnum:
			vector=new Vector<IssmDouble>(femmodel->nodes->NumberOfNodes());
			break;
		default:
			_error_("vector type: " << EnumToStringx(type) << " not supported yet!");
	}

	/*Look up in elements*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->FindParam(&domaintype,DomainTypeEnum);
		switch(domaintype){
			case Domain2DhorizontalEnum:
				element->GetVectorFromInputs(vector,name,type);
				break;
			case Domain3DEnum:
				if(!element->IsOnBase()) continue;
				element->GetVectorFromInputs(vector,name,type);
				break;
			default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
		}
	}
	vector->Assemble();
	/*Assign output pointers:*/
	*pvector=vector;
} /*}}}*/
void GetVectorFromInputsx(Vector<IssmDouble>** pvector,FemModel* femmodel,int name,int type,IssmDouble time){/*{{{*/

	Vector<IssmDouble>* vector=NULL;

	switch(type){
	case VertexPIdEnum: case VertexSIdEnum:
		vector=new Vector<IssmDouble>(femmodel->vertices->NumberOfVertices());
		break;
	case NodesEnum:case NodeSIdEnum:
		vector=new Vector<IssmDouble>(femmodel->nodes->NumberOfNodes());
		break;
	default:
			_error_("vector type: " << EnumToStringx(type) << " not supported yet!");
	}
	/*Look up in elements*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->GetVectorFromInputs(vector,name,type,time);
	}

	vector->Assemble();

	/*Assign output pointers:*/
	*pvector=vector;
}/*}}}*/
void GetVectorFromInputsx(IssmDouble** pvector,FemModel* femmodel,int name, int type){/*{{{*/

	/*output: */
	IssmDouble* vector=NULL;

	/*intermediary: */
	Vector<IssmDouble>* vec_vector=NULL;

	GetVectorFromInputsx(&vec_vector,femmodel,name,type);
	vector=vec_vector->ToMPISerial();

	/*Free resources:*/
	delete vec_vector;

	/*Assign output pointers:*/
	*pvector=vector;
}/*}}}*/
void GetVectoronBaseFromInputsx(IssmDouble** pvector,FemModel* femmodel,int name, int type){/*{{{*/

	/*output: */
	IssmDouble* vector=NULL;

	/*intermediary: */
	Vector<IssmDouble>* vec_vector=NULL;

	GetVectoronBaseFromInputsx(&vec_vector,femmodel,name,type);
	vector=vec_vector->ToMPISerial();

	/*Free resources:*/
	delete vec_vector;

	/*Assign output pointers:*/
	*pvector=vector;
}/*}}}*/
void GetVectorFromInputsx(IssmDouble** pvector,int* pvector_size, FemModel* femmodel,int name){ /*{{{*/

	int interpolation_type;
	/*this one is special: we don't specify the type, but let the nature of the inputs dictace.
	 * P0 -> ElementSIdEnum, P1 ->VertexSIdEnum: */

	/*We go find the input of the first element, and query its interpolation type: */
	Element* element=xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(0));
	Input* input=element->GetInput(name);
	if (!input) _error_("could not find input: " << name);

	interpolation_type=input->GetInputInterpolationType();
	if(interpolation_type==P0Enum){
		*pvector_size=femmodel->elements->NumberOfElements();
		GetVectorFromInputsx(pvector,femmodel,name, ElementSIdEnum);
	}
	else if(interpolation_type==P1Enum){
		*pvector_size=femmodel->vertices->NumberOfVertices();
		GetVectorFromInputsx(pvector,femmodel,name, VertexSIdEnum);
	}
	else _error_("interpolation type : " << interpolation_type << " not supported yet!");
}/*}}}*/
