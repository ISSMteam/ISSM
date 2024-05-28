/*!\file SetControlInputsFromVectorx
 * \brief retrieve vector from inputs in elements
 */

#include "./SetControlInputsFromVectorx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void SetControlInputsFromVectorx(FemModel* femmodel,IssmDouble* vector){

	int  num_controls;
	int* control_type = NULL;
	int* M = NULL;
	int* N = NULL;

	/*Retrieve some parameters*/
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);

	int offset = 0;
	for(int i=0;i<num_controls;i++){
		/*Is the control a Param?*/
		if(IsParamEnum(control_type[i])){
			femmodel->parameters->SetControlFromVector(vector,control_type[i],M[i],N[i],offset);
		}
		else if(IsInputEnum(control_type[i])){
			for(Object* & object : femmodel->elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->SetControlInputsFromVector(vector,control_type[i],i,offset,M[i],N[i]);
			}
		}
		else{
			_error_("not supported yet");
		}
		offset += M[i]*N[i]; 
	}

	xDelete<int>(control_type);
	xDelete<int>(M);
	xDelete<int>(N);
}

void SetControlInputsFromVectorx(FemModel* femmodel,Vector<IssmDouble>* vector){

	IssmDouble* serial_vector=vector->ToMPISerial();
	SetControlInputsFromVectorx(femmodel,serial_vector);
	xDelete<IssmDouble>(serial_vector);
}
#ifdef _HAVE_AD_
void SetControlInputsFromVectorx(FemModel* femmodel,IssmPDouble* vector){

	/*Get total size and recast*/
	int  num_controls;
	int* M = NULL;
	int* N = NULL;
	femmodel->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	femmodel->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	femmodel->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);

	int size = 0;
	for(int i=0;i<num_controls;i++) size += M[i]*N[i]; 

	IssmDouble* serial_vector=xNew<IssmDouble>(size);
	for(int i=0;i<size;i++) serial_vector[i] = reCast<IssmDouble>(vector[i]);

	SetControlInputsFromVectorx(femmodel,serial_vector);

	xDelete<IssmDouble>(serial_vector);
	xDelete<int>(M);
	xDelete<int>(N);
}
#endif
