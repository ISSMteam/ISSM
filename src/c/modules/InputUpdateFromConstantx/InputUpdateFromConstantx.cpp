/*!\file InputUpdateFromConstantx
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromConstantx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/Inputs.h"

void InputUpdateFromConstantx(FemModel* femmodel,bool constant, int name){
	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->InputUpdateFromConstant(constant,name);
	}
}
void InputUpdateFromConstantx(FemModel* femmodel,int constant, int name){

	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(Object* & object : femmodel->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->InputUpdateFromConstant(constant,name);
	}
}
void InputUpdateFromConstantx(FemModel* femmodel,int constant, int name, int type){

	if(type==P0Enum) InputUpdateFromConstantx(femmodel, constant,name);
	else if(type==P1Enum){

		if(VerboseModule()) _printf0_("   Input updates from constant (P1 version)\n");

		/*Elements and loads drive the update: */
		if(IsInputEnum(name)){
			for(Object* & object : femmodel->elements->objects){
				Element* element = xDynamicCast<Element*>(object);
				element->InputUpdateFromConstant(constant,name,P1Enum);
			}
		}
		else{
			_error_("not supported yet");
		}
	}
	else _error_("InputUpdateFromConstantx error message: type not supported yet!");

}

void InputUpdateFromConstantx(FemModel* femmodel,IssmDouble constant, int name){

	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	if(IsInputEnum(name)){
		for(Object* & object : femmodel->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->InputUpdateFromConstant(constant,name);
		}
	}
	else if(IsParamEnum(name)){
		if(femmodel->parameters->Exist(name)){
			femmodel->parameters->SetParam(constant,name);
		}
		else{
			_error_("Param not set");
		}
	}
	else{
		_error_("not supported");
	}
}
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmDouble constant, int name){

	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->SetElementInput(inputs,name,constant);
	}
}

void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmDouble constant, int name,int type){

	if(type==P0Enum) InputUpdateFromConstantx(inputs, elements, constant,name);
	else if(type==P1Enum){

		if(VerboseModule()) _printf0_("   Input updates from constant (P1 version)\n");

		/*Elements and loads drive the update: */
		if(IsInputEnum(name)){
			for(Object* & object : elements->objects){
				Element* element = xDynamicCast<Element*>(object);
				element->InputUpdateFromConstant(constant,name,P1Enum);
			}
		}
		else{
			_error_("not supported yet");
		}
	}
	else _error_("InputUpdateFromConstantx error message: type not supported yet!");

}
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,bool constant, int name){

	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Elements and loads drive the update: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->SetBoolInput(inputs,name,constant);
	}
}
#ifdef _HAVE_AD_
void InputUpdateFromConstantx(Inputs* inputs,Elements* elements,IssmPDouble constant, int name){

	if(VerboseModule()) _printf0_("   Input updates from constant\n");

	/*Convert to active variable!*/
	IssmDouble constant2 = constant;

	/*Elements and loads drive the update: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->SetElementInput(inputs,name,constant2);
	}
}
#endif
