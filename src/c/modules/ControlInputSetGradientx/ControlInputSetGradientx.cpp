/*!\file ControlInputSetGradientx
 * \brief retrieve gradient from inputs in elements
 */

#include "./ControlInputSetGradientx.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,IssmDouble* gradient){

	/*Intermediaries*/
	int  num_controls;
	int *control_type = NULL;
	int* M_all = NULL;
	int* N_all = NULL;
	int* interp_all = NULL;

	/*Retrieve some parameters*/
	parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	parameters->FindParam(&M_all,NULL,ControlInputSizeMEnum);
	parameters->FindParam(&N_all,NULL,ControlInputSizeNEnum);
	parameters->FindParam(&interp_all,NULL,ControlInputInterpolationEnum);

	int offset = 0;
	for(int i=0;i<num_controls;i++){
		/*Is the control a Param?*/
		if(IsParamEnum(control_type[i])){
			parameters->SetGradientFromVector(gradient, control_type[i], M_all[i], N_all[i], offset);
		}
		else if(IsInputEnum(control_type[i])){
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->ControlInputSetGradient(gradient,control_type[i],i,offset,M_all[i],N_all[i],interp_all[i]);
			}
		}
		else{
			_error_("not supported yet");
		}
		offset+=M_all[i]*N_all[i];
	}

	/*Clean up and return*/
	xDelete<int>(control_type);
	xDelete<int>(M_all);
	xDelete<int>(N_all);
	xDelete<int>(interp_all);

}
void ControlInputSetGradientx(Elements* elements,Nodes* nodes, Vertices* vertices, Loads* loads, Materials* materials, Parameters* parameters,Vector<IssmDouble>* gradient){

	/*Serialize gradient*/
	IssmDouble* serial_gradient=gradient->ToMPISerial();

	ControlInputSetGradientx(elements,nodes,vertices, loads, materials, parameters,serial_gradient);

	/*Clean up and return*/
	xDelete<IssmDouble>(serial_gradient);
}
