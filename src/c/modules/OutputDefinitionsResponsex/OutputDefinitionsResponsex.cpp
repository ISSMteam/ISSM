/*!\file OutputDefinitionsResponsex
 * \brief retrieve vector from inputs in elements
 */

#include "./OutputDefinitionsResponsex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"

int OutputDefinitionsResponsex(IssmDouble* presponse, FemModel* femmodel,const char* output_string){

	/*Ok, go find the output definitions dataset in the parameters, where our responses are hiding: */
	DataSet* output_definitions=((DataSetParam*)femmodel->parameters->FindParamObject(OutputdefinitionEnum))->value;

	/*Now, go through the output definitions, and retrieve the object which corresponds to our requested response, output_string: */
	for(Object* & object : output_definitions->objects){
		Definition* definition=dynamic_cast<Definition*>(object);

		char* name = definition->Name();
		if(strcmp(name,output_string)==0){

			/*This is the object that we have been chasing for. compute the response and return: */
			*presponse = definition->Response(femmodel);

			/*cleanup and return*/
			xDelete<char>(name);
			return 0;
		}
		xDelete<char>(name);
	}

	/*If we are here, did not find the definition for this response, not good!: */
	_printf0_("=================================================================\n");
	_printf0_("WARNING: Could not find the output \"" << output_string << "\"\n");
	_printf0_("         - either this output is unavailable for this model run \n");
	_printf0_("         - or there may be a spelling error in your requested_outputs \n");
	_printf0_("         - or there may be a spelling error in an outputdefinition \n");
	_printf0_("           object name (unlikely)\n");
	_printf0_("=================================================================\n");

	return 1;
}

int OutputDefinitionsResponsex(IssmDouble* presponse, FemModel* femmodel,int output_enum){

	/*Ok, go find the output definitions dataset in the parameters, where our responses are hiding: */
	DataSet* output_definitions=((DataSetParam*)femmodel->parameters->FindParamObject(OutputdefinitionEnum))->value;

	/*Now, go through the output definitions, and retrieve the object which corresponds to our requested response, output_enum: */
	for(Object* & object : output_definitions->objects){
		Definition* definition=dynamic_cast<Definition*>(object);

		int en = definition->DefinitionEnum();
		if(en==output_enum){

			/*This is the object that we have been chasing for. compute the response and return: */
			*presponse = definition->Response(femmodel);
			return 0;
		}
	}

	/*If we are here, did not find the definition for this response, not good!: */
	_printf0_("=================================================================\n");
	_printf0_("WARNING: Could not find the output \"" << EnumToStringx(output_enum)<< "\"\n");
	_printf0_("         - either this output is unavailable for this model run \n");
	_printf0_("         - or there may be a spelling error in your requested_outputs \n");
	_printf0_("         - or there may be a spelling error in an outputdefinition \n");
	_printf0_("           object name (unlikely)\n");
	_printf0_("=================================================================\n");

	return 1;
}
