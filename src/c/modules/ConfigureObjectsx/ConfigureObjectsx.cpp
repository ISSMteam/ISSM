/*!\file ConfigureObjectsx
 * \brief: configure objects in elements and loads to link in with nodes
 */

#include "./ConfigureObjectsx.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"

int	ConfigureObjectsx( Elements* elements, Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters,Inputs* inputs){

	if(VerboseMProcessor()) _printf0_("      Configuring elements...\n");
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->Configure(elements,loads,nodes,vertices,materials,parameters,inputs);
	}
	if(VerboseMProcessor()) _printf0_("      Configuring loads...\n");
	for(Object* object : loads->objects){
		Load* load=(Load*)object;
		load->Configure(elements,loads,nodes,vertices,materials,parameters);
	}
	if(VerboseMProcessor()) _printf0_("      Configuring materials...\n");
	for(Object* & object : materials->objects){
		Material* material=(Material*)object;
		material->Configure(elements);
	}
	if(VerboseMProcessor()) _printf0_("      Configuring inputs...\n");
	inputs->Configure(parameters);

	return 1;
}
