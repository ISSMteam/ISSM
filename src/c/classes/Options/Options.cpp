/*
 * \file Options.cpp
 * \brief: Implementation of Options class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <vector>
#include <algorithm>
#include <cstring>

#include "./Options.h"
#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"
/*}}}*/

/*Object constructors and destructor*/
Options::Options(){/*{{{*/
	return;
}
/*}}}*/
Options::~Options(){/*{{{*/
	return;
}
/*}}}*/

/*Object management*/
int  Options::AddOption(Option* in_option){/*{{{*/

	char* name=NULL;

	vector<Object*>::iterator object;
	Option* option=NULL;

	/*In debugging mode, check that the option is not a NULL pointer*/
	_assert_(in_option);

	/*Also, check the option name*/
	name=in_option->Name();

	if(!name) _error_("input option has an empty name");
	if(strchr(name,'.')) _error_("Option \"" << name << "\" has a protected character \".\"");
	if(strchr(name,'[')) _error_("Option \"" << name << "\" has a protected character \"[\"");
	if(strchr(name,']')) _error_("Option \"" << name << "\" has a protected character \"]\"");

	/*Finally, check that no option of the same name already exists in the dataset*/
	for(object=objects.begin();object<objects.end();object++){

		option=xDynamicCast<Option*>(*object);
		if (!strcmp(option->Name(),name)){
			_error_("Options \"" << name << "\" found multiple times");
			break;
		}
	}

	/*OK, all checks went well, add option to dataset*/
	this->AddObject(in_option);

	return 1;
}
/*}}}*/
Option* Options::GetOption(const char* name){/*{{{*/

	vector<Object*>::iterator object;
	Option* option=NULL;

	/*Go through options and find option: */
	for ( object=objects.begin() ; object < objects.end(); object++ ){

		option=xDynamicCast<Option*>(*object);
		//option=(Option*)(*object); //C-like cast
		/*There is a crash on some machines (Einar Olason) that needs to be fixed*/
		if(!option){
			_printf_("The dynamic_cast from Object* to Option* is failing.\n");
			_printf_("\n");
			_printf_("A quick workaround consists of using a C-like cast\n");
			_printf_("\n");
			_printf_("Open Options.cpp and change the dynamic_cast in Options::GetOption by a C-like cast\n");
			//_printf_("Open Options.h and replace the dynamic_cast of all the Get functions to C-like cats\n");
			_printf_("\n");
			_error_("Make the fix above and recompile ISSM");
		}

		if (!strncmp(name,option->Name(),strlen(option->Name()))){

			/*OK, now do we have a complete name? If not, it is a cell or a structure, we need to go further*/
			if(!strcmp(name,option->Name())){
				return option;
			}
			else{
				_error_("Cannot recover field \"" << name << "\" for an option of type " << EnumToStringx(option->ObjectEnum()));
			}
		}
	}

	/*Option not found return NULL pointer*/
	return NULL;
}
/*}}}*/
