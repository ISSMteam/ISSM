/*
 * \file Materials.cpp
 * \brief: Implementation of Materials class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Materials.h"
#include "./Material.h"
#include "../../shared/shared.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Materials::Materials(){/*{{{*/
	enum_type=MaterialsEnum;
	return;
}
/*}}}*/
Materials::~Materials(){/*{{{*/
	return;
}
/*}}}*/

/*Object management*/
void Materials::Configure(Elements* elements,Loads* loads, Nodes* nodes, Vertices* vertices, Materials* materials,Parameters* parameters){/*{{{*/

	vector<Object*>::iterator object;
	Material* material=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		material=xDynamicCast<Material*>(*object);
		material->Configure(elements);

	}

}
/*}}}*/
void Materials::ResetHooks(){/*{{{*/

	vector<Object*>::iterator object;
	Material* material=NULL;

	for ( object=objects.begin() ; object < objects.end(); object++ ){

		material=xDynamicCast<Material*>((*object));
		material->ResetHooks();

	}

}
/*}}}*/
