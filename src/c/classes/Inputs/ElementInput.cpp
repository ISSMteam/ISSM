/*!\file ElementInput.c
 * \brief: implementation of the ElementInput object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../../shared/shared.h"
#include "./ElementInput.h"

/*ElementInput constructors and destructor*/
ElementInput::ElementInput(){/*{{{*/
	this->interpolation  = -1;
	this->M              = -1;
	this->N              = -1;
	this->isserved       = false;
	this->element_values = NULL;
	this->values         = NULL;
}
/*}}}*/
ElementInput::~ElementInput(){/*{{{*/
	if(this->element_values) xDelete<IssmDouble>(this->element_values);
	if(this->values)         xDelete<IssmDouble>(this->values);
}
/*}}}*/

/*Numerics*/
int ElementInput::GetInputInterpolationType(void){/*{{{*/

	return this->interpolation;

}/*}}}*/
