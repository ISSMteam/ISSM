/*
 * \file Constraints.cpp
 * \brief: Implementation of Constraints class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Constraints.h"
#include "./Constraint.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

using namespace std;
/*}}}*/

/*Numerics: */
void Constraints::ActivatePenaltyMethod(int in_analysis){/*{{{*/

	for(Object* & object: this->objects){
		Constraint* constraint=(Constraint*)object;
		constraint->ActivatePenaltyMethod();
	}

}
/*}}}*/
