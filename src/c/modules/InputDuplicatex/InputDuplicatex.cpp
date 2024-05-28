/*!\file InputDuplicatex
 * \brief: duplicte  an input inside the elements, onto another, and wipe it off.
 */

#include "./InputDuplicatex.h"
#include "../../shared/shared.h"
#include "../../classes/classes.h"
#include "../../toolkits/toolkits.h"

void InputDuplicatex(FemModel* femmodel,int original_enum, int new_enum){
	femmodel->inputs->DuplicateInput(original_enum,new_enum);
}
