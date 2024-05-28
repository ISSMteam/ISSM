/*!\file Damagex
 * \brief: compute damage
 */

#include "./Damagex.h"
#include "../InputDuplicatex/InputDuplicatex.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

void Damagex(FemModel* femmodel){

	/*Recover Damage law Enum*/
	int damagelaw;
	femmodel->parameters->FindParam(&damagelaw,DamageLawEnum);

	/*Calculate damage*/
	switch(damagelaw){
		case 0:
			if(VerboseModule()) _printf0_("   computing damage analytically\n");
			InputDuplicatex(femmodel,DamageDEnum,DamageDOldEnum);
			InputDuplicatex(femmodel,DamageDbarEnum,DamageDbarOldEnum);
			femmodel->ElementOperationx(&Element::ComputeNewDamage);
			break;
		case 1:
		case 2:
			if(VerboseModule()) _printf0_("   computing damage using source term in advection scheme\n");
			/* Damage calculated using source term in DamageEvolutionAnalysis */
			break;
		case 3:
			if(VerboseModule()) _printf0_("   computing damage using source term in advection scheme\n");
			/* Damage calculated using source term in DamageEvolutionAnalysis */
			break;
		default:
			_error_("Damage law "<<EnumToStringx(damagelaw)<<" not implemented yet");
	}
}
