/*!\file UpdateMmesx
 * \brief: update Mmes prior to running any core
 */

#include "./UpdateMmesx.h"
#include "../MmeToInputFromIdx/MmeToInputFromIdx.h"
#include "../MmeToInputx/MmeToInputx.h"
#include "../InputUpdateFromDakotax/InputUpdateFromDakotax.h"

void  UpdateMmesx(FemModel* femmodel){ 
	

	if(VerboseModule()) _printf0_("   Updating Mmes\n");

	/*Deal with solid earth external: {{{*/

	if (femmodel->inputs->Exist(SolidearthExternalDisplacementUpRateEnum) && femmodel->inputs->GetInputObjectEnum(SolidearthExternalDisplacementUpRateEnum)==DatasetInputEnum){

		int horiz=0;
		femmodel->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
		
		int modelid=0;
		femmodel->parameters->FindParam(&modelid,SolidearthExternalModelidEnum);
		//_printf_("modelid: " << modelid << "\n");
		
		/*replace dataset of forcings with only one, the modelid'th:*/
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementUpRateEnum, P1Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalGeoidRateEnum, P1Enum);

		if (horiz){
			MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementNorthRateEnum, P1Enum);
			MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementEastRateEnum, P1Enum);
		}
		
	} /*}}}*/
	/*Deal with dsl: {{{*/
	if (femmodel->inputs->Exist(OceantransportSpcdslEnum) && femmodel->inputs->GetInputObjectEnum(OceantransportSpcdslEnum)==DatasetInputEnum){
		
		int modelid;
		femmodel->parameters->FindParam(&modelid,DslModelidEnum);

		/*replace dataset of forcings with only one, the modelid'th:*/
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,OceantransportSpcdslEnum, P1Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,OceantransportSpcstrEnum, P0Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,OceantransportSpcbottompressureEnum, P1Enum);
	} /*}}}*/
	/*Deal with solid earth ice loads: {{{*/
	if (femmodel->inputs->Exist(MmemasstransportThicknessEnum) && femmodel->inputs->GetInputObjectEnum(MmemasstransportThicknessEnum)==DatasetInputEnum){ 
	

		int nids;
		IssmDouble* modelids=NULL; 
		IssmDouble* partition = NULL;

		/*retrieve partition vector and model ids necessary to resolve our thickness:*/
		femmodel->parameters->FindParam(&modelids,&nids,NULL,MmemasstransportModelidsEnum);
		femmodel->parameters->FindParam(&partition,NULL,NULL,MmemasstransportPartitionEnum);
		//_printf_("modelids: " << nids << "|ids:" << modelids[0] << "|" << modelids[1] << "|" << modelids[2] << "|\n");

		MmeToInputx(femmodel,modelids,partition,nids,MmemasstransportThicknessEnum, P0Enum);

		if (femmodel->inputs->GetInputObjectEnum(MmemasstransportMaskIceLevelsetEnum)==DatasetInputEnum){
			MmeToInputx(femmodel,modelids,partition,nids,MmemasstransportMaskIceLevelsetEnum, P1Enum);
			/*delete MaskIceLevelsetEnum which will be replaced by MmemasstransportMaskIceLevelsetEnum for each time step: */
			femmodel->inputs->DeleteInput(MaskIceLevelsetEnum);
		}

		if (femmodel->inputs->GetInputObjectEnum(MmemasstransportMaskOceanLevelsetEnum)==DatasetInputEnum){
			MmeToInputx(femmodel,modelids,partition,nids,MmemasstransportMaskOceanLevelsetEnum, P1Enum);
			/*delete MaskOceanLevelsetEnum which will be replaced by MmemasstransportMaskOceanLevelsetEnum for each time step: */
			femmodel->inputs->DeleteInput(MaskOceanLevelsetEnum);
		}

		/*free ressources:*/
		xDelete<IssmDouble>(modelids);
		xDelete<IssmDouble>(partition);


	} /*}}}*/
} 
	
