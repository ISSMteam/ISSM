/*!\file InputUpdateFromDakotax
 * \brief: update datasets using  parameter inputs
 */

#include "./InputUpdateFromDakotax.h"
#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/Inputs/TransientInput.h"
#include "../../classes/Inputs/DatasetInput.h"
#include "../../classes/Inputs/TriaInput.h"
#include "../InputUpdateFromMatrixDakotax/InputUpdateFromMatrixDakotax.h"
#include "../InputUpdateFromConstantx/InputUpdateFromConstantx.h"
#include "../InputUpdateFromVectorDakotax/InputUpdateFromVectorDakotax.h"
#include "../UpdateMmesx/UpdateMmesx.h"
#include "../MmeToInputx/MmeToInputx.h"

void  InputUpdateFromDakotax(FemModel* femmodel,double* variables,char* *variables_descriptors,int numdakotavariables){ /*{{{*/

	int     i,j,k,l;

	IssmDouble **variable_partitions         = NULL;
	IssmDouble * variable_partition         = NULL;
	int * variable_partitions_npart         = NULL;
	int * variable_partitions_nt         = NULL;
	int          variable_partitions_num;
	int          npart;
	int          nt;
	int variablecount=0;

	double *distributed_values = NULL;
	double *parameter          = NULL;
	char   *descriptor         = NULL;
	char    root[50]; //root name of variable, ex: DragCoefficent, RhoIce, etc ...

	if (VerboseQmu())_printf0_("dakota variables updates\n");

	/*retrieve parameters: */
	femmodel->parameters->FindParam(&variable_partitions,&variable_partitions_num,NULL,NULL,QmuVariablePartitionsEnum);
	femmodel->parameters->FindParam(&variable_partitions_npart,NULL,NULL,QmuVariablePartitionsNpartEnum);
	femmodel->parameters->FindParam(&variable_partitions_nt,NULL,NULL,QmuVariablePartitionsNtEnum);

	/*Go through all dakota descriptors, ex: "rho_ice","thermal_conductivity","thickness1","thickness2", etc ..., and
	 * for each descriptor, take the variable value and plug it into the inputs (more or less :)):
	 * We also start with distributed and standard values , as they tend to be used to pluck data from a multi-modle ensemble (mme)
	 * which can then be scaled. Doing the scaling first would be impractical, as the entire mme would have to be scaled,
	 * which is a waste of time:*/

	variablecount=0;
	for(i=0;i<numdakotavariables;i++){ //these are the dakota variables, for all partitions.

		descriptor=variables_descriptors[i];

		/*From descriptor, figure out if the variable is scaled, indexed, distributed or just a simple variable: */
		if (strncmp(descriptor,"scaled_",7)==0){
			/*we are skipping these for now.*/
			npart=variable_partitions_npart[variablecount];
			nt=variable_partitions_nt[variablecount];

			/*increment i to skip the distributed values just collected: */
			i+=npart*nt-1; //careful, the for loop will add 1.
		}
		else if (strncmp(descriptor,"indexed_",8)==0){
			/*we are skipping these for now.*/
		}
		else if (strncmp(descriptor,"nodal_",8)==0){
			/*we are skipping these for now.*/
		}

		else if (strncmp(descriptor,"distributed_",12)==0){
			if (VerboseQmu())_printf0_("   updating variable " << descriptor << "\n");

			/*recover partition vector: */
			variable_partition=variable_partitions[variablecount];
			npart=variable_partitions_npart[variablecount];

			/*Variable is distributed. Determine root name of variable (ex: distributed_DragCoefficient_1 -> DragCoefficient).
			 * Allocate distributed_values and fill the distributed_values with the next npart variables: */

			memcpy(root,strstr(descriptor,"_")+1,(strlen(strstr(descriptor,"_")+1)+1)*sizeof(char));
			*strstr(root,"_")='\0';

			distributed_values=xNew<double>(npart);
			for(j=0;j<npart;j++){
				distributed_values[j]=variables[i+j];
			}

			//for (int j=0;j<npart;j++)_printf_(j << ":" << distributed_values[j] << "\n");

			//Call specialty code:
			InputUpdateSpecialtyCode(femmodel,distributed_values,variable_partition,npart,root);

			/*increment i to skip the distributed values just collected: */
			i+=npart-1; //careful, the for loop will add 1.

			/*Free allocations: */
			xDelete<double>(parameter);
			xDelete<double>(distributed_values);
		}
		else{
			/*Ok, standard variable, just update inputs using the variable: */
			if (VerboseQmu())_printf0_("   updating variable " << descriptor << "\n");
			InputUpdateFromConstantx(femmodel,variables[i],StringToEnumx(descriptor));
		}
		variablecount++;
	}

	variablecount=0;
	/*now deal with scaled variabes:*/
	for(i=0;i<numdakotavariables;i++){ //these are the dakota variables, for all partitions.

		descriptor=variables_descriptors[i];

		/*From descriptor, figure out if the variable is scaled, indexed, distributed or just a simple variable: */
		if (strncmp(descriptor,"scaled_",7)==0){

			if (VerboseQmu())_printf0_("   updating variable " << descriptor << "\n");

			/*recover partition vector: */
			variable_partition=variable_partitions[variablecount];
			npart=variable_partitions_npart[variablecount];
			nt=variable_partitions_nt[variablecount];

			/* Variable is scaled, determine its root name (ex: scaled_DragCoefficient_1 -> DragCoefficient). Allocate distributed_values and fill the
			 * distributed_values with the next npart variables coming from Dakota: */
			memcpy(root,strstr(descriptor,"_")+1,(strlen(strstr(descriptor,"_")+1)+1)*sizeof(char));
			*strstr(root,"_")='\0';

			distributed_values=xNew<double>(npart*nt);
			for(j=0;j<npart*nt;j++){
				distributed_values[j]=variables[i+j];
			}

			/*Scale variable inside the inputs:*/
			InputScaleFromDakotax(femmodel, distributed_values, variable_partition,npart, nt, StringToEnumx(root));

			/*increment i to skip the distributed values just collected: */
			i+=npart*nt-1; //careful, the for loop will add 1.

			/*Free allocations: */
			xDelete<double>(parameter);
			xDelete<double>(distributed_values);
		}
		variablecount++;
	}

	/*Resole Mmes now that we have updated model ids: */
	UpdateMmesx(femmodel);

	/*Save results:*/
	femmodel->results->AddResult(new GenericExternalResult<IssmPDouble*>(femmodel->results->Size()+1,"uq_variables",variables,numdakotavariables,1,1,0));

	/*Free resources:*/
	for(i=0;i<variable_partitions_num;i++){
		IssmDouble* matrix=variable_partitions[i];
		xDelete<IssmDouble>(matrix);
	}
	xDelete<IssmDouble*>(variable_partitions);
	xDelete<int>(variable_partitions_npart);
	xDelete<int>(variable_partitions_nt);

} /*}}}*/
void  InputUpdateSpecialtyCode(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* variable_partition,int npart,char* root){ //{{{

	/*Here, we put all the code that cannot be handled any other place: */
	if (strncmp(root,"MmemasstransportThickness",25)==0){ //surface load in solid earth class {{{

		if(VerboseQmu()){
			_printf0_("      Masstransport Thickness MME with ids: ");
			for (int i=0;i<npart;i++)_printf0_((int)distributed_values[i] << " ");
			_printf0_("\n");
		}

		if (femmodel->inputs->GetInputObjectEnum(MmemasstransportThicknessEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MmemasstransportThicknessEnum, P0Enum);

		if (femmodel->inputs->GetInputObjectEnum(MmemasstransportMaskIceLevelsetEnum)==DatasetInputEnum){
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MmemasstransportMaskIceLevelsetEnum, P1Enum);
			/*delete MaskIceLevelsetEnum which will be replaced by MmemasstransportMaskIceLevelsetEnum for each time step: */
			femmodel->inputs->DeleteInput(MaskIceLevelsetEnum);
		}

		if (femmodel->inputs->GetInputObjectEnum(MmemasstransportMaskOceanLevelsetEnum)==DatasetInputEnum){
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MmemasstransportMaskOceanLevelsetEnum, P1Enum);
			
			/*delete MaskOceanLevelsetEnum which will be replaced by MmemasstransportMaskOceanLevelsetEnum for each time step: */
			femmodel->inputs->DeleteInput(MaskOceanLevelsetEnum);
		}

	} /*}}}*/ 
	else if (strncmp(root,"SurfaceloadModelid",18)==0){ //surface load in solid earth class {{{

		if(VerboseQmu()){
			_printf0_("      SurfaceloadModelid MME, with ids: ");
			for (int i=0;i<npart;i++)_printf0_((int)distributed_values[i] << " ");
			_printf0_("\n");
		}

		if (femmodel->inputs->GetInputObjectEnum(MasstransportSpcthicknessEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MasstransportSpcthicknessEnum, P0Enum);

		if (femmodel->inputs->GetInputObjectEnum(MaskIceLevelsetEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MaskIceLevelsetEnum, P1Enum);

		if (femmodel->inputs->GetInputObjectEnum(MaskOceanLevelsetEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,MaskOceanLevelsetEnum, P1Enum);

	} /*}}}*/
	else if (strncmp(root,"SolidearthExternalModelid",18)==0){ //external solid earth solution in solid earth class {{{

		if(VerboseQmu()){
			_printf0_("      SolidearthExternalModelid MME, with ids: ");
			for (int i=0;i<npart;i++)_printf0_((int)distributed_values[i] << " ");
			_printf0_("\n");
		}

		if (femmodel->inputs->GetInputObjectEnum(SolidearthExternalDisplacementEastRateEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,SolidearthExternalDisplacementEastRateEnum, P1Enum);

		if (femmodel->inputs->GetInputObjectEnum(SolidearthExternalDisplacementUpRateEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,SolidearthExternalDisplacementUpRateEnum, P1Enum);

		if (femmodel->inputs->GetInputObjectEnum(SolidearthExternalDisplacementNorthRateEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,SolidearthExternalDisplacementNorthRateEnum, P1Enum);

		if (femmodel->inputs->GetInputObjectEnum(SolidearthExternalGeoidRateEnum)==DatasetInputEnum)
			MmeToInputx(femmodel,distributed_values,variable_partition,npart,SolidearthExternalGeoidRateEnum, P1Enum);

		//if (femmodel->inputs->GetInputObjectEnum(SolidearthExternalBarystaticSeaLevelRateEnum)==DatasetInputEnum)
		//	MmeToInput(femmodel,distributed_values,variable_partition,npart,SolidearthExternalBarystaticSeaLevelRateEnum, P1Enum);
	} /*}}}*/
	else _error_("InputUpdateSpecialtyCode error message: " << root << " not supported yet!");

}	//}}}
void  InputScaleFromDakotax(FemModel* femmodel,IssmDouble* distributed_values,IssmDouble* partition, int npart, int nt, int name){ /*{{{*/

	/*Copy input:*/
	femmodel->inputs->DuplicateInput(name,DummyEnum);

	/*Go through elements, copy input name to dummy, and scale it using the distributed_values and the partition vector:*/
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->InputScaleFromDakota(distributed_values,partition,npart,nt,name);
	}

	/*We created a dummy input, which was a scaled copy of the name input. Now wipe
	 * out the name input with the new input:*/
	femmodel->inputs->ChangeEnum(DummyEnum,name);

	/*Some specialty code:*/
	switch(name){
		case MaterialsRheologyBEnum:
			femmodel->inputs->DuplicateInput(name,MaterialsRheologyBbarEnum);
			break;
	}

} /*}}}*/
