/*!\file Cflevelsetmisfit.cpp
 * \brief: Cflevelsetmisfit Object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./classes.h"
#include "./ExternalResults/ExternalResult.h"
#include "./ExternalResults/Results.h"
#include "../datastructures/datastructures.h"
#include "./Elements/Element.h"
#include "./Elements/Elements.h"
#include "./FemModel.h"
#include "../modules/SurfaceAreax/SurfaceAreax.h"
#include "../classes/Params/Parameters.h"
#include "../classes/gauss/Gauss.h"
#include "./Inputs/DatasetInput.h"
/*}}}*/

/*Cflevelsetmisfit constructors, destructors :*/
Cflevelsetmisfit::Cflevelsetmisfit(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->model_enum = UNDEF;
	this->datatime=0.;
	this->timepassedflag = false;
	this->J = 0.;
}
/*}}}*/
Cflevelsetmisfit::Cflevelsetmisfit(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->model_enum=in_model_enum;
	this->datatime=in_datatime;
	this->timepassedflag=false;
	this->J = 0.;
}
/*}}}*/
Cflevelsetmisfit::Cflevelsetmisfit(char* in_name, int in_definitionenum, int in_model_enum, IssmDouble in_datatime, bool in_timepassedflag, IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->model_enum=in_model_enum;
	this->datatime=in_datatime;
	this->timepassedflag=in_timepassedflag;
	this->J = in_J;
}
/*}}}*/
Cflevelsetmisfit::~Cflevelsetmisfit(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Cflevelsetmisfit::copy() {/*{{{*/
	Cflevelsetmisfit* mf = new Cflevelsetmisfit(this->name,this->definitionenum, this->model_enum,this->datatime,this->timepassedflag, this->J);
	return (Object*) mf;
}
/*}}}*/
void Cflevelsetmisfit::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cflevelsetmisfit::Echo(void){/*{{{*/
	_printf_(" Cflevelsetmisfit: " << name << " " << this->definitionenum << "\n");
	_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
	_printf_("    datatime: " << datatime << "\n");
	_printf_("	  timepassedflag: "<<timepassedflag<<"\n");
}
/*}}}*/
int Cflevelsetmisfit::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cflevelsetmisfit::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	/*ok, marshall operations: */
   int object_enum=CflevelsetmisfitEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->model_enum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->datatime);
	marshallhandle->call(this->timepassedflag);
	marshallhandle->call(this->J);
} 
/*}}}*/
int Cflevelsetmisfit::ObjectEnum(void){/*{{{*/
	return CflevelsetmisfitEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Cflevelsetmisfit::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cflevelsetmisfit::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cflevelsetmisfit::Response(FemModel* femmodel){/*{{{*/
	 /*diverse: */
	 IssmDouble time;

	 /*recover time parameters: */
	 femmodel->parameters->FindParam(&time,TimeEnum);
	 if(datatime<=time && !timepassedflag){

		 IssmDouble J_part = 0.;
		 IssmDouble J_sum  = 0.;

		 for(Object* & object : femmodel->elements->objects){
			 Element* element=xDynamicCast<Element*>(object);
			 J_part+=this->Cflevelsetmisfit_Calculation(element,model_enum);
		 }

		 ISSM_MPI_Allreduce ( (void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		 ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		 this->timepassedflag = true;
		 this->J = J_sum;
	 }

	 return this->J;
 }/*}}}*/
IssmDouble Cflevelsetmisfit::Cflevelsetmisfit_Calculation(Element* element, int model_enum){/*{{{*/

	int        domaintype,numcomponents;
	IssmDouble Jelem=0.;
	IssmDouble misfit,Jdet;
	IssmDouble model,obs,weight;
	IssmDouble modelLevel,obsLevel;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnSurface()) return 0.;

	/*If on water, return 0: */
//	if(!element->IsIceInElement()) return 0.;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   numcomponents   = 1; break;
		case Domain3DEnum:           numcomponents   = 2; break;
		case Domain2DhorizontalEnum: numcomponents   = 2; break;
		default: _error_("not supported yet");
	}

	/*Spawn surface element*/
	Element* topelement = element->SpawnTopElement();

	/* Get node coordinates*/
	topelement->GetVerticesCoordinates(&xyz_list);

	/*Retrieve all inputs we will be needing: */
	DatasetInput *datasetinput = topelement->GetDatasetInput(definitionenum); _assert_(datasetinput);
	Input        *model_input  = topelement->GetInput(model_enum);            _assert_(model_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		datasetinput->GetInputValue(&weight,gauss,WeightsLevelsetObservationEnum);
		model_input->GetInputValue(&model,gauss);
		datasetinput->GetInputValue(&obs,gauss,LevelsetObservationEnum);

		/*Compute Levelset misfit:
		 *  J = || H(\phi) - H(\phi_obs)||^2
		 *                                           */
//		modelLevel = this->Heaviside(model);
//		obsLevel = this->Heaviside(obs);
		modelLevel = model;
		obsLevel = obs;

		misfit=0.5*(modelLevel-obsLevel)*(modelLevel-obsLevel);

//		if (modelLevel*obsLevel<0.0) {
		/*Add to cost function*/
			Jelem+=misfit*weight*Jdet*gauss->weight;
//		}
	}

	/*clean up and Return: */
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}/*}}}*/
IssmDouble Cflevelsetmisfit::Heaviside(IssmDouble x){/*{{{*/
	if (x>0) return 1.;
	else if (x<0) return -1.;
	else return 0.;
}
	/*}}}*/
