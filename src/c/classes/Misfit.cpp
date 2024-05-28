/*!\file Misfit.cpp
 * \brief: Misfit object
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
#include "../modules/OutputDefinitionsResponsex/OutputDefinitionsResponsex.h"
#include "../modules/GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "../classes/Params/Parameters.h"
#include "../classes/gauss/Gauss.h"
/*}}}*/

/*Misfit constructors, destructors :*/
Misfit::Misfit(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->model_enum = UNDEF;
	this->observation_enum = UNDEF;
	this->weights_enum = UNDEF;
	this->timeinterpolation=NULL;
	this->local=1;
	this->misfit=0;
	this->lock=0;

}
/*}}}*/
Misfit::Misfit(char* in_name, int in_definitionenum, int in_model_enum, int in_observation_enum, char* in_timeinterpolation, int in_local, int in_weights_enum){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->timeinterpolation = xNew<char>(strlen(in_timeinterpolation)+1);
	xMemCpy<char>(this->timeinterpolation,in_timeinterpolation,strlen(in_timeinterpolation)+1);

	this->model_enum=in_model_enum;
	this->observation_enum=in_observation_enum;
	this->weights_enum=in_weights_enum;
	this->local=in_local;

	this->misfit=0;
	this->lock=0;
}
/*}}}*/
Misfit::~Misfit(){/*{{{*/
	if(this->name)xDelete(this->name);
	if(this->timeinterpolation)xDelete(this->timeinterpolation);
	this->misfit=0;
	this->lock=0;
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Misfit::copy() {/*{{{*/
	Misfit* mf = new Misfit(this->name,this->definitionenum, this->model_enum,this->observation_enum,this->timeinterpolation,this->local,this->weights_enum);
	mf->misfit=this->misfit;
	mf->lock=this->lock;
	return (Object*) mf;
}
/*}}}*/
void Misfit::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Misfit::Echo(void){/*{{{*/
	_printf_(" Misfit: " << name << " " << this->definitionenum << "\n");
	_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
	_printf_("    observation_enum: " << observation_enum << " " << EnumToStringx(observation_enum) << "\n");
	_printf_("    weights_enum: " << weights_enum << " " << EnumToStringx(weights_enum) << "\n");
	_printf_("    timeinterpolation: " << timeinterpolation << "\n");
	_printf_("    local: " << local << "\n");
}
/*}}}*/
int Misfit::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Misfit::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	_error_("not implemented yet!");
}
/*}}}*/
int Misfit::ObjectEnum(void){/*{{{*/
	return MisfitEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Misfit::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Misfit::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Misfit::Response(FemModel* femmodel){/*{{{*/

	 /*diverse: */
	 IssmDouble time,starttime,finaltime;
	 IssmDouble dt;

	 /*recover time parameters: */
	 femmodel->parameters->FindParam(&starttime,TimesteppingStartTimeEnum);
	 femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	 femmodel->parameters->FindParam(&time,TimeEnum);
	 femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);

	 if (this->local==1){ /*area integration using elements: {{{*/

		 IssmDouble misfit_t=0.;
		 IssmDouble all_misfit_t=0.;
		 IssmDouble area_t=0.;
		 IssmDouble all_area_t;

		 /*If we are locked, return time average: */
		 if(this->lock)return misfit/(time-starttime);

		 for(Object* & object : femmodel->elements->objects){
			 Element* element=xDynamicCast<Element*>(object);
			 misfit_t+=element->Misfit(model_enum,observation_enum,weights_enum);
			 area_t+=element->MisfitArea(weights_enum);
		 }

		 ISSM_MPI_Allreduce ( (void*)&misfit_t,(void*)&all_misfit_t,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		 ISSM_MPI_Allreduce ( (void*)&area_t,(void*)&all_area_t,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		 area_t=all_area_t;
		 misfit_t=all_misfit_t;

		 /*Divide by surface area if not nill!: */
		 if (area_t!=0) misfit_t=misfit_t/area_t;

		 /*Add this time's contribution to curent misfit: */
		 misfit+=dt*misfit_t;

		 /*Do we lock? i.e. are we at final_time? :*/
		 if(time==finaltime)this->lock=1;

		 /*What we return is the value of misfit / time if transient*/
		 if(time!=0.) return misfit/(time-starttime);
		 return misfit;
	 } /*}}}*/
	 else if (this->local==2){ /*vertex by vertex computation: {{{*/

		 IssmDouble* model = NULL;
		 IssmDouble* observation= NULL;
		 IssmDouble* weights= NULL;
		 int msize,osize,wsize;

		 /*Are we transient?:*/
		 if (time==0){
			 IssmDouble misfit_t=0.;

			 /*get global vectors: */
			 GetVectorFromInputsx(&model,&msize,femmodel,model_enum);
			 GetVectorFromInputsx(&observation,&osize,femmodel,observation_enum);_assert_(msize==osize);
			 GetVectorFromInputsx(&weights,&wsize,femmodel,weights_enum); _assert_(wsize==msize);

			 int count=0;
			 for (int i=0;i<msize;i++){
				 misfit_t += pow(model[i]-observation[i],2)*weights[i];
				 if (weights[i]!=0)count++;
			 }
			 misfit=sqrt(misfit_t/count);

			 /*Free resources:*/
			 xDelete<IssmDouble>(model);
			 xDelete<IssmDouble>(observation);
			 xDelete<IssmDouble>(weights);

			 /*return value: */
			 return misfit;
		 }
		 else{

			 IssmDouble misfit_t=0.;
			 IssmDouble all_misfit_t=0.;

			 /*If we are locked, return time average: */
			 if(this->lock)return misfit/(time-starttime);

			 /*get global vectors: */
			 GetVectorFromInputsx(&model,&msize,femmodel,model_enum);
			 GetVectorFromInputsx(&observation,&osize,femmodel,observation_enum);_assert_(msize==osize);
			 GetVectorFromInputsx(&weights,&wsize,femmodel,weights_enum); _assert_(wsize==msize);

			 int count=0;
			 for (int i=0;i<msize;i++){
				 misfit_t += pow(model[i]-observation[i],2)*weights[i];
				 if (weights[i]!=0)count++;
			 }
			 misfit=sqrt(misfit_t/count);

			 /*Add this time's contribution to curent misfit: */
			 misfit=sqrt(misfit_t)/count;
			 misfit+=dt*misfit_t;

			 /*Do we lock? i.e. are we at final_time? :*/
			 if(time==finaltime)this->lock=1;

			 /*Free resources:*/
			 xDelete<IssmDouble>(model);
			 xDelete<IssmDouble>(observation);
			 xDelete<IssmDouble>(weights);

			 /*What we return is the value of misfit / time: */
			 return misfit/(time-starttime);
		 }

	 } /*}}}*/
	 else{ /*global computation: {{{ */

		 IssmDouble model, observation;
		 int ierr;

		 /*If we are locked, return time average: */
		 if(this->lock) return misfit/(time-starttime);

		 /*First, the global  model response: */
		 ierr = OutputDefinitionsResponsex(&model, femmodel, this->model_enum);
		 if(ierr) _error_("could not evaluate response");

		 /*Now, the observation is buried inside the elements, go fish it in the first element (cludgy, needs fixing): */
		 Element* element = (Element*)femmodel->elements->GetObjectByOffset(0); _assert_(element);
		 Input*  input   = element->GetInput(observation_enum); _assert_(input);
		 input->GetInputAverage(&observation);

		 /*Add this time's contribution to curent misfit: */
		 misfit+=dt*(model-observation);

		 /*Do we lock? i.e. are we at final_time? :*/
		 if(time==finaltime)this->lock=1;

		 /*What we return is the value of misfit / time: */
		 return misfit/(time-starttime);
	 } /*}}}*/

 }
	/*}}}*/
