/*!\file Cfsurfacesquaretransient.cpp
 * \brief: Cfsurfacesquaretransient Object
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

/*Cfsurfacesquaretransient constructors, destructors :*/
Cfsurfacesquaretransient::Cfsurfacesquaretransient(){/*{{{*/

	this->definitionenum = -1;
	this->name           = NULL;
	this->model_enum     = UNDEF;
	this->datatimes      = NULL;
	this->passedflags    = NULL;
	this->J              = 0.;
}
/*}}}*/
Cfsurfacesquaretransient::Cfsurfacesquaretransient(char* in_name, int in_definitionenum, int in_model_enum, int in_num_datatimes, IssmDouble* in_datatimes, bool* in_passedflags, IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name = xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->model_enum=in_model_enum;
	this->num_datatimes = in_num_datatimes;

	/*Allocate arrays*/
	_assert_(this->num_datatimes>0);
	this->datatimes   = xNew<IssmDouble>(this->num_datatimes);
	this->passedflags = xNew<bool>(this->num_datatimes);
	xMemCpy<IssmDouble>(this->datatimes,in_datatimes,this->num_datatimes);
	xMemCpy<bool>(this->passedflags,in_passedflags,this->num_datatimes);

	#ifdef _ISSM_DEBUG_ 
	for(int i=0;i<this->num_datatimes-1;i++){
		if(this->datatimes[i+1]<=this->datatimes[i]){
			_error_("time series is not in chronological order");
		}
	}
	#endif

	this->J = in_J;
}
/*}}}*/
Cfsurfacesquaretransient::Cfsurfacesquaretransient(char* in_name, int in_definitionenum, int in_model_enum, int in_num_datatimes, IssmDouble* in_datatimes){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name = xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->model_enum=in_model_enum;
	this->num_datatimes = in_num_datatimes;

	/*Allocate arrays*/
	_assert_(this->num_datatimes>0);
	this->datatimes   = xNew<IssmDouble>(this->num_datatimes);
	this->passedflags = xNew<bool>(this->num_datatimes);
	xMemCpy<IssmDouble>(this->datatimes,in_datatimes,this->num_datatimes);

	/*initialize passedtimes to false*/
	for(int i=0;i<this->num_datatimes;i++) this->passedflags[i]= false;
	this->J = 0;
}
/*}}}*/
Cfsurfacesquaretransient::~Cfsurfacesquaretransient(){/*{{{*/
	if(this->name) xDelete(this->name);
	if(this->datatimes) xDelete(this->datatimes);
	if(this->passedflags) xDelete(this->passedflags);
}
/*}}}*/

/*Object virtual function resolutoin: */
Object* Cfsurfacesquaretransient::copy() {/*{{{*/
	Cfsurfacesquaretransient* output = new Cfsurfacesquaretransient(this->name,this->definitionenum, this->model_enum, this->num_datatimes, this->datatimes,this->passedflags, this->J);
	return (Object*)output;
}
/*}}}*/
void Cfsurfacesquaretransient::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cfsurfacesquaretransient::Echo(void){/*{{{*/
	_printf_(" Cfsurfacesquaretransient: " << name << " " << this->definitionenum << "\n");
	_printf_("    model_enum: " << model_enum << " " << EnumToStringx(model_enum) << "\n");
	_error_("not implemented yet");
}
/*}}}*/
int Cfsurfacesquaretransient::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cfsurfacesquaretransient::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int object_enum=CfsurfacesquaretransientEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->model_enum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->num_datatimes);
	marshallhandle->call(this->datatimes,this->num_datatimes);
	marshallhandle->call(this->passedflags,this->num_datatimes);
	marshallhandle->call(this->J);
} 
/*}}}*/
int Cfsurfacesquaretransient::ObjectEnum(void){/*{{{*/
	return CfsurfacesquaretransientEnum;
}
/*}}}*/

/*Definition virtual function resolutoin: */
int Cfsurfacesquaretransient::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cfsurfacesquaretransient::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cfsurfacesquaretransient::Response(FemModel* femmodel){/*{{{*/

	/*recover model time parameters: */
	IssmDouble time;
	femmodel->parameters->FindParam(&time,TimeEnum);

	/*Find closest datatime that is less than time*/
	int pos=-1;
	for(int i=0;i<this->num_datatimes;i++){
		if(this->datatimes[i]<=time){
			pos = i;
		}
		else{
			break;
		}
	}

	/*if pos=-1, time is earlier than the first data observation in this dataset*/
	if(pos==-1){
		_assert_(this->J==0.);
		return 0.;
	}

	/*Check that we have not yet calculated this cost function*/
	if(this->passedflags[pos]){
		return this->J;
	}

	/*Calculate cost function for this time slice*/
	IssmDouble J_part=0.;
	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		J_part+=this->Cfsurfacesquaretransient_Calculation(element,model_enum);
	}

	/*Sum across partition*/
	IssmDouble J_sum;
	ISSM_MPI_Allreduce((void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Record this cost function so that we do not recalculate it later*/
	this->passedflags[pos]= true;
	this->J += J_sum;

	/*Return full cost function this far*/
	return this->J;
}/*}}}*/
IssmDouble Cfsurfacesquaretransient::Cfsurfacesquaretransient_Calculation(Element* element, int model_enum){/*{{{*/

	IssmDouble Jelem=0.;
	IssmDouble misfit,Jdet;
	IssmDouble model,obs,weight;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnSurface()) return 0.;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

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
		datasetinput->GetInputValue(&weight,gauss,WeightsSurfaceObservationEnum);
		datasetinput->GetInputValue(&obs,gauss,SurfaceObservationEnum);
		model_input->GetInputValue(&model,gauss);

		/*Compute Misfit
		 *     
		 *       1  [           2 ]
		 *  J = --- | (x - x   )  |
		 *       2  [       obs   ]
		 **/
		misfit=0.5*(model-obs)*(model-obs);

		/*Add to cost function*/
		Jelem+=misfit*weight*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}/*}}}*/
