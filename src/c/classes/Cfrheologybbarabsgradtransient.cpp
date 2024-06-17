/*!\file Cfrheologybbarabsgradtransient.cpp
 * \brief: Cfrheologybbarabsgradtransient Object
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

/*Cfrheologybbarabsgradtransient constructors, destructors :*/
Cfrheologybbarabsgradtransient::Cfrheologybbarabsgradtransient(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->datatimes         = NULL;
	this->passedflags   = NULL;
	this->J = 0.;
}
/*}}}*/
Cfrheologybbarabsgradtransient::Cfrheologybbarabsgradtransient(char* in_name, int in_definitionenum, int in_num_datatimes, IssmDouble* in_datatimes, bool* in_passedflags, IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->num_datatimes = in_num_datatimes;

	/*Allocate arrays*/
	_assert_(this->num_datatimes>0);
	this->datatimes   = xNew<IssmDouble>(this->num_datatimes);
	this->passedflags = xNew<bool>(this->num_datatimes);
	xMemCpy<IssmDouble>(this->datatimes,in_datatimes,this->num_datatimes);
	xMemCpy<bool>(this->passedflags,in_passedflags,this->num_datatimes);

	this->J = in_J;
}
/*}}}*/
Cfrheologybbarabsgradtransient::Cfrheologybbarabsgradtransient(char* in_name, int in_definitionenum, int in_num_datatimes, IssmDouble* in_datatimes){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

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
Cfrheologybbarabsgradtransient::~Cfrheologybbarabsgradtransient(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Cfrheologybbarabsgradtransient::copy() {/*{{{*/
	Cfrheologybbarabsgradtransient* mf = new Cfrheologybbarabsgradtransient(this->name,this->definitionenum, this->num_datatimes, this->datatimes, this->passedflags, this->J);
	return (Object*) mf;
}
/*}}}*/
void Cfrheologybbarabsgradtransient::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cfrheologybbarabsgradtransient::Echo(void){/*{{{*/
	_printf_(" Cfrheologybbarabsgradtransient: " << name << " " << this->definitionenum << "\n");
	_error_("not implemented yet");
}
/*}}}*/
int Cfrheologybbarabsgradtransient::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cfrheologybbarabsgradtransient::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	/*ok, marshall operations: */
	int object_enum=CfrheologybbarabsgradtransientEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->num_datatimes);
   marshallhandle->call(this->datatimes,this->num_datatimes);
   marshallhandle->call(this->passedflags,this->num_datatimes);
   marshallhandle->call(this->J);
} 
/*}}}*/
int Cfrheologybbarabsgradtransient::ObjectEnum(void){/*{{{*/
	return CfrheologybbarabsgradtransientEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Cfrheologybbarabsgradtransient::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cfrheologybbarabsgradtransient::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cfrheologybbarabsgradtransient::Response(FemModel* femmodel){/*{{{*/

	/*recover time parameters: */
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
	IssmDouble J_sum=0.;

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		J_part+=this->Cfrheologybbarabsgradtransient_Calculation(element);
	}

	ISSM_MPI_Allreduce( (void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	this->passedflags[pos]= true;
	this->J += J_sum;

	return this->J;
}/*}}}*/
IssmDouble Cfrheologybbarabsgradtransient::Cfrheologybbarabsgradtransient_Calculation(Element* element){/*{{{*/

	int        domaintype,numcomponents;
	IssmDouble Jelem=0.;
	IssmDouble Jdet;
	IssmDouble dp[2],weight;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnBase()) return 0.;

	/*If on water, return 0: */
	if(!element->IsIceInElement()) return 0.;

	/*Get problem dimension*/
	element->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DverticalEnum:   numcomponents   = 1; break;
		case Domain3DEnum:           numcomponents   = 2; break;
		case Domain2DhorizontalEnum: numcomponents   = 2; break;
		default: _error_("not supported yet");
	}

	/*Spawn surface element*/
	Element* basalelement = element->SpawnBasalElement();

	/* Get node coordinates*/
	basalelement->GetVerticesCoordinates(&xyz_list);

	/*Get input if it already exists*/
	DatasetInput *datasetinput = basalelement->GetDatasetInput(definitionenum);  _assert_(datasetinput);
	Input* rheologyb_input=basalelement->GetInput(MaterialsRheologyBbarEnum);                  _assert_(rheologyb_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		datasetinput->GetInputValue(&weight,gauss,WeightsSurfaceObservationEnum);
		rheologyb_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

		/*Add to cost function*/
		Jelem+=weight*.5*dp[0]*dp[0]*Jdet*gauss->weight;
		if(numcomponents==2) Jelem+=weight*.5*dp[1]*dp[1]*Jdet*gauss->weight;
	}

	/*clean up and Return: */
	if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}/*}}}*/
