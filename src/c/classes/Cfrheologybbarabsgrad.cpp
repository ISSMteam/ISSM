/*!\file Cfrheologybbarabsgrad.cpp
 * \brief: Cfrheologybbarabsgrad Object
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

/*Cfrheologybbarabsgrad constructors, destructors :*/
Cfrheologybbarabsgrad::Cfrheologybbarabsgrad(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->J = 0.;
	this->firsttimepassed = false;
}
/*}}}*/
Cfrheologybbarabsgrad::Cfrheologybbarabsgrad(char* in_name, int in_definitionenum){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J=0;
	this->firsttimepassed = false;
}
/*}}}*/
Cfrheologybbarabsgrad::Cfrheologybbarabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J=in_J;
	this->firsttimepassed = false;
}
/*}}}*/
Cfrheologybbarabsgrad::Cfrheologybbarabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J, bool in_firsttimepassed){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J=in_J;
	this->firsttimepassed = in_firsttimepassed;
}
/*}}}*/
Cfrheologybbarabsgrad::~Cfrheologybbarabsgrad(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Cfrheologybbarabsgrad::copy() {/*{{{*/
	Cfrheologybbarabsgrad* mf = new Cfrheologybbarabsgrad(this->name,this->definitionenum, this->J, this->firsttimepassed);
	return (Object*) mf;
}
/*}}}*/
void Cfrheologybbarabsgrad::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cfrheologybbarabsgrad::Echo(void){/*{{{*/
	_printf_(" Cfrheologybbarabsgrad: " << name << " " << this->definitionenum << "\n");
}
/*}}}*/
int Cfrheologybbarabsgrad::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cfrheologybbarabsgrad::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	/*ok, marshall operations: */
	int object_enum=CfrheologybbarabsgradEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->J);
	marshallhandle->call(this->firsttimepassed);
} 
/*}}}*/
int Cfrheologybbarabsgrad::ObjectEnum(void){/*{{{*/
	return CfrheologybbarabsgradEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Cfrheologybbarabsgrad::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cfrheologybbarabsgrad::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cfrheologybbarabsgrad::Response(FemModel* femmodel){/*{{{*/

	/*recover parameters: */
	IssmDouble J_part=0.;
	IssmDouble J_sum=0.;

	if (!this->firsttimepassed){
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			J_part+=this->Cfrheologybbarabsgrad_Calculation(element);
		}

		ISSM_MPI_Allreduce( (void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		this->J = J_sum;

		this->firsttimepassed = true;
	}
	return this->J;
}/*}}}*/
IssmDouble Cfrheologybbarabsgrad::Cfrheologybbarabsgrad_Calculation(Element* element){/*{{{*/

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
