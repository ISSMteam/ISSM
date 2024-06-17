/*!\file Cfdragcoeffabsgrad.cpp
 * \brief: Cfdragcoeffabsgrad Object
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

/*Cfdragcoeffabsgrad constructors, destructors :*/
Cfdragcoeffabsgrad::Cfdragcoeffabsgrad(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->J = 0.;
	this->firsttimepassed = false;
}
/*}}}*/
Cfdragcoeffabsgrad::Cfdragcoeffabsgrad(char* in_name, int in_definitionenum){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J = 0.;
	this->firsttimepassed = false;
}
/*}}}*/
Cfdragcoeffabsgrad::Cfdragcoeffabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J = in_J;
	this->firsttimepassed = false;
}
/*}}}*/
Cfdragcoeffabsgrad::Cfdragcoeffabsgrad(char* in_name, int in_definitionenum, IssmDouble in_J, bool in_firsttimepassed){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->J = in_J;
	this->firsttimepassed = in_firsttimepassed;
}
/*}}}*/
Cfdragcoeffabsgrad::~Cfdragcoeffabsgrad(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Cfdragcoeffabsgrad::copy() {/*{{{*/
	Cfdragcoeffabsgrad* mf = new Cfdragcoeffabsgrad(this->name,this->definitionenum, this->J, this->firsttimepassed);
	return (Object*) mf;
}
/*}}}*/
void Cfdragcoeffabsgrad::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cfdragcoeffabsgrad::Echo(void){/*{{{*/
	_printf_(" Cfdragcoeffabsgrad: " << name << " " << this->definitionenum << "\n");
}
/*}}}*/
int Cfdragcoeffabsgrad::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cfdragcoeffabsgrad::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	/*ok, marshall operations: */
   int object_enum=CfdragcoeffabsgradEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->J);
	marshallhandle->call(this->firsttimepassed);
} 
/*}}}*/
int Cfdragcoeffabsgrad::ObjectEnum(void){/*{{{*/
	return CfdragcoeffabsgradEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Cfdragcoeffabsgrad::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cfdragcoeffabsgrad::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cfdragcoeffabsgrad::Response(FemModel* femmodel){/*{{{*/

	/*recover parameters: */
	IssmDouble J_part=0.;
	IssmDouble J_sum=0.;

	if (!this->firsttimepassed){
		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			J_part+=this->Cfdragcoeffabsgrad_Calculation(element);
		}

		ISSM_MPI_Allreduce ( (void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		this->J=J_sum;

		this->firsttimepassed = true;
	}
	return this->J;
}/*}}}*/
IssmDouble Cfdragcoeffabsgrad::Cfdragcoeffabsgrad_Calculation(Element* element){/*{{{*/

	int        domaintype,numcomponents,frictionlaw;
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
	Input        *drag_input   = NULL;

	/* get the friction law: if 2-Weertman, 11-Schoof, 14-RegularizedCoulomb 15-RegularizedCoulomb2, which has a different names of C */
	element->FindParam(&frictionlaw, FrictionLawEnum);
	switch(frictionlaw) {
		case 2:
		case 11:
		case 14:
		case 15:
			drag_input = basalelement->GetInput(FrictionCEnum); _assert_(drag_input);
			break;
		default:
			drag_input = basalelement->GetInput(FrictionCoefficientEnum); _assert_(drag_input);
	}

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=basalelement->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		basalelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		datasetinput->GetInputValue(&weight,gauss,WeightsSurfaceObservationEnum);
		drag_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

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
