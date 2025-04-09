/*!\file Cfsurfacelogvel.cpp
 * \brief: Cfsurfacelogvel Object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <float.h> /*defines DBL_EPSILON*/
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

/*Cfsurfacelogvel constructors, destructors :*/
Cfsurfacelogvel::Cfsurfacelogvel(){/*{{{*/

	this->definitionenum = -1;
	this->name = NULL;
	this->datatime=0.;
	this->timepassedflag = false;
	this->J = 0.;

}
/*}}}*/
Cfsurfacelogvel::Cfsurfacelogvel(char* in_name, int in_definitionenum, IssmDouble in_datatime){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->datatime=in_datatime;

	this->timepassedflag=false;
	this->J=0.;

}
/*}}}*/
Cfsurfacelogvel::Cfsurfacelogvel(char* in_name, int in_definitionenum, IssmDouble in_datatime, bool in_timepassedflag,IssmDouble in_J){/*{{{*/

	this->definitionenum=in_definitionenum;

	this->name		= xNew<char>(strlen(in_name)+1);
	xMemCpy<char>(this->name,in_name,strlen(in_name)+1);

	this->datatime=in_datatime;
	this->timepassedflag=in_timepassedflag;
	this->J=in_J;

}
/*}}}*/
Cfsurfacelogvel::~Cfsurfacelogvel(){/*{{{*/
	if(this->name)xDelete(this->name);
}
/*}}}*/
/*Object virtual function resolutoin: */
Object* Cfsurfacelogvel::copy() {/*{{{*/
	Cfsurfacelogvel* mf = new Cfsurfacelogvel(this->name,this->definitionenum,this->datatime,this->timepassedflag, this->J);
	return (Object*) mf;
}
/*}}}*/
void Cfsurfacelogvel::DeepEcho(void){/*{{{*/
	this->Echo();
}
/*}}}*/
void Cfsurfacelogvel::Echo(void){/*{{{*/
	_printf_(" Cfsurfacelogvel: " << name << " " << this->definitionenum << "\n");
	_printf_("    datatime: " << datatime << "\n");
	_printf_("	  timepassedflag: "<<timepassedflag<<"\n");
	_printf_("	  J: "<<J<<"\n");
}
/*}}}*/
int Cfsurfacelogvel::Id(void){/*{{{*/
	return -1;
}
/*}}}*/
void Cfsurfacelogvel::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	int object_enum=CfsurfacelogvelEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->definitionenum);
	marshallhandle->call(this->name);
	marshallhandle->call(this->datatime);
	marshallhandle->call(this->timepassedflag);
	marshallhandle->call(this->J);
} 
/*}}}*/
int Cfsurfacelogvel::ObjectEnum(void){/*{{{*/
	return CfsurfacelogvelEnum;
}
/*}}}*/
/*Definition virtual function resolutoin: */
int Cfsurfacelogvel::DefinitionEnum(){/*{{{*/
	return this->definitionenum;
}
/*}}}*/
char* Cfsurfacelogvel::Name(){/*{{{*/
	char* name2=xNew<char>(strlen(this->name)+1);
	xMemCpy(name2,this->name,strlen(this->name)+1);

	return name2;
}
/*}}}*/
IssmDouble Cfsurfacelogvel::Response(FemModel* femmodel){/*{{{*/

	/*recover time parameters: */
	IssmDouble time;
	femmodel->parameters->FindParam(&time,TimeEnum);

	if(this->datatime<=time && !this->timepassedflag){

		IssmDouble J_part=0.;
		IssmDouble J_sum=0.;

		for(Object* & object : femmodel->elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			J_part+=this->Cfsurfacelogvel_Calculation(element,definitionenum);
		}

		ISSM_MPI_Allreduce ( (void*)&J_part,(void*)&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
		ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		this->timepassedflag = true;
		this->J = J_sum;
	}

	return this->J;
}/*}}}*/
IssmDouble Cfsurfacelogvel::Cfsurfacelogvel_Calculation(Element* element, int definitionenum){/*{{{*/

	int        domaintype,numcomponents;
	IssmDouble Jelem=0.;
	IssmDouble epsvel=DBL_EPSILON;
	IssmDouble meanvel=3.170979198376458e-05; /*1000 m/yr*/
	IssmDouble velocity_mag,obs_velocity_mag;
	IssmDouble misfit,Jdet;
	IssmDouble vx,vy,vxobs,vyobs,weight;
	IssmDouble* xyz_list = NULL;

	/*Get basal element*/
	if(!element->IsOnSurface()) return 0.;

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
	Element* topelement = element->SpawnTopElement();

	/* Get node coordinates*/
	topelement->GetVerticesCoordinates(&xyz_list);

	/*Get model values*/
	Input *vx_input = topelement->GetInput(VxEnum); _assert_(vx_input);
	Input *vy_input = NULL;
	if(numcomponents==2){
		vy_input = topelement->GetInput(VyEnum); _assert_(vy_input);
	}

	/*Retrieve all inputs we will be needing: */
	DatasetInput *datasetinput = topelement->GetDatasetInput(definitionenum); _assert_(datasetinput);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=topelement->NewGauss(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		topelement->JacobianDeterminant(&Jdet,xyz_list,gauss);

		/*Get all parameters at gaussian point*/
		datasetinput->GetInputValue(&weight,gauss,WeightsSurfaceObservationEnum);
		vx_input->GetInputValue(&vx,gauss);
		datasetinput->GetInputValue(&vxobs,gauss,VxObsEnum);
		if(numcomponents==2){
			vy_input->GetInputValue(&vy,gauss);
			datasetinput->GetInputValue(&vyobs,gauss,VyObsEnum);
		}

		/*Compute SurfaceLogVelMisfit:
		 *        *                 [        vel + eps     ] 2
		 *               * J = 4 \bar{v}^2 | log ( -----------  ) |
		 *                          [       vel   + eps    ]
		 *                                      obs
		 */
		if(numcomponents==1){
			velocity_mag    =fabs(vx)+epsvel;
			obs_velocity_mag=fabs(vxobs)+epsvel;
		}
		else{
			velocity_mag    =sqrt(vx*vx+vy*vy)+epsvel;
			obs_velocity_mag=sqrt(vxobs*vxobs+vyobs*vyobs)+epsvel;
		}

		misfit=4*pow(meanvel,2)*pow(log(velocity_mag/obs_velocity_mag),2);

		/*Add to cost function*/
		Jelem+=misfit*weight*Jdet*gauss->weight;

	}

	/*clean up and Return: */
	if(topelement->IsSpawnedElement()){topelement->DeleteMaterials(); delete topelement;};
	xDelete<IssmDouble>(xyz_list);
	delete gauss;
	return Jelem;
}/*}}}*/
