/*!\file BarystaticContributions.c
 * \brief: implementation of the BarystaticContributions object
 */

/*Include files: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./BarystaticContributions.h" 
#include "../toolkits/toolkits.h"
#include "./classes.h"
/*}}}*/

/*Constructors and destructors:*/
BarystaticContributions::BarystaticContributions(IoModel* iomodel ){ /*{{{*/

	int nel;

	iomodel->FetchData(&nice,"md.solidearth.npartice");
	if(nice){
		iomodel->FetchData(&pice,&nel,NULL,"md.solidearth.partitionice");
		ice=new Vector<IssmDouble>(nice);
		cumice=new Vector<IssmDouble>(nice); cumice->Set(0); cumice->Assemble();
	}
	else{
		ice=new Vector<IssmDouble>(1);
		cumice=new Vector<IssmDouble>(1);
	}

	iomodel->FetchData(&nhydro,"md.solidearth.nparthydro");
	if(nhydro){
		iomodel->FetchData(&phydro,&nel,NULL,"md.solidearth.partitionhydro");
		hydro=new Vector<IssmDouble>(nhydro);
		cumhydro=new Vector<IssmDouble>(nhydro); cumhydro->Set(0); cumhydro->Assemble();
	}
	else{
		hydro=new Vector<IssmDouble>(1);
		cumhydro=new Vector<IssmDouble>(1);
	}
	iomodel->FetchData(&nocean,"md.solidearth.npartocean");
	if(nocean){
		iomodel->FetchData(&pocean,&nel,NULL,"md.solidearth.partitionocean");
		ocean=new Vector<IssmDouble>(nocean);
		cumocean=new Vector<IssmDouble>(nocean); cumocean->Set(0); cumocean->Assemble();
	}
	else{
		ocean=new Vector<IssmDouble>(1);
		cumocean=new Vector<IssmDouble>(1);
	}

} /*}}}*/
BarystaticContributions::~BarystaticContributions(){ /*{{{*/
	delete ice;   delete cumice;
	delete hydro; delete cumhydro;
	delete ocean; delete cumocean;
	if(nice)xDelete<IssmDouble>(pice);
	if(nhydro)xDelete<IssmDouble>(phydro);
	if(nocean)xDelete<IssmDouble>(pocean);
}; /*}}}*/

/*Support routines:*/
IssmDouble BarystaticContributions::Total(){ /*{{{*/

	IssmDouble  sumice,sumhydro,sumocean;

	ice->Assemble();
	hydro->Assemble();
	ocean->Assemble();

	ice->Sum(&sumice);
	hydro->Sum(&sumhydro);
	ocean->Sum(&sumocean);

	return sumice+sumhydro+sumocean;

} /*}}}*/
IssmDouble BarystaticContributions::CumTotal(){ /*{{{*/

	IssmDouble sumice,sumhydro,sumocean;

	cumice->Assemble();
	cumhydro->Assemble();
	cumocean->Assemble();

	cumice->Sum(&sumice);
	cumhydro->Sum(&sumhydro);
	cumocean->Sum(&sumocean);

	return sumice+sumhydro+sumocean;

} /*}}}*/
void BarystaticContributions::Cumulate(Parameters* parameters){ /*{{{*/

	cumice->AXPY(ice,1);
	cumocean->AXPY(ocean,1);
	cumhydro->AXPY(hydro,1);

} /*}}}*/
void BarystaticContributions::Set(int eid, IssmDouble icevalue, IssmDouble hydrovalue, IssmDouble oceanvalue){ /*{{{*/

	int id;
	if(nice){
		id=reCast<int>(pice[eid]);
		ice->SetValue(id,icevalue,ADD_VAL);
	}
	else{
		ice->SetValue(0,icevalue,ADD_VAL);
	}

	if(nhydro){
		id=reCast<int>(phydro[eid]);
		hydro->SetValue(id,hydrovalue,ADD_VAL);
	}
	else{
		hydro->SetValue(0,hydrovalue,ADD_VAL);
	}

	if(nocean){
		id=reCast<int>(pocean[eid]);
		ocean->SetValue(id,oceanvalue,ADD_VAL);
	}
	else{
		ocean->SetValue(0,oceanvalue,ADD_VAL);
	}

} /*}}}*/
void BarystaticContributions::Reset(){ /*{{{*/

	ice->Set(0.);
	hydro->Set(0.);
	ocean->Set(0.);

} /*}}}*/
void BarystaticContributions::Save(Results* results, Parameters* parameters, IssmDouble oceanarea){ /*{{{*/

	int        step;
	IssmDouble time;
	IssmDouble rho_water;

	IssmDouble* cumice_serial=NULL;
	IssmDouble* cumhydro_serial=NULL;
	IssmDouble* cumocean_serial=NULL;

	IssmDouble sumice,sumhydro,sumocean;

	parameters->FindParam(&step,StepEnum);
	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&rho_water,MaterialsRhoSeawaterEnum);

	ice->Sum(&sumice); hydro->Sum(&sumhydro); ocean->Sum(&sumocean);
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,BslcEnum,this->Total()/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,BslcIceEnum,sumice/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,BslcHydroEnum,sumice/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,BslcOceanEnum,sumocean/oceanarea/rho_water,step,time));

	cumice->Sum(&sumice); cumhydro->Sum(&sumhydro); cumocean->Sum(&sumocean);
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,CumBslcEnum,this->CumTotal()/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,CumBslcIceEnum,sumice/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,CumBslcHydroEnum,sumhydro/oceanarea/rho_water,step,time));
	results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,CumBslcOceanEnum,sumocean/oceanarea/rho_water,step,time));

	if(nice){
		cumice_serial=this->cumice->ToMPISerial0(); for (int i=0;i<nice;i++)cumice_serial[i]=cumice_serial[i]/oceanarea/rho_water;
		results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,CumBslcIcePartitionEnum,cumice_serial,nice,1,step,time));
	}
	if(nhydro){
		cumhydro_serial=this->cumhydro->ToMPISerial0(); for (int i=0;i<nhydro;i++)cumhydro_serial[i]=cumhydro_serial[i]/oceanarea/rho_water;
		results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,CumBslcHydroPartitionEnum,cumhydro_serial,nhydro,1,step,time));
	}
	if(nocean){
		cumocean_serial=this->cumocean->ToMPISerial0(); for (int i=0;i<nocean;i++)cumocean_serial[i]=cumocean_serial[i]/oceanarea/rho_water;
		results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,CumBslcOceanPartitionEnum,cumocean_serial,nocean,1,step,time));
	}

	if(IssmComm::GetRank()==0){
		if(nice)xDelete<IssmDouble>(cumice_serial);
		if(nhydro)xDelete<IssmDouble>(cumhydro_serial);
		if(nocean)xDelete<IssmDouble>(cumocean_serial);
	}
	return;

} /*}}}*/
