/*!\file OceanExchangeDatax
 * \brief: exchange of data with ocean model 
 */

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../modules.h"
#include "./OceanExchangeDatax.h"

void OceanExchangeDatax(FemModel* femmodel, bool init_stage){

	#ifndef _HAVE_AD_
	if(VerboseSolution()) _printf0_("   ocean coupling: exchanging information\n");
	int my_rank;
	ISSM_MPI_Comm tomitgcmcomm;
	ISSM_MPI_Status status;

	my_rank=IssmComm::GetRank();
	GenericParam<ISSM_MPI_Comm>* parcom = dynamic_cast<GenericParam<ISSM_MPI_Comm>*>(femmodel->parameters->FindParamObject(ToMITgcmCommEnum));
	if(!parcom)_error_("TransferForcing error message: could not find ToMITgcmCommEnum communicator");
	tomitgcmcomm=parcom->GetParameterValue();

	int oceangridnxsize,oceangridnysize,ngrids_ocean,nels_ocean,isoceancoupling;
	IssmDouble  oceantime,coupling_time,time,yts;
	IssmDouble rho_ice;
	IssmDouble *oceanmelt         = NULL;
	IssmDouble *oceangridx;
	IssmDouble *oceangridy;
	IssmDouble *icethickness_oceangrid = NULL;
	IssmDouble *icemask_oceangrid = NULL;
	IssmDouble* x_ice             = NULL;
	IssmDouble* y_ice             = NULL;
	IssmDouble* lat_ice           = NULL;
	IssmDouble* lon_ice           = NULL;
	IssmDouble* icethickness      = NULL;
	IssmDouble* icemask           = NULL;
	IssmDouble* melt_mesh         = NULL;
	int*        index_ice         = NULL;
	int*        index_ocean       = NULL;
	int         ngrids_ice=femmodel->vertices->NumberOfVertices();
	int         nels_ice=femmodel->elements->NumberOfElements();

	/*Recover fixed parameters and store them*/
	femmodel->parameters->FindParam(&coupling_time,TimesteppingCouplingTimeEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&isoceancoupling,TransientIsoceancouplingEnum);

	/*Exchange or recover mesh and inputs needed*/
	if(init_stage==true){
		if(my_rank==0){
			ISSM_MPI_Send(&coupling_time,1,ISSM_MPI_DOUBLE,0,10001000,tomitgcmcomm);
			ISSM_MPI_Recv(&oceangridnxsize,1,ISSM_MPI_INT,0,10001003,tomitgcmcomm,&status);
			ISSM_MPI_Recv(&oceangridnysize,1,ISSM_MPI_INT,0,10001004,tomitgcmcomm,&status);
		}
		ngrids_ocean=oceangridnxsize*oceangridnysize;
		ISSM_MPI_Bcast(&oceangridnxsize,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(&oceangridnysize,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(&ngrids_ocean,1,ISSM_MPI_INT,0,IssmComm::GetComm());
		femmodel->parameters->SetParam(oceangridnxsize,OceanGridNxEnum);
		femmodel->parameters->SetParam(oceangridnysize,OceanGridNyEnum);
		oceangridx=xNew<IssmDouble>(ngrids_ocean);
		oceangridy=xNew<IssmDouble>(ngrids_ocean);
		if(my_rank==0){
			ISSM_MPI_Recv(oceangridx,ngrids_ocean,ISSM_MPI_DOUBLE,0,10001005,tomitgcmcomm,&status);
			ISSM_MPI_Recv(oceangridy,ngrids_ocean,ISSM_MPI_DOUBLE,0,10001006,tomitgcmcomm,&status);
		}

		ISSM_MPI_Bcast(oceangridx,ngrids_ocean,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(oceangridy,ngrids_ocean,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		femmodel->parameters->SetParam(oceangridx,ngrids_ocean,OceanGridXEnum);
		femmodel->parameters->SetParam(oceangridy,ngrids_ocean,OceanGridYEnum);
	}
	else{
		/*Recoved ocean grid from parameters*/
		femmodel->parameters->FindParam(&oceangridx,&ngrids_ocean,OceanGridXEnum);
		femmodel->parameters->FindParam(&oceangridy,&ngrids_ocean,OceanGridYEnum);
	}

	/*Interpolate ice thickness and mask onto ocean grid*/
	femmodel->GetMesh(femmodel->vertices,femmodel->elements,&x_ice,&y_ice,&index_ice);
	BamgTriangulatex(&index_ocean,&nels_ocean,oceangridx,oceangridy,ngrids_ocean);
	if(isoceancoupling==2){
		femmodel->vertices->LatLonList(&lat_ice,&lon_ice);
	}
	else{
		femmodel->vertices->XYList(&lon_ice,&lat_ice);
	}
	GetVectorFromInputsx(&icethickness,femmodel,ThicknessEnum,VertexSIdEnum);
	Options* options = new Options();
	GenericOption<double> *odouble = new GenericOption<double>();
	const char* name = "default";
	odouble->name =xNew<char>(strlen(name)+1);
	memcpy(odouble->name,name,(strlen(name)+1)*sizeof(char));
	odouble->value=+9999.;
	odouble->size[0]=1;
	odouble->size[1]=1;
	options->AddOption(odouble);
	InterpFromMeshToMesh2dx(&icethickness_oceangrid,index_ice,lon_ice,lat_ice,ngrids_ice,nels_ice,
					icethickness,ngrids_ice,1,oceangridx,oceangridy,ngrids_ocean,options);
	delete options;
	xDelete<IssmDouble>(icethickness);

	GetVectorFromInputsx(&icemask,femmodel,MaskIceLevelsetEnum,VertexSIdEnum);
	Options* options2 = new Options();
	GenericOption<double> *odouble2 = new GenericOption<double>();
	const char* name2 = "default";
	odouble2->name =xNew<char>(strlen(name2)+1);
	memcpy(odouble2->name,name2,(strlen(name2)+1)*sizeof(char));
	odouble2->value=+1.;
	odouble2->size[0]=1;
	odouble2->size[1]=1;
	options2->AddOption(odouble2);
	InterpFromMeshToMesh2dx(&icemask_oceangrid,index_ice,lon_ice,lat_ice,ngrids_ice,nels_ice,
				icemask,ngrids_ice,1,oceangridx,oceangridy,ngrids_ocean,options2);
	delete options2;
	xDelete<IssmDouble>(icemask);

	/*Put +9999 for places where there is no ice!*/
	femmodel->parameters->FindParam(&rho_ice,MaterialsRhoIceEnum);
	for(int i=0;i<ngrids_ocean;i++) icethickness_oceangrid[i]=icethickness_oceangrid[i]*rho_ice; //ocean needs ice mass in kg/m^2
	for(int i=0;i<ngrids_ocean;i++) if(icemask_oceangrid[i]>0.) icethickness_oceangrid[i]=+9999.;
	xDelete<IssmDouble>(icemask_oceangrid);

	if(init_stage==true){ //just send icethickness
		if(my_rank==0){
			ISSM_MPI_Send(icethickness_oceangrid,ngrids_ocean,ISSM_MPI_DOUBLE,0,10001008,tomitgcmcomm);
		}
	}
	else{ //send and receive exchanged data
		femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
		if(my_rank==0){
			ISSM_MPI_Send(icethickness_oceangrid,ngrids_ocean,ISSM_MPI_DOUBLE,0,10001008,tomitgcmcomm);
			ISSM_MPI_Send(&time,1,ISSM_MPI_DOUBLE,0,10001001,tomitgcmcomm);
			ISSM_MPI_Recv(&oceantime,1,ISSM_MPI_DOUBLE,0,10001002,tomitgcmcomm,&status);
			if((oceantime - time > 0.1*yts) & (oceantime - time < -0.1*yts)) _error_("Ocean and ice time are starting to diverge");
			oceanmelt = xNew<IssmDouble>(ngrids_ocean);
			ISSM_MPI_Recv(oceanmelt,ngrids_ocean,ISSM_MPI_DOUBLE,0,10001007,tomitgcmcomm,&status);
		}
		ISSM_MPI_Bcast(&oceantime,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
		if(my_rank!=0) oceanmelt=xNew<IssmDouble>(ngrids_ocean);
		ISSM_MPI_Bcast(oceanmelt,ngrids_ocean,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		/*Interp melt onto ice grid*/
		InterpFromMeshToMesh2dx(&melt_mesh,index_ocean,oceangridx,oceangridy,ngrids_ocean,nels_ocean,
					oceanmelt,ngrids_ocean,1,
					lon_ice,lat_ice,ngrids_ice,NULL);

		for(int i=0;i<ngrids_ice;i++) melt_mesh[i]=-melt_mesh[i]/rho_ice; //heat flux provided by ocean is in kg/m^2/s
		InputUpdateFromVectorx(femmodel,melt_mesh,BasalforcingsFloatingiceMeltingRateEnum,VertexSIdEnum);
	}

	/*Delete*/
	xDelete<int>(index_ice);
	xDelete<int>(index_ocean);
	xDelete<IssmDouble>(lat_ice);
	xDelete<IssmDouble>(lon_ice);
	xDelete<IssmDouble>(x_ice);
	xDelete<IssmDouble>(y_ice);
	xDelete<IssmDouble>(icethickness_oceangrid);
	xDelete<IssmDouble>(oceangridx);
	xDelete<IssmDouble>(oceangridy);
	xDelete<IssmDouble>(melt_mesh);
	xDelete<IssmDouble>(oceanmelt);
	#else
	_error_("not supported");
	#endif
}
