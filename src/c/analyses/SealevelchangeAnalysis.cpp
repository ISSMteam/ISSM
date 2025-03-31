#include "./SealevelchangeAnalysis.h"
#include <math.h>
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Inputs/TransientInput.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../cores/cores.h"

/*Model processing*/
void SealevelchangeAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void SealevelchangeAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void SealevelchangeAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,SealevelchangeAnalysisEnum,P1Enum);
}/*}}}*/
int  SealevelchangeAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SealevelchangeAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int isexternal=0;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Create inputs: */
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);

	iomodel->FetchDataToInput(inputs,elements,"md.solidearth.transfercount",CouplingTransferCountEnum);

	/*external solidearthsolution: solid-Earth model*/
	iomodel->FetchData(&isexternal,"md.solidearth.isexternal");

	if(isexternal){
		iomodel->FetchDataToInput(inputs,elements,"md.solidearth.external.displacementeast",SolidearthExternalDisplacementEastRateEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.solidearth.external.displacementnorth",SolidearthExternalDisplacementNorthRateEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.solidearth.external.displacementup",SolidearthExternalDisplacementUpRateEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.solidearth.external.geoid",SolidearthExternalGeoidRateEnum);

		/*Resolve Mmes using the modelid, if necessary:*/
		if (inputs->GetInputObjectEnum(SolidearthExternalDisplacementUpRateEnum)==DatasetInputEnum){
			int modelid;

			/*retrieve model id: */
			iomodel->FetchData(&modelid,"md.solidearth.external.modelid");

			/*replace dataset of forcings with only one, the modelid'th:*/
			MmeToInputFromIdx(inputs,elements,NULL,modelid,SolidearthExternalDisplacementNorthRateEnum, P1Enum);
			MmeToInputFromIdx(inputs,elements,NULL,modelid,SolidearthExternalDisplacementEastRateEnum, P1Enum);
			MmeToInputFromIdx(inputs,elements,NULL,modelid,SolidearthExternalDisplacementUpRateEnum, P1Enum);
			MmeToInputFromIdx(inputs,elements,NULL,modelid,SolidearthExternalGeoidRateEnum, P1Enum);
		}
	}

	/*Initialize solid earth motion and sea level: */
	iomodel->ConstantToInput(inputs,elements,0.,BedEastEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0.,BedNorthEnum,P1Enum);
    iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum);

	/*Initialize loads: no! this should be done by the corresponding mass transports!*/
	iomodel->ConstantToInput(inputs,elements,0.,DeltaTwsEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0.,DeltaIceThicknessEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0.,DeltaBottomPressureEnum,P1Enum);

}/*}}}*/
void SealevelchangeAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         nl;
	int         ntimesteps;
	IssmDouble* love_h=NULL;
	IssmDouble* love_k=NULL;
	IssmDouble* love_l=NULL;
	IssmDouble* love_th=NULL;
	IssmDouble* love_tk=NULL;
	IssmDouble* love_tl=NULL;
	IssmDouble* love_pmtf_colinear=NULL;
	IssmDouble* love_pmtf_ortho=NULL;
	IssmDouble* love_timefreq=NULL;
	bool        love_istime=true;
	int         externalnature=0;
	int         isexternal=0;

	IssmDouble* G_gravi = NULL;
	IssmDouble* G_gravi_local = NULL;
	IssmDouble* G_viscoelastic = NULL;
	IssmDouble* G_viscoelastic_interpolated = NULL;
	IssmDouble* G_viscoelastic_local = NULL;
	IssmDouble* U_viscoelastic = NULL;
	IssmDouble* U_viscoelastic_interpolated = NULL;
	IssmDouble* U_viscoelastic_local = NULL;
	IssmDouble* H_viscoelastic = NULL;
	IssmDouble* H_viscoelastic_interpolated= NULL;
	IssmDouble* H_viscoelastic_local = NULL;
	IssmDouble* Pmtf_col_interpolated = NULL;
	IssmDouble* Pmtf_ortho_interpolated = NULL;
	IssmDouble* Pmtf_z_interpolated = NULL;
	IssmDouble* Love_th2_interpolated = NULL;
	IssmDouble* Love_tk2_interpolated = NULL;
	IssmDouble* Love_tl2_interpolated = NULL;

	int         M,m,lower_row,upper_row;
	IssmDouble  degacc=.01;
	IssmDouble  timeacc;
	IssmDouble  planetradius=0;
	IssmDouble  planetarea=0;
	IssmDouble  constant=0;
	IssmDouble  rho_earth;
	int		isgrd=0;
	bool		selfattraction=false;
	bool		elastic=false;
	bool		viscous=false;
	bool		rotation=false;
	int         ndeg;
	int         horiz;

	bool istime=true;
	IssmDouble start_time,final_time;
	int  nt,precomputednt;

	int     numoutputs;
	char**  requestedoutputs = NULL;
	int* recvcounts = NULL;
	int* displs=NULL;

	/*transition vectors: */
	IssmDouble **transitions    = NULL;
	int         *transitions_M    = NULL;
	int         *transitions_N    = NULL;
	int          ntransitions;
	IssmDouble*  partitionice=NULL;
	IssmDouble*  partitionhydro=NULL;
	IssmDouble*  partitionocean=NULL;
	IssmDouble*  bslcice_partition=NULL;
	IssmDouble*  bslchydro_partition=NULL;
	IssmDouble*  bslcocean_partition=NULL;
	int          npartice,nparthydro,npartocean,nel;
	int          grdmodel;

	/*some constant parameters: */
	parameters->AddObject(iomodel->CopyConstantObject("md.dsl.model",DslModelEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.isexternal",SolidearthIsExternalEnum));
	iomodel->FetchData(&isexternal,"md.solidearth.isexternal");
	if(isexternal) parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.external.nature",SolidearthExternalNatureEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.runfrequency",SolidearthSettingsRunFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.degacc",SolidearthSettingsDegreeAccuracyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.reltol",SolidearthSettingsReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.abstol",SolidearthSettingsAbstolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.maxiter",SolidearthSettingsMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.selfattraction",SolidearthSettingsSelfAttractionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.horiz",SolidearthSettingsHorizEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.elastic",SolidearthSettingsElasticEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.viscous",SolidearthSettingsViscousEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.rotation",SolidearthSettingsRotationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.equatorialmoi",RotationalEquatorialMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.polarmoi",RotationalPolarMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.angularvelocity",RotationalAngularVelocityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.grdocean",SolidearthSettingsGrdOceanEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.ocean_area_scaling",SolidearthSettingsOceanAreaScalingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.sealevelloading",SolidearthSettingsSealevelLoadingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.isgrd",SolidearthSettingsGRDEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.compute_bp_grd",SolidearthSettingsComputeBpGrdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.planetradius",SolidearthPlanetRadiusEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.cross_section_shape",SolidearthSettingsCrossSectionShapeEnum));

	parameters->AddObject(new DoubleParam(CumBslcEnum,0.0));
	parameters->AddObject(new DoubleParam(CumBslcIceEnum,0.0));
	parameters->AddObject(new DoubleParam(CumBslcHydroEnum,0.0));
	parameters->AddObject(new DoubleParam(CumGmtslcEnum,0.0));
	/*compute planet area and plug into parameters:*/
	iomodel->FetchData(&planetradius,"md.solidearth.planetradius");
	iomodel->FetchData(&rho_earth,"md.materials.earth_density");
	planetarea=4*M_PI*planetradius*planetradius;
	parameters->AddObject(new DoubleParam(SolidearthPlanetAreaEnum,planetarea));

	/*Deal with partition of the barystatic contribution:*/
	iomodel->FetchData(&npartice,"md.solidearth.npartice");
	parameters->AddObject(new IntParam(SolidearthNpartIceEnum,npartice));
	if(npartice){
		iomodel->FetchData(&partitionice,&nel,NULL,"md.solidearth.partitionice");
		parameters->AddObject(new DoubleMatParam(SolidearthPartitionIceEnum,partitionice,nel,1));
		parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.npartice",SolidearthNpartIceEnum));
		bslcice_partition=xNewZeroInit<IssmDouble>(npartice);
		parameters->AddObject(new DoubleMatParam(CumBslcIcePartitionEnum,bslcice_partition,npartice,1));
		xDelete<IssmDouble>(partitionice);
	}

	iomodel->FetchData(&nparthydro,"md.solidearth.nparthydro");
	parameters->AddObject(new IntParam(SolidearthNpartHydroEnum,nparthydro));
	if(nparthydro){
		iomodel->FetchData(&partitionhydro,&nel,NULL,"md.solidearth.partitionhydro");
		parameters->AddObject(new DoubleMatParam(SolidearthPartitionHydroEnum,partitionhydro,nel,1));
		parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.nparthydro",SolidearthNpartHydroEnum));
		bslchydro_partition=xNewZeroInit<IssmDouble>(nparthydro);
		parameters->AddObject(new DoubleMatParam(CumBslcHydroPartitionEnum,bslchydro_partition,nparthydro,1));
		xDelete<IssmDouble>(partitionhydro);
	}

	iomodel->FetchData(&npartocean,"md.solidearth.npartocean");
	parameters->AddObject(new IntParam(SolidearthNpartOceanEnum,npartocean));
	if(npartocean){
		iomodel->FetchData(&partitionocean,&nel,NULL,"md.solidearth.partitionocean");
		parameters->AddObject(new DoubleMatParam(SolidearthPartitionOceanEnum,partitionocean,nel,1));
		parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.npartocean",SolidearthNpartOceanEnum));
		bslcocean_partition=xNewZeroInit<IssmDouble>(npartocean);
		parameters->AddObject(new DoubleMatParam(CumBslcOceanPartitionEnum,bslcocean_partition,npartocean,1));
		xDelete<IssmDouble>(partitionocean);
	}
	/*New optimized code:*/
	ToolkitsOptionsFromAnalysis(parameters,SealevelchangeAnalysisEnum); //this is requested by the BarystaticContributions class inner vectors.
	BarystaticContributions* barystaticcontributions=new BarystaticContributions(iomodel);
	parameters->AddObject(new GenericParam<BarystaticContributions*>(barystaticcontributions,BarystaticContributionsEnum));

	/*Deal with external multi-model ensembles: {{{*/
	if(isexternal){
		iomodel->FetchData(&externalnature,"md.solidearth.external.nature");
		if(externalnature>=3){
			IssmDouble modelid; 
			int nummodels;

			/*create double param, not int param, because Dakota will be updating it as a double potentially: */
			iomodel->FetchData(&modelid,"md.solidearth.external.modelid");
			parameters->AddObject(new DoubleParam(SolidearthExternalModelidEnum,modelid));
			parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.external.nummodels",SolidearthExternalNummodelsEnum));
			iomodel->FetchData(&nummodels,"md.solidearth.external.nummodels");

			/*quick checks: */
			if(nummodels<=0)_error_("mme solidearth solution object in  md.solidearth.external field should contain at least 1 ensemble model!");
			if(modelid<=0 || modelid>nummodels)_error_("modelid field in md.solidearth.external field should be between 1 and the number of ensemble runs!");
		} /*}}}*/
	}

	/*Indicate we have not yet run the Geometry Core module: */
	parameters->AddObject(new BoolParam(SealevelchangeGeometryDoneEnum,false));

	parameters->FindParam(&grdmodel,GrdModelEnum);
	parameters->FindParam(&isgrd,SolidearthSettingsGRDEnum);
	if(grdmodel==ElasticEnum && isgrd){
		/*Deal with elasticity {{{*/
		iomodel->FetchData(&selfattraction,"md.solidearth.settings.selfattraction");
		iomodel->FetchData(&elastic,"md.solidearth.settings.elastic");
		iomodel->FetchData(&viscous,"md.solidearth.settings.viscous");
		iomodel->FetchData(&rotation,"md.solidearth.settings.rotation");
		iomodel->FetchData(&horiz,"md.solidearth.settings.horiz");

		if(selfattraction){
			/*compute green functions for a range of angles*/
			iomodel->FetchData(&degacc,"md.solidearth.settings.degacc");
			M=reCast<int,IssmDouble>(180./degacc+1.);
		}

		//default values
		nt=1;
		ntimesteps=1;
		/*love numbers: */
		if(viscous || elastic){
			int dummy;

			iomodel->FetchData(&timeacc,"md.solidearth.settings.timeacc");
			iomodel->FetchData(&start_time,"md.timestepping.start_time");
			iomodel->FetchData(&final_time,"md.timestepping.final_time");
			iomodel->FetchData(&love_istime,"md.solidearth.lovenumbers.istime");
			iomodel->FetchData(&love_timefreq,&precomputednt,&dummy,"md.solidearth.lovenumbers.timefreq");
			iomodel->FetchData(&love_h,&ndeg,&precomputednt,"md.solidearth.lovenumbers.h");
			iomodel->FetchData(&love_k,&ndeg,&precomputednt,"md.solidearth.lovenumbers.k");
			iomodel->FetchData(&love_l,&ndeg,&precomputednt,"md.solidearth.lovenumbers.l");

			parameters->AddObject(new DoubleParam(SolidearthSettingsTimeAccEnum,timeacc));
			//parameters->AddObject(new DoubleMatParam(LoadLoveHEnum,love_h,ndeg,precomputednt));
			//parameters->AddObject(new DoubleMatParam(LoadLoveKEnum,love_k,ndeg,precomputednt));
			//parameters->AddObject(new DoubleMatParam(LoadLoveLEnum,love_l,ndeg,precomputednt));

			if (rotation){
				iomodel->FetchData(&love_th,&ndeg,&precomputednt,"md.solidearth.lovenumbers.th");
				iomodel->FetchData(&love_tk,&ndeg,&precomputednt,"md.solidearth.lovenumbers.tk");
				iomodel->FetchData(&love_tl,&ndeg,&precomputednt,"md.solidearth.lovenumbers.tl");
				iomodel->FetchData(&love_pmtf_colinear,&dummy,&precomputednt,"md.solidearth.lovenumbers.pmtf_colinear");
				iomodel->FetchData(&love_pmtf_ortho,&dummy,&precomputednt,"md.solidearth.lovenumbers.pmtf_ortho");

				parameters->AddObject(new DoubleMatParam(LovePolarMotionTransferFunctionColinearEnum,love_pmtf_colinear,1,precomputednt));
				parameters->AddObject(new DoubleMatParam(LovePolarMotionTransferFunctionOrthogonalEnum,love_pmtf_ortho,1,precomputednt));
				//parameters->AddObject(new DoubleMatParam(TidalLoveHEnum,love_th,ndeg,precomputednt));
				//parameters->AddObject(new DoubleMatParam(TidalLoveKEnum,love_tk,ndeg,precomputednt));
				//parameters->AddObject(new DoubleMatParam(TidalLoveLEnum,love_tl,ndeg,precomputednt));
			}

			parameters->AddObject(new DoubleMatParam(LoveTimeFreqEnum,love_timefreq,precomputednt,1));
			parameters->AddObject(new BoolParam(LoveIsTimeEnum,love_istime));

			// AD performance is sensitive to calls to ensurecontiguous.
			// // Providing "t" will cause ensurecontiguous to be called.
			if(viscous){
				IssmDouble* viscoustimes=NULL;
				ntimesteps=precomputednt; 
				nt=reCast<int,IssmDouble>((final_time-start_time)/timeacc)+1;

				parameters->AddObject(new IntParam(SealevelchangeViscousNumStepsEnum,nt));
				/*Initialize viscous stack times:*/
				viscoustimes=xNew<IssmDouble>(nt);
				for(int t=0;t<nt;t++){
					viscoustimes[t]=start_time+timeacc*t;
				}
				parameters->AddObject(new DoubleVecParam(SealevelchangeViscousTimesEnum,viscoustimes,nt));
				parameters->AddObject(new IntParam(SealevelchangeViscousIndexEnum,0));
				xDelete<IssmDouble>(viscoustimes);
				if (rotation){
					IssmDouble* viscouspolarmotion=NULL;
					viscouspolarmotion=xNewZeroInit<IssmDouble>(3*nt);
					parameters->AddObject(new DoubleMatParam(SealevelchangeViscousPolarMotionEnum,viscouspolarmotion,3,nt));
					xDelete<IssmDouble>(viscouspolarmotion);
				}
			}
#ifdef _HAVE_AD_
			U_viscoelastic=xNew<IssmDouble>(M*ntimesteps,"t");
			if(horiz)H_viscoelastic=xNew<IssmDouble>(M*ntimesteps,"t");
#else
			U_viscoelastic=xNew<IssmDouble>(M*ntimesteps);
			if(horiz)H_viscoelastic=xNew<IssmDouble>(M*ntimesteps);
#endif
		}
		if(selfattraction){
			/*compute combined legendre + love number (elastic green function):*/
			m=DetermineLocalSize(M,IssmComm::GetComm());
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,IssmComm::GetComm());
#ifdef _HAVE_AD_
			G_gravi=xNew<IssmDouble>(M,"t");
			G_gravi_local=xNew<IssmDouble>(m,"t");
			G_viscoelastic=xNew<IssmDouble>(M*ntimesteps,"t");
			G_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps,"t");
#else
			G_gravi=xNew<IssmDouble>(M);
			G_gravi_local=xNew<IssmDouble>(m);
			G_viscoelastic=xNew<IssmDouble>(M*ntimesteps);
			G_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps);
#endif
		}
		if(viscous | elastic){
#ifdef _HAVE_AD_
			U_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps,"t");
			if(horiz)H_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps,"t");
#else
			U_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps);
			if(horiz)H_viscoelastic_local=xNew<IssmDouble>(m*ntimesteps);
#endif
		}

		if(rotation) parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.lovenumbers.tk2secular",TidalLoveK2SecularEnum));
		constant=3/rho_earth/planetarea;
		if(selfattraction){
			for(int i=lower_row;i<upper_row;i++){
				IssmDouble alpha,x;
				alpha= reCast<IssmDouble>(i)*degacc * M_PI / 180.0;
				G_gravi_local[i-lower_row]= constant*.5/sin(alpha/2.0);
			}
			if(viscous | elastic){
				for(int i=lower_row;i<upper_row;i++){
					IssmDouble alpha,x;
					alpha= reCast<IssmDouble>(i)*degacc * M_PI / 180.0;

					for(int t=0;t<ntimesteps;t++){
						G_viscoelastic_local[(i-lower_row)*ntimesteps+t]= (1.0+love_k[(ndeg-1)*precomputednt+t]-love_h[(ndeg-1)*precomputednt+t])*G_gravi_local[i-lower_row];
						U_viscoelastic_local[(i-lower_row)*ntimesteps+t]= (love_h[(ndeg-1)*precomputednt+t])*G_gravi_local[i-lower_row];
						if(horiz)H_viscoelastic_local[(i-lower_row)*ntimesteps+t]= 0; 
					}

					IssmDouble Pn = 0.; 
					IssmDouble Pn1 = 0.; 
					IssmDouble Pn2 = 0.; 
					IssmDouble Pn_p = 0.; 
					IssmDouble Pn_p1 = 0.; 
					IssmDouble Pn_p2 = 0.; 

					for (int n=0;n<ndeg;n++) {

						/*compute legendre polynomials: P_n(cos\theta) & d P_n(cos\theta)/ d\theta: */
						if(n==0){
							Pn=1; 
							Pn_p=0; 
						}
						else if(n==1){ 
							Pn = cos(alpha); 
							Pn_p = 1; 
						}
						else{
							Pn = ( (2*n-1)*cos(alpha)*Pn1 - (n-1)*Pn2 ) /n;
							Pn_p = ( (2*n-1)*(Pn1+cos(alpha)*Pn_p1) - (n-1)*Pn_p2 ) /n;
						}
						Pn2=Pn1; Pn1=Pn;
						Pn_p2=Pn_p1; Pn_p1=Pn_p;

						for(int t=0;t<ntimesteps;t++){
							IssmDouble deltalove_G;
							IssmDouble deltalove_U;

							deltalove_G = (love_k[n*precomputednt+t]-love_k[(ndeg-1)*precomputednt+t]-love_h[n*precomputednt+t]+love_h[(ndeg-1)*precomputednt+t]);
							deltalove_U = (love_h[n*precomputednt+t]-love_h[(ndeg-1)*precomputednt+t]);

							G_viscoelastic_local[(i-lower_row)*ntimesteps+t] += constant*deltalove_G*Pn;		                // gravitational potential 
							U_viscoelastic_local[(i-lower_row)*ntimesteps+t] += constant*deltalove_U*Pn;		                // vertical (up) displacement 
							if(horiz)H_viscoelastic_local[(i-lower_row)*ntimesteps+t] += constant*sin(alpha)*love_l[n*precomputednt+t]*Pn_p;		// horizontal displacements 
						}
					}
				}
			}
			else { //just copy G_gravi into G_viscoelastic
				for(int i=lower_row;i<upper_row;i++){
					for(int t=0;t<ntimesteps;t++){
						G_viscoelastic_local[(i-lower_row)*ntimesteps+t]= G_gravi_local[i-lower_row];
					}
				}
			}
			/*merge G_viscoelastic_local into G_viscoelastic; U_viscoelastic_local into U_viscoelastic; H_viscoelastic_local to H_viscoelastic:{{{*/
			recvcounts=xNew<int>(IssmComm::GetSize());
			displs=xNew<int>(IssmComm::GetSize());
			int  rc;
			int  offset;

			//deal with selfattraction first: 
			ISSM_MPI_Allgather(&m,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());
			/*displs: */
			ISSM_MPI_Allgather(&lower_row,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());
			/*All gather:*/
			ISSM_MPI_Allgatherv(G_gravi_local, m, ISSM_MPI_DOUBLE, G_gravi, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());

			rc=m*ntimesteps;
			offset=lower_row*ntimesteps;
			ISSM_MPI_Allgather(&rc,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgather(&offset,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());
			ISSM_MPI_Allgatherv(G_viscoelastic_local, m*ntimesteps, ISSM_MPI_DOUBLE, G_viscoelastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
			if(elastic){
				ISSM_MPI_Allgatherv(U_viscoelastic_local, m*ntimesteps, ISSM_MPI_DOUBLE, U_viscoelastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
				if(horiz)ISSM_MPI_Allgatherv(H_viscoelastic_local, m*ntimesteps, ISSM_MPI_DOUBLE, H_viscoelastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
			}

			/*free resources: */
			xDelete<int>(recvcounts);
			xDelete<int>(displs);

			/*Avoid singularity at 0: */
			G_gravi[0]=G_gravi[1];
			for(int t=0;t<ntimesteps;t++){
				G_viscoelastic[t]=G_viscoelastic[ntimesteps+t];
			}
			if(elastic){
				for(int t=0;t<ntimesteps;t++){
					U_viscoelastic[t]=U_viscoelastic[ntimesteps+t];
					if(horiz)H_viscoelastic[t]=H_viscoelastic[ntimesteps+t];
				}
			}

			/*Reinterpolate viscoelastic green kernels onto a regular gridded time 
			 *with steps equal to timeacc:*/
			if(viscous){
				nt=reCast<int,IssmDouble>((final_time-start_time)/timeacc)+1;
#ifdef _HAVE_AD_
				G_viscoelastic_interpolated=xNew<IssmDouble>(M*nt,"t");
				U_viscoelastic_interpolated=xNew<IssmDouble>(M*nt,"t");
				if(horiz) H_viscoelastic_interpolated=xNew<IssmDouble>(M*nt,"t");
				if(rotation){
					Pmtf_col_interpolated=xNew<IssmDouble>(nt,"t");
					Pmtf_ortho_interpolated=xNew<IssmDouble>(nt,"t");
					Pmtf_z_interpolated=xNew<IssmDouble>(nt,"t");
					Love_tk2_interpolated=xNew<IssmDouble>(nt,"t");
					Love_th2_interpolated=xNew<IssmDouble>(nt,"t");
					if (horiz) Love_tl2_interpolated=xNew<IssmDouble>(nt,"t");
				}
#else
				G_viscoelastic_interpolated=xNew<IssmDouble>(M*nt);
				U_viscoelastic_interpolated=xNew<IssmDouble>(M*nt);
				if(horiz) H_viscoelastic_interpolated=xNew<IssmDouble>(M*nt);
				if(rotation){
					Pmtf_col_interpolated=xNew<IssmDouble>(nt);
					Pmtf_ortho_interpolated=xNew<IssmDouble>(nt);
					Pmtf_z_interpolated=xNew<IssmDouble>(nt);
					Love_tk2_interpolated=xNew<IssmDouble>(nt);
					Love_th2_interpolated=xNew<IssmDouble>(nt);
					if (horiz) Love_tl2_interpolated=xNew<IssmDouble>(nt);
				}
#endif
				for(int t=0;t<nt;t++){
					IssmDouble lincoeff;
					IssmDouble viscoelastic_time=t*timeacc;
					int        timeindex2=-1;
					/*Find a way to interpolate precomputed Gkernels to our solution time stepping:*/
					if(t!=0){
						for(int t2=0;t2<ntimesteps;t2++){
							if (viscoelastic_time<love_timefreq[t2]){
								timeindex2=t2-1;
								if(timeindex2<0)_error_("Temporal Love numbers are computed  with a time accuracy superior to the requested solution time step!");
								lincoeff=(viscoelastic_time-love_timefreq[t2-1])/(love_timefreq[t2]-love_timefreq[t2-1]);
								break;
							}
						}
						if(timeindex2==-1)_error_("Temporal love numbers should be extended in time to encompass the requested solution time interval!");
					}
					else{
						timeindex2=0;
						lincoeff=0;
					}

					for(int index=0;index<M;index++){
						int timeindex=index*nt+t;
						int timepreindex= index*ntimesteps+timeindex2;
						G_viscoelastic_interpolated[timeindex]=(1-lincoeff)*G_viscoelastic[timepreindex]+lincoeff*G_viscoelastic[timepreindex+1];
						U_viscoelastic_interpolated[timeindex]=(1-lincoeff)*U_viscoelastic[timepreindex]+lincoeff*U_viscoelastic[timepreindex+1];
						if(horiz)H_viscoelastic_interpolated[timeindex]=(1-lincoeff)*H_viscoelastic[timepreindex]+lincoeff*H_viscoelastic[timepreindex+1];
					}

					if(rotation){
						int timepreindex= 2*precomputednt+timeindex2;
						Pmtf_col_interpolated[t]=(1.0-lincoeff)*love_pmtf_colinear[timeindex2]+lincoeff*love_pmtf_colinear[timeindex2+1];
						Pmtf_ortho_interpolated[t]=(1.0-lincoeff)*love_pmtf_ortho[timeindex2]+lincoeff*love_pmtf_ortho[timeindex2+1];
						Pmtf_z_interpolated[t]=1.0+(1.0-lincoeff)*love_k[timepreindex]+lincoeff*love_k[timepreindex+1];
						Love_tk2_interpolated[t]=(1.0-lincoeff)*love_tk[timepreindex]+lincoeff*love_tk[timepreindex+1];
						Love_th2_interpolated[t]=(1.0-lincoeff)*love_th[timepreindex]+lincoeff*love_th[timepreindex+1];
						if (horiz) Love_tl2_interpolated[t]=(1.0-lincoeff)*love_tl[timepreindex]+lincoeff*love_tl[timepreindex+1];
					}
				}
			}
			else {

				nt=1; //in elastic, or if we run only selfattraction, we need only one step
#ifdef _HAVE_AD_
				G_viscoelastic_interpolated=xNew<IssmDouble>(M,"t");
#else
				G_viscoelastic_interpolated=xNew<IssmDouble>(M);
#endif
				xMemCpy<IssmDouble>(G_viscoelastic_interpolated,G_viscoelastic,M);

				if(elastic){
#ifdef _HAVE_AD_
					U_viscoelastic_interpolated=xNew<IssmDouble>(M,"t");
					if (horiz) H_viscoelastic_interpolated=xNew<IssmDouble>(M,"t");
#else
					U_viscoelastic_interpolated=xNew<IssmDouble>(M);
					if (horiz) H_viscoelastic_interpolated=xNew<IssmDouble>(M);
#endif
					xMemCpy<IssmDouble>(U_viscoelastic_interpolated,U_viscoelastic,M);
					if (horiz) xMemCpy<IssmDouble>(H_viscoelastic_interpolated,H_viscoelastic,M);

					if(rotation){ //if this cpu handles degree 2
#ifdef _HAVE_AD_
						Pmtf_col_interpolated=xNew<IssmDouble>(1,"t");
						Pmtf_ortho_interpolated=xNew<IssmDouble>(1,"t");
						Pmtf_z_interpolated=xNew<IssmDouble>(1,"t");
						Love_tk2_interpolated=xNew<IssmDouble>(1,"t");
						Love_th2_interpolated=xNew<IssmDouble>(1,"t");
						if (horiz) Love_tl2_interpolated=xNew<IssmDouble>(1,"t");
#else
						Pmtf_col_interpolated=xNew<IssmDouble>(1);
						Pmtf_ortho_interpolated=xNew<IssmDouble>(1);
						Pmtf_z_interpolated=xNew<IssmDouble>(1);
						Love_tk2_interpolated=xNew<IssmDouble>(1);
						Love_th2_interpolated=xNew<IssmDouble>(1);
						if (horiz) Love_tl2_interpolated=xNew<IssmDouble>(1);
#endif

						Pmtf_col_interpolated[0]=love_pmtf_colinear[0];
						Pmtf_ortho_interpolated[0]=love_pmtf_ortho[0];
						Pmtf_z_interpolated[0]=1.0+love_k[2];
						Love_tk2_interpolated[0]=love_tk[2];
						Love_th2_interpolated[0]=love_th[2];
						if (horiz) Love_tl2_interpolated[0]=love_tl[2];
					}

				}
			}	
			/*Save our precomputed tables into parameters*/
			parameters->AddObject(new DoubleVecParam(SealevelchangeGSelfAttractionEnum,G_gravi,M));
			parameters->AddObject(new DoubleVecParam(SealevelchangeGViscoElasticEnum,G_viscoelastic_interpolated,M*nt));
			if(viscous || elastic){
				parameters->AddObject(new DoubleVecParam(SealevelchangeUViscoElasticEnum,U_viscoelastic_interpolated,M*nt));
				if(horiz)parameters->AddObject(new DoubleVecParam(SealevelchangeHViscoElasticEnum,H_viscoelastic_interpolated,M*nt));
				if(rotation){
					parameters->AddObject(new DoubleVecParam(SealevelchangePolarMotionTransferFunctionColinearEnum,Pmtf_col_interpolated,nt));
					parameters->AddObject(new DoubleVecParam(SealevelchangePolarMotionTransferFunctionOrthogonalEnum,Pmtf_ortho_interpolated,nt));
					parameters->AddObject(new DoubleVecParam(SealevelchangePolarMotionTransferFunctionZEnum,Pmtf_z_interpolated,nt));
					parameters->AddObject(new DoubleVecParam(SealevelchangeTidalH2Enum,Love_th2_interpolated,nt));
					parameters->AddObject(new DoubleVecParam(SealevelchangeTidalK2Enum,Love_tk2_interpolated,nt));
					if (horiz) parameters->AddObject(new DoubleVecParam(SealevelchangeTidalL2Enum,Love_tl2_interpolated,nt));
				}
			}
			/*free resources: */
			xDelete<IssmDouble>(G_gravi);
			xDelete<IssmDouble>(G_gravi_local);
			xDelete<IssmDouble>(G_viscoelastic);
			xDelete<IssmDouble>(G_viscoelastic_local);
			xDelete<IssmDouble>(G_viscoelastic_interpolated);
			if(elastic){
				xDelete<IssmDouble>(love_timefreq);
				xDelete<IssmDouble>(love_h);
				xDelete<IssmDouble>(love_k);
				xDelete<IssmDouble>(love_l);
				xDelete<IssmDouble>(love_th);
				xDelete<IssmDouble>(love_tk);
				xDelete<IssmDouble>(love_tl);

				xDelete<IssmDouble>(G_viscoelastic_interpolated);
				xDelete<IssmDouble>(U_viscoelastic);
				xDelete<IssmDouble>(U_viscoelastic_interpolated);
				xDelete<IssmDouble>(U_viscoelastic_local);
				xDelete<IssmDouble>(U_viscoelastic_interpolated);
				if(horiz){
					xDelete<IssmDouble>(H_viscoelastic);
					xDelete<IssmDouble>(H_viscoelastic_interpolated);
					xDelete<IssmDouble>(H_viscoelastic_local);
					xDelete<IssmDouble>(H_viscoelastic_interpolated);
				}
				if(rotation){
					xDelete<IssmDouble>(love_pmtf_colinear);
					xDelete<IssmDouble>(love_pmtf_ortho);
					xDelete<IssmDouble>(Pmtf_col_interpolated);
					xDelete<IssmDouble>(Pmtf_ortho_interpolated);
					xDelete<IssmDouble>(Pmtf_z_interpolated);
					xDelete<IssmDouble>(Love_tk2_interpolated);
					xDelete<IssmDouble>(Love_th2_interpolated);
					if (horiz) xDelete<IssmDouble>(Love_tl2_interpolated);
				}
			}
		} /*}}}*/

		/*Indicate we have not yet run the Geometry Core module: */
		parameters->AddObject(new BoolParam(SealevelchangeGeometryDoneEnum,false));
		/*}}}*/
	}

	/*Transitions:{{{ */
	iomodel->FetchData(&transitions,&transitions_M,&transitions_N,&ntransitions,"md.solidearth.transitions");
	if(transitions){
		parameters->AddObject(new DoubleMatArrayParam(SealevelchangeTransitionsEnum,transitions,ntransitions,transitions_M,transitions_N));

		for(int i=0;i<ntransitions;i++){
			IssmDouble* transition=transitions[i];
			xDelete<IssmDouble>(transition);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_M);
		xDelete<int>(transitions_N);
	} /*}}}*/
	/*Requested outputs {{{*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.solidearth.requested_outputs");
	if(numoutputs)parameters->AddObject(new StringArrayParam(SealevelchangeRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.solidearth.requested_outputs");
	/*}}}*/

}/*}}}*/

/*Finite Element Analysis*/
void           SealevelchangeAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           SealevelchangeAnalysis::PreCore(FemModel* femmodel){/*{{{*/

	int isuq=0;
	int modelid=0;

	/*Resolve Mmes using the modelid, if necessary: meaning if we are running a transient model and that UQ computations have not been triggered:*/
	femmodel->parameters->FindParam(&isuq,QmuIsdakotaEnum);
	if (!isuq && femmodel->inputs->GetInputObjectEnum(SolidearthExternalDisplacementEastRateEnum)==DatasetInputEnum){
		femmodel->parameters->FindParam(&modelid,SolidearthExternalModelidEnum);

		/*replace dataset of forcings with only one, the modelid'th:*/
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementNorthRateEnum, P1Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementEastRateEnum, P1Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalDisplacementUpRateEnum, P1Enum);
		MmeToInputFromIdx(femmodel->inputs,femmodel->elements,femmodel->parameters,modelid-1,SolidearthExternalGeoidRateEnum, P1Enum);

	}		

	/*run sea level change core geometry only once, after the Model Processor is done:*/
	sealevelchange_initialgeometry(femmodel);

}/*}}}*/
ElementVector* SealevelchangeAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* SealevelchangeAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* SealevelchangeAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* SealevelchangeAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           SealevelchangeAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           SealevelchangeAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           SealevelchangeAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemeneted yet!");

}/*}}}*/
void           SealevelchangeAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
