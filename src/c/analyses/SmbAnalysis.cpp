#include "./SmbAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

// FIX
#include "./shared/io/Print/Print.h"

/*Model processing*/
void SmbAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void SmbAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void SmbAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,SmbAnalysisEnum,P1Enum);
}/*}}}*/
int  SmbAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SmbAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    smb_model;
	bool   isdelta18o,ismungsm,isd18opd,issetpddfac,isprecipscaled,istemperaturescaled,isfirnwarming,isstochastic,ismappedforcing;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Figure out smb model: */
	iomodel->FindConstant(&smb_model,"md.smb.model");
	iomodel->FindConstant(&isstochastic,"md.stochasticforcing.isstochasticforcing");
	InputUpdateFromConstantx(inputs,elements,false,SmbIsInitializedEnum);
	switch(smb_model){
		case SMBforcingEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",SmbMassBalanceEnum,0.);
			if(isstochastic){
				iomodel->FetchDataToInput(inputs,elements,"md.stochasticforcing.default_id",StochasticForcingDefaultIdEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.mass_balance",BaselineSmbMassBalanceEnum,0.);
			}
			break;
		case SMBgembEnum:
			iomodel->FindConstant(&ismappedforcing,"md.smb.ismappedforcing");
			if (!ismappedforcing){
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Ta",SmbTaEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.V",SmbVEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.dswrf",SmbDswrfEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.dswdiffrf",SmbDswdiffrfEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.dlwrf",SmbDlwrfEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.P",SmbPEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.eAir",SmbEAirEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.pAir",SmbPAirEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Tmean",SmbTmeanEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Vmean",SmbVmeanEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.C",SmbCEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Tz",SmbTzEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Vz",SmbVzEnum);
			} else {
				iomodel->FetchDataToInput(inputs,elements,"md.smb.mappedforcingpoint",SmbMappedforcingpointEnum);
			}

			iomodel->FetchDataToInput(inputs,elements,"md.smb.zTop",SmbZTopEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dzTop",SmbDzTopEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dzMin",SmbDzMinEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.zY",SmbZYEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.zMax",SmbZMaxEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.zMin",SmbZMinEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Dzini",SmbDziniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Dini",SmbDiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Reini",SmbReiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Gdnini",SmbGdniniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Gspini",SmbGspiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.ECini",SmbECiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Wini",SmbWiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Aini",SmbAiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Adiffini",SmbAdiffiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Tini",SmbTiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.Sizeini",SmbSizeiniEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.aValue",SmbAValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dulwrfValue",SmbDulwrfValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.teValue",SmbTeValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.szaValue",SmbSzaValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.cotValue",SmbCotValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.ccsnowValue",SmbCcsnowValueEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.cciceValue",SmbCciceValueEnum);
			break;
		case SMBpddEnum:
			iomodel->FindConstant(&isdelta18o,"md.smb.isdelta18o");
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FindConstant(&issetpddfac,"md.smb.issetpddfac");
			iomodel->FetchDataToInput(inputs,elements,"md.thermal.spctemperature",ThermalSpctemperatureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0p",SmbS0pEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0t",SmbS0tEnum);
			if(isdelta18o || ismungsm){
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.temperatures_lgm",SmbTemperaturesLgmEnum);
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.temperatures_presentday",SmbTemperaturesPresentdayEnum);
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitations_presentday",SmbPrecipitationsPresentdayEnum);
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitations_lgm",SmbPrecipitationsLgmEnum);
			}else{
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitation",SmbPrecipitationEnum);
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.monthlytemperatures",SmbMonthlytemperaturesEnum);
			}
			if(issetpddfac){
				iomodel->FetchDataToInput(inputs,elements,"md.smb.pddfac_snow",SmbPddfacSnowEnum,-1.);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.pddfac_ice",SmbPddfacIceEnum,-1.);
			}
			break;
		case SMBpddSicopolisEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0p",SmbS0pEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0t",SmbS0tEnum);
			iomodel->FindConstant(&isfirnwarming,"md.smb.isfirnwarming");
			iomodel->FetchDataToInput(inputs,elements,"md.smb.smb_corr",SmbSmbCorrEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.precipitation_anomaly",SmbPrecipitationsAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.temperature_anomaly",SmbTemperaturesAnomalyEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.monthlytemperatures",SmbMonthlytemperaturesEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitation",SmbPrecipitationEnum);
			break;
		case SMBpddGCMEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.enhance_factor",SmbEnhanceFactorEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.lapserates",SmbGCMLapseratesEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.ref_surf",SmbGCMRefSurfaceEnum);
			break;
		case SMBd18opddEnum:
			iomodel->FindConstant(&istemperaturescaled,"md.smb.istemperaturescaled");
			iomodel->FindConstant(&isprecipscaled,"md.smb.isprecipscaled");
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FindConstant(&isd18opd,"md.smb.isd18opd");
			iomodel->FindConstant(&issetpddfac,"md.smb.issetpddfac");
			iomodel->FetchDataToInput(inputs,elements,"md.thermal.spctemperature",ThermalSpctemperatureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0p",SmbS0pEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0t",SmbS0tEnum);
			if(isd18opd){
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.temperatures_presentday",SmbTemperaturesPresentdayEnum);
				iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitations_presentday",SmbPrecipitationsPresentdayEnum);
				if(!istemperaturescaled){
					/*Fetch array*/
					IssmDouble* doublearray = NULL;
					int         M,N;
					iomodel->FetchData(&doublearray,&M,&N,"md.smb.temperatures_reconstructed");
					if(M!=iomodel->numberofvertices+1) _error_("md.smb.temperatures_reconstructed should have nbv+1 rows");
					if(N%12!=0) _error_("md.smb.temperatures_reconstructed should have a multiple of 12 columns (since it is monthly)");

					/*Build indices*/
					int* ids = xNew<int>(N); for(int i=0;i<N;i++) ids[i] = i;

					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputCreate(doublearray,M-1,N,ids,N,inputs,iomodel,SmbTemperaturesReconstructedEnum);
					}
					xDelete<int>(ids);
					iomodel->DeleteData(doublearray,"md.smb.temperatures_reconstructed");
				}
				if(!isprecipscaled){
					/*Fetch array*/
					IssmDouble* doublearray = NULL;
					int         M,N;
					iomodel->FetchData(&doublearray,&M,&N,"md.smb.precipitations_reconstructed");
					if(M!=iomodel->numberofvertices+1) _error_("md.smb.precipitations_reconstructed should have nbv+1 rows");
					if(N%12!=0) _error_("md.smb.precipitations_reconstructed should have a multiple of 12 columns (since it is monthly)");

					/*Build indices*/
					int* ids = xNew<int>(N); for(int i=0;i<N;i++) ids[i] = i;

					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputCreate(doublearray,M-1,N,ids,N,inputs,iomodel,SmbPrecipitationsReconstructedEnum);
					}
					xDelete<int>(ids);
					iomodel->DeleteData(doublearray,"md.smb.precipitations_reconstructed");
				}
			}
			if(issetpddfac){
				iomodel->FetchDataToInput(inputs,elements,"md.smb.pddfac_snow",SmbPddfacSnowEnum,-1.);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.pddfac_ice",SmbPddfacIceEnum,-1.);
			}
			break;
		case SMBgradientsEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.href",SmbHrefEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.smbref",SmbSmbrefEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_pos",SmbBPosEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_neg",SmbBNegEnum);
			break;
		case SMBgradientselaEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.ela",SmbElaEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_pos",SmbBPosEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_neg",SmbBNegEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_max",SmbBMaxEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.b_min",SmbBMinEnum);
			break;
		case SMBarmaEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.basin_id",SmbBasinsIdEnum);
			break;
		case SMBhenningEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.smbref",SmbSmbrefEnum,0.);
			break;
		case SMBcomponentsEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.accumulation",SmbAccumulationEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.evaporation",SmbEvaporationEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.runoff",SmbRunoffEnum,0.);
			break;
		case SMBmeltcomponentsEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.smb.accumulation",SmbAccumulationEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.evaporation",SmbEvaporationEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.melt",SmbMeltEnum,0.);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.refreeze",SmbRefreezeEnum,0.);
			break;
		case SMBgradientscomponentsEnum:
			/* Nothing to add to input */
			break;
		case SMBsemicEnum:
			int ismethod;
			//if(VerboseSolution()) _printf0_("   smb semic: UpdateElements.\n");
			//iomodel->FetchDataToInput(inputs,elements,"md.thermal.spctemperature",ThermalSpctemperatureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0gcm",SmbS0gcmEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailysnowfall",SmbDailysnowfallEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailyrainfall",SmbDailyrainfallEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailydsradiation",SmbDailydsradiationEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailydlradiation",SmbDailydlradiationEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailywindspeed",SmbDailywindspeedEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailypressure",SmbDailypressureEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailyairdensity",SmbDailyairdensityEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailyairhumidity",SmbDailyairhumidityEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dailytemperature",SmbDailytemperatureEnum);
			// assign initial SEMIC temperature from initialization class.
			if(VerboseSolution()) _printf0_("   smb semic: UpdateElements - temperature.\n");
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.temperature",TemperatureEnum);

			iomodel->FindConstant(&ismethod,"md.smb.ismethod");
			if (ismethod == 1){
				if(VerboseSolution()) _printf0_("   smb semic: UpdateElements - albedo.\n");
				iomodel->FetchDataToInput(inputs,elements,"md.smb.albedo",SmbAlbedoInitEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.albedo_snow",SmbAlbedoSnowInitEnum);
				if(VerboseSolution()) _printf0_("   smb semic: UpdateElements - Hice/Hsnow.\n");
				iomodel->FetchDataToInput(inputs,elements,"md.smb.hice",SmbHIceInitEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.hsnow",SmbHSnowInitEnum);

				// initial Temperature amplitude.
				if(VerboseSolution()) _printf0_("   smb semic: UpdateElements - Tamp.\n");
				iomodel->FetchDataToInput(inputs,elements,"md.smb.Tamp",SmbTampEnum);

				// assign masking 
				iomodel->FetchDataToInput(inputs,elements,"md.smb.mask",SmbMaskEnum);
				iomodel->FetchDataToInput(inputs,elements,"md.smb.qmr",SmbSemicQmrInitEnum);
				if(VerboseSolution()) _printf0_("   smb semic: UpdateElements - done.\n");
			}
			break;
		case SMBdebrisEvattEnum:
			iomodel->FetchDataToInput(inputs,elements,"md.initialization.debris",DebrisThicknessEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.s0t",SmbS0tEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.snowheight",SmbSnowheightEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.temperature",SmbMonthlytemperaturesEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.precipitation",SmbPrecipitationEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.dsradiation",SmbMonthlydsradiationEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.dlradiation",SmbMonthlydlradiationEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.windspeed",SmbMonthlywindspeedEnum);
			iomodel->FetchDataToDatasetInput(inputs,elements,"md.smb.airhumidity",SmbMonthlyairhumidityEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.precipitation_anomaly",SmbPrecipitationsAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.temperature_anomaly",SmbTemperaturesAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dsradiation_anomaly",SmbDsradiationAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.dlradiation_anomaly",SmbDlradiationAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.windspeed_anomaly",SmbWindspeedAnomalyEnum);
			iomodel->FetchDataToInput(inputs,elements,"md.smb.airhumidity_anomaly",SmbAirhumidityAnomalyEnum);
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

}/*}}}*/
void SmbAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;
	bool    isdelta18o,ismungsm,isd18opd,issetpddfac,interp,cycle,isfirnwarming,ismappedforcing;
	int     smb_model, smbslices, averaging;
	IssmDouble *temp = NULL;
	int         N,M,Nt,Nx,Ny;

	parameters->AddObject(iomodel->CopyConstantObject("md.smb.model",SmbEnum));

	iomodel->FindConstant(&smb_model,"md.smb.model");
	iomodel->FindConstant(&interp,"md.timestepping.interp_forcing");
	iomodel->FindConstant(&cycle,"md.timestepping.cycle_forcing");

	iomodel->FindConstant(&smbslices,"md.smb.steps_per_step");
	parameters->AddObject(new IntParam(SmbStepsPerStepEnum,smbslices));

	parameters->AddObject(iomodel->CopyConstantObject("md.smb.averaging",SmbAveragingEnum));

	switch(smb_model){
		case SMBforcingEnum:
		case SMBgradientsEnum:
		case SMBgradientselaEnum:
		case SMBhenningEnum:
		case SMBcomponentsEnum:
		case SMBmeltcomponentsEnum:
			break;
			//case SMBarmaEnum:
		case SMBarmaEnum:
			/*Add parameters that are not in standard nbvertices format*/
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.num_basins",SmbNumBasinsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.num_params",SmbNumParamsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.num_breaks",SmbNumBreaksEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ar_order",SmbARMAarOrderEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ma_order",SmbARMAmaOrderEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.arma_timestep",SmbARMATimestepEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.num_bins",SmbNumElevationBinsEnum));
			iomodel->FetchData(&temp,&M,&N,"md.smb.datebreaks");
			parameters->AddObject(new DoubleMatParam(SmbARMAdatebreaksEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.polynomialparams");
			parameters->AddObject(new DoubleMatParam(SmbARMApolyparamsEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.arlag_coefs");
			parameters->AddObject(new DoubleMatParam(SmbARMAarlagcoefsEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.malag_coefs");
			parameters->AddObject(new DoubleMatParam(SmbARMAmalagcoefsEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.lapserates");
			parameters->AddObject(new DoubleMatParam(SmbLapseRatesEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.elevationbins");
			parameters->AddObject(new DoubleMatParam(SmbElevationBinsEnum,temp,M,N));
			xDelete<IssmDouble>(temp);
			iomodel->FetchData(&temp,&M,&N,"md.smb.refelevation");
			parameters->AddObject(new DoubleVecParam(SmbRefElevationEnum,temp,N));
			xDelete<IssmDouble>(temp);
			break;
		case SMBpddSicopolisEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isfirnwarming",SmbIsfirnwarmingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.pdd_fac_ice",PddfacIceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.pdd_fac_snow",PddfacSnowEnum));
			break;
		case SMBdebrisEvattEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.qlaps",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dsgrad",SmbSWgradEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dlgrad",SmbLWgradEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.windspeedgrad",SmbWindspeedgradEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.humiditygrad",SmbHumiditygradEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.icealbedo",SmbIcealbedoEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.snowalbedo",SmbSnowalbedoEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.debrisalbedo",SmbDebrisalbedoEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isAnderson",SmbDebrisIsAndersonEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.iscryokarst",SmbDebrisIsCryokarstEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.AndersonD0",SmbDebrisAndersonD0Enum));
			break;
		case SMBgembEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aIce",SmbAIceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aSnow",SmbASnowEnum));
			iomodel->FindConstant(&ismappedforcing,"md.smb.ismappedforcing");
			if (ismappedforcing){
				iomodel->FetchData(&temp,&M,&N,"md.smb.Ta"); _assert_(M>=1 && N>=1); 
				parameters->AddObject(new TransientArrayParam(SmbTaParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.Ta");
				iomodel->FetchData(&temp,&M,&N,"md.smb.V"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbVParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.V");
				iomodel->FetchData(&temp,&M,&N,"md.smb.dswrf"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbDswrfParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.dswrf");
				iomodel->FetchData(&temp,&M,&N,"md.smb.dswdiffrf"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbDswdiffrfParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.dswdiffrf");
				iomodel->FetchData(&temp,&M,&N,"md.smb.dlwrf"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbDlwrfParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.dlwrf");
				iomodel->FetchData(&temp,&M,&N,"md.smb.P"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbPParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.P");
				iomodel->FetchData(&temp,&M,&N,"md.smb.eAir"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbEAirParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.eAir");
				iomodel->FetchData(&temp,&M,&N,"md.smb.pAir"); _assert_(M>=1 && N>=1);
				parameters->AddObject(new TransientArrayParam(SmbPAirParamEnum,temp,&temp[N*(M-1)],interp,cycle,N,M));
				iomodel->DeleteData(temp,"md.smb.pAir");

				iomodel->FetchData(&temp,&M,&N,"md.smb.Tmean"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbTmeanParamEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.Tmean");
				iomodel->FetchData(&temp,&M,&N,"md.smb.Vmean"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbVmeanParamEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.Vmean");
				iomodel->FetchData(&temp,&M,&N,"md.smb.C"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbCParamEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.C");
				iomodel->FetchData(&temp,&M,&N,"md.smb.Tz"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbTzParamEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.Tz");
				iomodel->FetchData(&temp,&M,&N,"md.smb.Vz"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbVzParamEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.Vz");
				iomodel->FetchData(&temp,&M,&N,"md.smb.mappedforcingelevation"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbMappedforcingelevationEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.mappedforcingelevation");

				iomodel->FetchData(&temp,&M,&N,"md.smb.lapseTaValue"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbLapseTaValueEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.lapseTaValue");
				iomodel->FetchData(&temp,&M,&N,"md.smb.lapsedlwrfValue"); _assert_(N==1);
				parameters->AddObject(new DoubleVecParam(SmbLapsedlwrfValueEnum,&temp[0],M));
				iomodel->DeleteData(temp,"md.smb.lapsedlwrfValue");

			}
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aIdx",SmbAIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.eIdx",SmbEIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.tcIdx",SmbTcIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.swIdx",SmbSwIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.denIdx",SmbDenIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dsnowIdx",SmbDsnowIdxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.cldFrac",SmbCldFracEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.t0wet",SmbT0wetEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.t0dry",SmbT0dryEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.K",SmbKEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aSnow",SmbASnowEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.aIce",SmbAIceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dt",SmbDtEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isgraingrowth",SmbIsgraingrowthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isalbedo",SmbIsalbedoEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isshortwave",SmbIsshortwaveEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isthermal",SmbIsthermalEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isaccumulation",SmbIsaccumulationEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismelt",SmbIsmeltEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdensification",SmbIsdensificationEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isturbulentflux",SmbIsturbulentfluxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isconstrainsurfaceT",SmbIsconstrainsurfaceTEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdeltaLWup",SmbIsdeltaLWupEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismappedforcing",SmbIsmappedforcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isprecipforcingremapped",SmbIsprecipforcingremappedEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.InitDensityScaling",SmbInitDensityScalingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ThermoDeltaTScaling",SmbThermoDeltaTScalingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.adThresh",SmbAdThreshEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.teThresh",SmbTeThreshEnum));
			break;
		case SMBpddEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlapslgm",SmbRlapslgmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdelta18o",SmbIsdelta18oEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismungsm",SmbIsmungsmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.issetpddfac",SmbIssetpddfacEnum));
			iomodel->FindConstant(&isdelta18o,"md.smb.isdelta18o");
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");

			if(ismungsm){
				iomodel->FetchData(&temp,&N,&M,"md.smb.Pfac"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbPfacEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.Pfac");

				iomodel->FetchData(&temp,&N,&M,"md.smb.Tdiff"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbTdiffEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.Tdiff");

				iomodel->FetchData(&temp,&N,&M,"md.smb.sealev"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbSealevEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.sealev");
			}
			if(isdelta18o){
				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.delta18o");

				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o_surface"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oSurfaceEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.delta18o_surface");
			}

			break;
		case SMBpddGCMEnum:
			iomodel->FetchData(&temp,&M,&N,"md.smb.x_grid"); _assert_(N==1); Nx = M;
			parameters->AddObject(new DoubleVecParam(SmbGCMXgridEnum,&temp[0],M));
			iomodel->DeleteData(temp,"md.smb.lat");
			iomodel->FetchData(&temp,&M,&N,"md.smb.y_grid"); _assert_(N==1); Ny = M;
			parameters->AddObject(new DoubleVecParam(SmbGCMYgridEnum,&temp[0],M));
			iomodel->DeleteData(temp,"md.smb.lon");
			iomodel->FetchData(&temp,&M,&N,"md.smb.precipitation"); _assert_(N>=1 && M==Nx*Ny+1);
			parameters->AddObject(new TransientGriddedFieldParam(SmbGCMPrecipitationEnum,temp,&temp[N*(M-1)],interp,cycle,Nx,Ny,N));
			iomodel->DeleteData(temp,"md.smb.precipitation");
			iomodel->FetchData(&temp,&M,&N,"md.smb.temperature"); _assert_(N>=1 && M==Nx*Ny+1);
			parameters->AddObject(new TransientGriddedFieldParam(SmbGCMTemperatureEnum,temp,&temp[N*(M-1)],interp,cycle,Nx,Ny,N));
			iomodel->DeleteData(temp,"md.smb.temperature");

			parameters->AddObject(iomodel->CopyConstantObject("md.smb.allsolidtemperature",SmbAllSolidTempEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.allliquidtemperature",SmbAllLiquidTempEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ddf_snow",SmbDdfSnowEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ddf_ice",SmbDdfIceEnum));

			break;
		case SMBd18opddEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlapslgm",SmbRlapslgmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.dpermil",SmbDpermilEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismungsm",SmbIsmungsmEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.isd18opd",SmbIsd18opdEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.issetpddfac",SmbIssetpddfacEnum));
			iomodel->FindConstant(&ismungsm,"md.smb.ismungsm");
			iomodel->FindConstant(&isd18opd,"md.smb.isd18opd");
			iomodel->FindConstant(&issetpddfac,"md.smb.issetpddfac");
			if(isd18opd){
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.f",SmbFEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.istemperaturescaled",SmbIstemperaturescaledEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.isprecipscaled",SmbIsprecipscaledEnum));
				iomodel->FetchData(&temp,&N,&M,"md.smb.delta18o"); _assert_(N==2);
				parameters->AddObject(new TransientParam(SmbDelta18oEnum,&temp[0],&temp[M],interp,cycle,M));
				iomodel->DeleteData(temp,"md.smb.delta18o");

				IssmDouble yts;
				bool istemperaturescaled,isprecipscaled;
				iomodel->FindConstant(&yts,"md.constants.yts");
				iomodel->FindConstant(&istemperaturescaled,"md.smb.istemperaturescaled");
				iomodel->FindConstant(&isprecipscaled,"md.smb.isprecipscaled");
				if(!istemperaturescaled){
					/*Fetch array*/
					IssmDouble* doublearray = NULL;
					int         M,N;
					iomodel->FetchData(&doublearray,&M,&N,"md.smb.temperatures_reconstructed");
					if(M!=iomodel->numberofvertices+1) _error_("md.smb.temperatures_reconstructed should have nbv+1 rows");
					if(N%12!=0) _error_("md.smb.temperatures_reconstructed should have a multiple of 12 columns (since it is monthly)");
					int numyears = N/12; _assert_(numyears*12==N);

					/*Check times*/
					#ifdef _ISSM_DEBUG_
					for(int i=0;i<numyears;i++){
						for(int j=1;j<12;j++){
							//_assert_(floor(doublearray[(M-1)*N+i*12+j]/yts)==floor(doublearray[(M-1)*N+i*12]/yts));
							_assert_(doublearray[(M-1)*N+i*12+j]>doublearray[(M-1)*N+i*12+j-1]);
						}
					}
					#endif

					/*Build time*/
					IssmDouble* times = xNew<IssmDouble>(numyears); for(int i=0;i<numyears;i++) times[i] = doublearray[(M-1)*N+i*12];
					parameters->AddObject(new DoubleVecParam(SmbTemperaturesReconstructedYearsEnum,times,numyears));
					xDelete<IssmDouble>(times);
					iomodel->DeleteData(doublearray,"md.smb.temperatures_reconstructed");
				}
				if(!isprecipscaled){
					/*Fetch array*/
					IssmDouble* doublearray = NULL;
					int         M,N;
					iomodel->FetchData(&doublearray,&M,&N,"md.smb.precipitations_reconstructed");
					if(M!=iomodel->numberofvertices+1) _error_("md.smb.precipitations_reconstructed should have nbv+1 rows");
					if(N%12!=0) _error_("md.smb.precipitations_reconstructed should have a multiple of 12 columns (since it is monthly)");
					int numyears = N/12; _assert_(numyears*12==N);

					/*Check times*/
					#ifdef _ISSM_DEBUG_
					for(int i=0;i<numyears;i++){
						for(int j=1;j<12;j++){
							//_assert_(floor(doublearray[(M-1)*N+i*12+j]/yts)==floor(doublearray[(M-1)*N+i*12]/yts));
							_assert_(doublearray[(M-1)*N+i*12+j]>doublearray[(M-1)*N+i*12+j-1]);
						}
					}
					#endif

					/*Build time*/
					IssmDouble* times = xNew<IssmDouble>(numyears); for(int i=0;i<numyears;i++) times[i] = doublearray[(M-1)*N+i*12];
					parameters->AddObject(new DoubleVecParam(SmbPrecipitationsReconstructedYearsEnum,times,numyears));
					xDelete<IssmDouble>(times);
					iomodel->DeleteData(doublearray,"md.smb.precipitations_reconstructed");
				}
			}

			break;
		case SMBgradientscomponentsEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.accualti",SmbAccualtiEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.runoffalti",SmbRunoffaltiEnum));

			iomodel->FetchData(&temp,&N,&M,"md.smb.accugrad"); _assert_(N==2);
			parameters->AddObject(new TransientParam(SmbAccugradEnum,&temp[0],&temp[M],interp,cycle,M));
			iomodel->DeleteData(temp,"md.smb.accugrad");
			iomodel->FetchData(&temp,&N,&M,"md.smb.runoffgrad"); _assert_(N==2);
			parameters->AddObject(new TransientParam(SmbRunoffgradEnum,&temp[0],&temp[M],interp,cycle,M));
			iomodel->DeleteData(temp,"md.smb.runoffgrad");

			iomodel->FetchData(&temp,&N,&M,"md.smb.accuref"); _assert_(N==2);
			parameters->AddObject(new TransientParam(SmbAccurefEnum,&temp[0],&temp[M],interp,cycle,M));
			iomodel->DeleteData(temp,"md.smb.accuref");
			iomodel->FetchData(&temp,&N,&M,"md.smb.runoffref"); _assert_(N==2);
			parameters->AddObject(new TransientParam(SmbRunoffrefEnum,&temp[0],&temp[M],interp,cycle,M));
			iomodel->DeleteData(temp,"md.smb.runoffref");
			break;
		case SMBsemicEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfac",SmbDesfacEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rlaps",SmbRlapsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.rdl",SmbRdlEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.smb.ismethod",SmbSemicMethodEnum));
			int ismethod;
			parameters->FindParam(&ismethod,SmbSemicMethodEnum);
			if (ismethod == 1){
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.desfacElevation",SmbDesfacElevEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.hcrit",SmbSemicHcritEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.rcrit",SmbSemicRcritEnum));
				/*Define albedo parameters.*/
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.albedo_scheme",SmbAlbedoSchemeEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.alb_smax",SmbAlbedoSnowMaxEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.alb_smin",SmbAlbedoSnowMinEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.albi",SmbAlbedoIceEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.albl",SmbAlbedoLandEnum));

				//albedo parameter - slatter
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.tmax",SmbSemicTmaxEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.tmin",SmbSemicTminEnum));

				//albedo parameter - isba & denby
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.mcrit",SmbSemicMcritEnum));
				//albedo parameter - isba
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.tau_a",SmbSemicTauAEnum)); 
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.tau_f",SmbSemicTauFEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.wcrit",SmbSemicWcritEnum));
				//albedo parameter - alex
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.tmid",SmbSemicTmidEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.afac",SmbSemicAfacEnum));

				/* Set specific options*/
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.isdesertification",SmbSemicIsDesertificationEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.smb.isLWDcorrect",SmbSemicIsLWDcorrectEnum));
			}
			/*Nothing to add to parameters*/
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.smb.requested_outputs");
	parameters->AddObject(new IntParam(SmbNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SmbRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.smb.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           SmbAnalysis::Core(FemModel* femmodel){/*{{{*/

	int    smb_model;

	/*Figure out smb model: */
	femmodel->parameters->FindParam(&smb_model,SmbEnum);

	/*branch to correct module*/
	switch(smb_model){
		case SMBforcingEnum:
			SmbForcingx(femmodel);
			break;
		case SMBgembEnum:
			Gembx(femmodel);
			break;
		case SMBpddEnum:
			bool isdelta18o,ismungsm;
			femmodel->parameters->FindParam(&isdelta18o,SmbIsdelta18oEnum);
			femmodel->parameters->FindParam(&ismungsm,SmbIsmungsmEnum);
			if(isdelta18o){
				if(VerboseSolution()) _printf0_("   call Delta18oParameterization module\n");
				Delta18oParameterizationx(femmodel);
			}
			if(ismungsm){
				if(VerboseSolution()) _printf0_("   call MungsmtpParameterization module\n");
				MungsmtpParameterizationx(femmodel);
			}
			if(VerboseSolution()) _printf0_("   call positive degree day module\n");
			PositiveDegreeDayx(femmodel);
			break;
		case SMBpddSicopolisEnum:
			if(VerboseSolution()) _printf0_("   call SICOPOLIS positive degree day module\n");
			PositiveDegreeDaySicopolisx(femmodel);
			break;
		case SMBpddGCMEnum:
			if(VerboseSolution()) _printf0_("   call positive degree day module based on downsacling GCM data\n");
			PositiveDegreeDayGCMx(femmodel);
			break;
		case SMBd18opddEnum:
			bool isd18opd;
			femmodel->parameters->FindParam(&isd18opd,SmbIsd18opdEnum);
			if(isd18opd){
				if(VerboseSolution()) _printf0_("   call Delta18opdParameterization module\n");
				Delta18opdParameterizationx(femmodel);
				if(VerboseSolution()) _printf0_("   call positive degree day module\n");
				PositiveDegreeDayx(femmodel);
			}
			break;
		case SMBgradientsEnum:
			if(VerboseSolution())_printf0_("	call smb gradients module\n");
			SmbGradientsx(femmodel);
			break;
		case SMBgradientselaEnum:
			if(VerboseSolution())_printf0_("	call smb gradients ela module\n");
			SmbGradientsElax(femmodel);
			break;
		case SMBarmaEnum:
			if(VerboseSolution())_printf0_("    call smb arma module\n");
			Smbarmax(femmodel);
			break;
		case SMBhenningEnum:
			if(VerboseSolution())_printf0_("  call smb Henning module\n");
			SmbHenningx(femmodel);
			break;
		case SMBcomponentsEnum:
			if(VerboseSolution())_printf0_("  call smb Components module\n");
			SmbComponentsx(femmodel);
			break;
		case SMBmeltcomponentsEnum:
			if(VerboseSolution())_printf0_("  call smb Melt Components module\n");
			SmbMeltComponentsx(femmodel);
			break;
		case SMBgcmEnum:
			/*Nothing to be done*/
			break;
		case SMBgradientscomponentsEnum:
			if(VerboseSolution())_printf0_("	call smb gradients components module\n");
			SmbGradientsComponentsx(femmodel);
			break;
		case SMBsemicEnum:
			#ifdef _HAVE_SEMIC_
			if(VerboseSolution())_printf0_("   call smb SEMIC module\n");
			int ismethod;
			femmodel->parameters->FindParam(&ismethod,SmbSemicMethodEnum);
			SmbSemicx(femmodel,ismethod);
			#else
			_error_("SEMIC not installed");
			#endif //_HAVE_SEMIC_
			break;
		case SMBdebrisEvattEnum:
			if(VerboseSolution())_printf0_("        call smb Evatt debris module\n");
			SmbDebrisEvattx(femmodel);
			break;
		default:
			_error_("Surface mass balance model "<<EnumToStringx(smb_model)<<" not supported yet");
	}

}/*}}}*/
void           SmbAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* SmbAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementMatrix* SmbAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* SmbAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* SmbAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           SmbAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void           SmbAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
