/*!\file: CreateParameters.cpp
 * \brief general driver for creating parameters dataset
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "../ParseToolkitsOptionsx/ParseToolkitsOptionsx.h"
#include "./ModelProcessorx.h"

void CreateParameters(Parameters* parameters,IoModel* iomodel,char* rootpath,FILE* toolkitsoptionsfid,const int solution_type){

	int         i,j,m,k;
	int         numoutputs,basalforcing_model,timestepping_type;
	char**      requestedoutputs = NULL;
	char**      outputonnodes = NULL;
	char*       fieldname = NULL;
	IssmDouble  time;

	/*parameters for mass flux:*/
	int          mass_flux_num_profiles     = 0;
	bool         qmu_mass_flux_present      = false;
	bool         autodiff_mass_flux_present = false;
	bool         mass_flux_present          = false;
	bool         interp,cycle;
	IssmDouble **array                      = NULL;
	int         *mdims_array                = NULL;
	int         *ndims_array                = NULL;
	IssmDouble  *temp_matrix                = NULL;
	int          temp_m,temp_n;
	IssmDouble  *matrix                     = NULL;
	int          count;

	IssmDouble *temp = NULL;
	IssmDouble *transparam = NULL;
	IssmDouble  yts;
	int         N,M;

	/*Copy some constants from iomodel */
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.domain_type",DomainTypeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.domain_dimension",DomainDimensionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.output_frequency",SettingsOutputFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.sb_coupling_frequency",SettingsSbCouplingFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.checkpoint_frequency",SettingsCheckpointFrequencyEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.constants.yts",ConstantsYtsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.debug.profiling",DebugProfilingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.average_vertex_connectivity",MeshAverageVertexConnectivityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.waitonlock",SettingsWaitonlockEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberofvertices",MeshNumberofverticesEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberofelements",MeshNumberofelementsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.io_gather",SettingsIoGatherEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.settings.solver_residue_threshold",SettingsSolverResidueThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.isautodiff",AutodiffIsautodiffEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.qmu.isdakota",QmuIsdakotaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.inversion.iscontrol",InversionIscontrolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.inversion.type",InversionTypeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.calving.law",CalvingLawEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.frontalforcings.parameterization",FrontalForcingsParamEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.isgrd",SolidearthSettingsGRDEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.grdmodel",GrdModelEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.runfrequency",SolidearthSettingsRunFrequencyEnum));
	parameters->AddObject(new IntParam(SealevelchangeRunCountEnum,1));

	  {/*This is specific to ice...*/
		parameters->AddObject(iomodel->CopyConstantObject("md.mesh.elementtype",MeshElementtypeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.steadystate.reltol",SteadystateReltolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.steadystate.maxiter",SteadystateMaxiterEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.migration",GroundinglineMigrationEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.friction_interpolation",GroundinglineFrictionInterpolationEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.melt_interpolation",GroundinglineMeltInterpolationEnum));
		//parameters->AddObject(iomodel->CopyConstantObject("md.groundingline.intrusion_distance",GroundinglineIntrusionDistanceEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isstressbalance",TransientIsstressbalanceEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ismasstransport",TransientIsmasstransportEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ismmemasstransport",TransientIsmmemasstransportEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isoceantransport",TransientIsoceantransportEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isage",TransientIsageEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.issmb",TransientIssmbEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isthermal",TransientIsthermalEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isgroundingline",TransientIsgroundinglineEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isesa",TransientIsesaEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isdamageevolution",TransientIsdamageevolutionEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ishydrology",TransientIshydrologyEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.ismovingfront",TransientIsmovingfrontEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isslc",TransientIsslcEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isoceancoupling",TransientIsoceancouplingEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.amr_frequency",TransientAmrFrequencyEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.isdebris",TransientIsdebrisEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.transient.issampling",TransientIssamplingEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.isstochasticforcing",StochasticForcingIsStochasticForcingEnum));

		/*For stress balance only*/
		parameters->AddObject(iomodel->CopyConstantObject("md.flowequation.isFS",FlowequationIsFSEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.rift_penalty_threshold",StressbalanceRiftPenaltyThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.maxiter",StressbalanceMaxiterEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.restol",StressbalanceRestolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.reltol",StressbalanceReltolEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.stressbalance.abstol",StressbalanceAbstolEnum));
		if(iomodel->domaintype==Domain3DEnum)
		 parameters->AddObject(iomodel->CopyConstantObject("md.mesh.numberoflayers",MeshNumberoflayersEnum));
	  }

	/*amr properties*/
	int amrtype,amr_frequency;
	iomodel->FindConstant(&amr_frequency,"md.transient.amr_frequency");
	if(solution_type==TransientSolutionEnum && amr_frequency){
		/*Load common amr parameters*/
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.type",AmrTypeEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.groundingline_distance",AmrGroundingLineDistanceEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.icefront_distance",AmrIceFrontDistanceEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.thicknesserror_threshold",AmrThicknessErrorThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.thicknesserror_groupthreshold",AmrThicknessErrorGroupThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.thicknesserror_maximum",AmrThicknessErrorMaximumEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.deviatoricerror_threshold",AmrDeviatoricErrorThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.deviatoricerror_groupthreshold",AmrDeviatoricErrorGroupThresholdEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.deviatoricerror_maximum",AmrDeviatoricErrorMaximumEnum));
		parameters->AddObject(iomodel->CopyConstantObject("md.amr.restart",AmrRestartEnum));
		/*Load specific amr parameters*/
		iomodel->FindConstant(&amrtype,"md.amr.type");
		switch(amrtype){
			#ifdef _HAVE_NEOPZ_
			case AmrNeopzEnum:
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.level_max",AmrLevelMaxEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.gradation",AmrGradationEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.lag",AmrLagEnum));
				break;
			#endif

			#ifdef _HAVE_BAMG_
			case AmrBamgEnum:
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.hmin",AmrHminEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.hmax",AmrHmaxEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.err",AmrErrEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.keepmetric",AmrKeepMetricEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.gradation",AmrGradationEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.groundingline_resolution",AmrGroundingLineResolutionEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.icefront_resolution",AmrIceFrontResolutionEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.thicknesserror_resolution",AmrThicknessErrorResolutionEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.amr.deviatoricerror_resolution",AmrDeviatoricErrorResolutionEnum));
				/*Convert fieldname to enum and put it in params*/
				iomodel->FindConstant(&fieldname,"md.amr.fieldname");
				parameters->AddObject(new IntParam(AmrFieldEnum,StringToEnumx(fieldname)));
				xDelete<char>(fieldname);
				break;
			#endif

			default:
				_error_("Adaptive mesh refinement "<<EnumToStringx(amrtype)<<" not implemented yet");
		}
	}

	/*Basal forcing parameters*/
	parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.model",BasalforcingsEnum));
	iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
	switch(basalforcing_model){
		case FloatingMeltRateEnum:
			/*Nothing to add to parameters*/
			break;
		case LinearFloatingMeltRateEnum:
			iomodel->FindConstant(&interp,"md.timestepping.interp_forcing");
			iomodel->FindConstant(&cycle,"md.timestepping.cycle_forcing");
			iomodel->FetchData(&transparam,&N,&M,"md.basalforcings.deepwater_melting_rate");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(BasalforcingsDeepwaterMeltingRateEnum,transparam[0]));
			}
			else{
				_assert_(N==2);
				parameters->AddObject(new TransientParam(BasalforcingsDeepwaterMeltingRateEnum,&transparam[0],&transparam[M],interp,cycle,M));
			}
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&N,&M,"md.basalforcings.upperwater_melting_rate");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(BasalforcingsUpperwaterMeltingRateEnum,transparam[0]));
			}
			else{
				_assert_(N==2);
				parameters->AddObject(new TransientParam(BasalforcingsUpperwaterMeltingRateEnum,&transparam[0],&transparam[M],interp,cycle,M));
			}
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&N,&M,"md.basalforcings.deepwater_elevation");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(BasalforcingsDeepwaterElevationEnum,transparam[0]));
			}
			else{
				_assert_(N==2);
				parameters->AddObject(new TransientParam(BasalforcingsDeepwaterElevationEnum,&transparam[0],&transparam[M],interp,cycle,M));
			}
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&N,&M,"md.basalforcings.upperwater_elevation");
			if(N==1){
				_assert_(M==1);
				parameters->AddObject(new DoubleParam(BasalforcingsUpperwaterElevationEnum,transparam[0]));
			}
			else{
				_assert_(N==2);
				parameters->AddObject(new TransientParam(BasalforcingsUpperwaterElevationEnum,&transparam[0],&transparam[M],interp,cycle,M));
			}
			xDelete<IssmDouble>(transparam);
			break;
		case SpatialLinearFloatingMeltRateEnum:
			/*Nothing to add to parameters:*/
			break;
		case MismipFloatingMeltRateEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.threshold_thickness",BasalforcingsThresholdThicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.upperdepth_melt",BasalforcingsUpperdepthMeltEnum));
			break;
		case MantlePlumeGeothermalFluxEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.mantleconductivity",BasalforcingsMantleconductivityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.nusselt",BasalforcingsNusseltEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.dtbg",BasalforcingsDtbgEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumeradius",BasalforcingsPlumeradiusEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.topplumedepth",BasalforcingsTopplumedepthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.bottomplumedepth",BasalforcingsBottomplumedepthEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumex",BasalforcingsPlumexEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.plumey",BasalforcingsPlumeyEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.crustthickness",BasalforcingsCrustthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.uppercrustthickness",BasalforcingsUppercrustthicknessEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.uppercrustheat",BasalforcingsUppercrustheatEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.lowercrustheat",BasalforcingsLowercrustheatEnum));
			break;
		case BasalforcingsPicoEnum:
			iomodel->FindConstant(&interp,"md.timestepping.interp_forcing");
			iomodel->FindConstant(&cycle,"md.timestepping.cycle_forcing");
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.num_basins",BasalforcingsPicoNumBasinsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.maxboxcount",BasalforcingsPicoMaxboxcountEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.gamma_T",BasalforcingsPicoGammaTEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.isplume",BasalforcingsPicoIsplumeEnum));
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.farocean_temperature");
			_assert_(M>=1 && N>=1);
			parameters->AddObject(new TransientArrayParam(BasalforcingsPicoFarOceantemperatureEnum,transparam,&transparam[N*(M-1)],interp,cycle,N,M));
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.farocean_salinity");
			_assert_(M>=1 && N>=1);
			parameters->AddObject(new TransientArrayParam(BasalforcingsPicoFarOceansalinityEnum,transparam,&transparam[N*(M-1)],interp,cycle,N,M));
			xDelete<IssmDouble>(transparam);
			break;
		case BasalforcingsIsmip6Enum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.num_basins",BasalforcingsIsmip6NumBasinsEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.gamma_0",BasalforcingsIsmip6Gamma0Enum));
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.islocal",BasalforcingsIsmip6IsLocalEnum));
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.delta_t");
			parameters->AddObject(new DoubleVecParam(BasalforcingsIsmip6DeltaTEnum,transparam,N));
			xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.tf_depths");
			parameters->AddObject(new DoubleVecParam(BasalforcingsIsmip6TfDepthsEnum,transparam,N));
			xDelete<IssmDouble>(transparam);
			break;
		case BeckmannGoosseFloatingMeltRateEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.isthermalforcing",BasalforcingsIsThermalForcingEnum));
			break;
		case LinearFloatingMeltRatearmaEnum:
			/*Add parameters that are not in standard nbvertices format*/
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.num_basins",BasalforcingsLinearNumBasinsEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.num_breaks",BasalforcingsLinearNumBreaksEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.num_params",BasalforcingsLinearNumParamsEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.ar_order",BasalforcingsARMAarOrderEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.ma_order",BasalforcingsARMAmaOrderEnum));
         parameters->AddObject(iomodel->CopyConstantObject("md.basalforcings.arma_timestep",BasalforcingsARMATimestepEnum));
         iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.datebreaks");
         parameters->AddObject(new DoubleMatParam(BasalforcingsARMAdatebreaksEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.polynomialparams");
         parameters->AddObject(new DoubleMatParam(BasalforcingsARMApolyparamsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.arlag_coefs");
         parameters->AddObject(new DoubleMatParam(BasalforcingsARMAarlagcoefsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.malag_coefs");
         parameters->AddObject(new DoubleMatParam(BasalforcingsARMAmalagcoefsEnum,transparam,M,N));
         xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.upperwater_melting_rate");
         parameters->AddObject(new DoubleVecParam(BasalforcingsUpperwaterMeltingRateEnum,transparam,N));
         xDelete<IssmDouble>(transparam);
         iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.upperwater_elevation");
         parameters->AddObject(new DoubleVecParam(BasalforcingsUpperwaterElevationEnum,transparam,N));
         xDelete<IssmDouble>(transparam);
			iomodel->FetchData(&transparam,&M,&N,"md.basalforcings.deepwater_elevation");
         parameters->AddObject(new DoubleVecParam(BasalforcingsDeepwaterElevationEnum,transparam,N));
         xDelete<IssmDouble>(transparam);
			break;
		default:
			_error_("Basal forcing model "<<EnumToStringx(basalforcing_model)<<" not supported yet");
	}

	/*some parameters that did not come with the iomodel: */
	parameters->AddObject(new IntParam(SolutionTypeEnum,solution_type));

	/*Time stepping*/
	parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.type",TimesteppingTypeEnum));
	iomodel->FindConstant(&timestepping_type,"md.timestepping.type");
	switch(timestepping_type){
		case FixedTimesteppingEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.start_time",TimesteppingStartTimeEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.final_time",TimesteppingFinalTimeEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.time_step",TimesteppingTimeStepEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.interp_forcing",TimesteppingInterpForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.average_forcing",TimesteppingAverageForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.cycle_forcing",TimesteppingCycleForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.coupling_time",TimesteppingCouplingTimeEnum));
			break;
		case AdaptiveTimesteppingEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.start_time",TimesteppingStartTimeEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.final_time",TimesteppingFinalTimeEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.time_step_min",TimesteppingTimeStepMinEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.time_step_max",TimesteppingTimeStepMaxEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.cfl_coefficient",TimesteppingCflCoefficientEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.interp_forcing",TimesteppingInterpForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.average_forcing",TimesteppingAverageForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.cycle_forcing",TimesteppingCycleForcingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.timestepping.coupling_time",TimesteppingCouplingTimeEnum));
			break;
		default:
			_error_("Time stepping \""<<EnumToStringx(timestepping_type)<<"\" not supported yet");
	}
	iomodel->FindConstant(&time,"md.timestepping.start_time");
	parameters->AddObject(new DoubleParam(TimeEnum,time));
	parameters->AddObject(new IntParam(StepEnum,0));

	/*By default, save all results*/
	parameters->AddObject(new BoolParam(SaveResultsEnum,true));

	/*Should we output results on nodes?*/
	iomodel->FindConstant(&outputonnodes,&numoutputs,"md.settings.results_on_nodes");
	parameters->AddObject(new IntParam(SettingsNumResultsOnNodesEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SettingsResultsOnNodesEnum,outputonnodes,numoutputs));
	iomodel->DeleteData(&outputonnodes,numoutputs,"md.settings.results_on_nodes");

	/*Requested outputs */
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.transient.requested_outputs");
	parameters->AddObject(new IntParam(TransientNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(TransientRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.transient.requested_outputs");

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.steadystate.requested_outputs");
	parameters->AddObject(new IntParam(SteadystateNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SteadystateRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.steadystate.requested_outputs");

	int materialstype;
	iomodel->FindConstant(&materialstype,"md.materials.type");

	switch(materialstype){
		case MaticeEnum:
		case MatdamageiceEnum:
		case MatenhancediceEnum:
		case MatestarEnum:
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_ice",MaterialsRhoIceEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_water",MaterialsRhoSeawaterEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_freshwater",MaterialsRhoFreshwaterEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.mu_water",MaterialsMuWaterEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.heatcapacity",MaterialsHeatcapacityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.thermalconductivity",MaterialsThermalconductivityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.temperateiceconductivity",MaterialsTemperateiceconductivityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.effectiveconductivity_averaging",MaterialsEffectiveconductivityAveragingEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.latentheat",MaterialsLatentheatEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.beta",MaterialsBetaEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.meltingpoint",MaterialsMeltingpointEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.constants.referencetemperature",ConstantsReferencetemperatureEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.mixed_layer_capacity",MaterialsMixedLayerCapacityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.thermal_exchange_velocity",MaterialsThermalExchangeVelocityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.constants.g",ConstantsGEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.constants.gravitational_constant",ConstantsNewtonGravityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.rheology_law",MaterialsRheologyLawEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.earth_density",MaterialsEarthDensityEnum));

			break;
		case MaterialsEnum:{
			int nnat,dummy;
			int* nature=NULL;
			iomodel->FetchData(&nature,&nnat,&dummy,"md.materials.nature");
			for(int i=0;i<nnat;i++){
				switch(IoCodeToEnumNature(nature[i])){
					case MatlithoEnum:
						break;
					case MaticeEnum:
					case MatdamageiceEnum:
					case MatenhancediceEnum:
					case MatestarEnum:
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_ice",MaterialsRhoIceEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_water",MaterialsRhoSeawaterEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_freshwater",MaterialsRhoFreshwaterEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.mu_water",MaterialsMuWaterEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.heatcapacity",MaterialsHeatcapacityEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.thermalconductivity",MaterialsThermalconductivityEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.temperateiceconductivity",MaterialsTemperateiceconductivityEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.effectiveconductivity_averaging",MaterialsEffectiveconductivityAveragingEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.latentheat",MaterialsLatentheatEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.beta",MaterialsBetaEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.meltingpoint",MaterialsMeltingpointEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.constants.referencetemperature",ConstantsReferencetemperatureEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.mixed_layer_capacity",MaterialsMixedLayerCapacityEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.thermal_exchange_velocity",MaterialsThermalExchangeVelocityEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.constants.g",ConstantsGEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rheology_law",MaterialsRheologyLawEnum));
						/*slc:*/
						break;
					case MathydroEnum:
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_ice",MaterialsRhoIceEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_water",MaterialsRhoSeawaterEnum));
						parameters->AddObject(iomodel->CopyConstantObject("md.materials.rho_freshwater",MaterialsRhoFreshwaterEnum));
						break;
				}
			}
			parameters->AddObject(iomodel->CopyConstantObject("md.materials.earth_density",MaterialsEarthDensityEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.constants.gravitational_constant",ConstantsNewtonGravityEnum));
			/*Free rssources:*/
			xDelete<int>(nature);
			break;

	}
		default:
			_error_("Material "<< EnumToStringx(materialstype) <<" not supported yet");
	}

	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	parameters->AddObject(new BoolParam(HydrologyIsWaterPressureArmaEnum,false));
	if(hydrology_model==HydrologydcEnum){
		IssmDouble sedcomp, sedporo, watcomp, rhofresh, g;
		iomodel->FindConstant(&sedcomp,"md.hydrology.sediment_compressibility");
		iomodel->FindConstant(&sedporo,"md.hydrology.sediment_porosity");
		iomodel->FindConstant(&watcomp,"md.hydrology.water_compressibility");
		iomodel->FindConstant(&rhofresh,"md.materials.rho_freshwater");
		iomodel->FindConstant(&g,"md.constants.g");

		parameters->AddObject(new DoubleParam(HydrologydcSedimentLayerCompressibilityEnum,(watcomp + sedcomp/sedporo)));
		parameters->AddObject(new DoubleParam(HydrologydcSedimentPoreWaterMassEnum,(rhofresh*g*sedporo)));

		bool isefficientlayer;
		iomodel->FindConstant(&isefficientlayer,"md.hydrology.isefficientlayer");
		if(isefficientlayer){
			IssmDouble eplcomp, eplporo;
			iomodel->FindConstant(&eplcomp,"md.hydrology.epl_compressibility");
			iomodel->FindConstant(&eplporo,"md.hydrology.epl_porosity");
			parameters->AddObject(new DoubleParam(HydrologydcEplLayerCompressibilityEnum,(watcomp + eplcomp/eplporo)));
			parameters->AddObject(new DoubleParam(HydrologydcEplPoreWaterMassEnum,(rhofresh*g*eplporo)));

		}
	}
	else if(hydrology_model==HydrologyshreveEnum){
		/*Nothing to add*/
	}
	else if(hydrology_model==HydrologyshaktiEnum){
		/*Nothing to add*/
	}
	else if(hydrology_model==HydrologypismEnum){
		/*Nothing to add*/
	}
	else if(hydrology_model==HydrologyGlaDSEnum){
		/*Nothing to add*/
	}
	else if(hydrology_model==HydrologyTwsEnum){
		/*Nothing to add*/
	}
	else if(hydrology_model==HydrologyarmapwEnum){
		parameters->SetParam(true,HydrologyIsWaterPressureArmaEnum);
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.num_basins",HydrologyNumBasinsEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.num_breaks",HydrologyarmaNumBreaksEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.num_params",HydrologyarmaNumParamsEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ar_order",HydrologyarmaarOrderEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ma_order",HydrologyarmamaOrderEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.arma_timestep",HydrologyarmaTimestepEnum));
      iomodel->FetchData(&transparam,&M,&N,"md.hydrology.datebreaks");
      parameters->AddObject(new DoubleMatParam(HydrologyarmadatebreaksEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.hydrology.polynomialparams");
      parameters->AddObject(new DoubleMatParam(HydrologyarmapolyparamsEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.hydrology.arlag_coefs");
      parameters->AddObject(new DoubleMatParam(HydrologyarmaarlagcoefsEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.hydrology.malag_coefs");
      parameters->AddObject(new DoubleMatParam(HydrologyarmamalagcoefsEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.hydrology.monthlyfactors");
      parameters->AddObject(new DoubleMatParam(HydrologyarmaMonthlyFactorsEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
   }
	else{
		_error_("Hydrology model "<<EnumToStringx(hydrology_model)<<" not supported yet");
	}

	if(materialstype==MatdamageiceEnum){
		iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.damage.requested_outputs");
		parameters->AddObject(new IntParam(DamageEvolutionNumRequestedOutputsEnum,numoutputs));
		if(numoutputs)parameters->AddObject(new StringArrayParam(DamageEvolutionRequestedOutputsEnum,requestedoutputs,numoutputs));
		iomodel->DeleteData(&requestedoutputs,numoutputs,"md.damage.requested_outputs");
	}

	bool isstochasticforcing;
   parameters->FindParam(&isstochasticforcing,StochasticForcingIsStochasticForcingEnum);
	/*Stochastic Effective Pressure false by default*/
	parameters->AddObject(new BoolParam(StochasticForcingIsWaterPressureEnum,false));
   if(isstochasticforcing){
      int num_fields,num_tcov,stochastic_dim;
      char** fields;
      parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.num_fields",StochasticForcingNumFieldsEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.defaultdimension",StochasticForcingDefaultDimensionEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.stochastictimestep",StochasticForcingTimestepEnum));
      parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.num_timescovariance",StochasticForcingNumTimesCovarianceEnum));
      iomodel->FindConstant(&fields,&num_fields,"md.stochasticforcing.fields");
      if(num_fields<1) _error_("no stochasticforcing fields found");
      int* stochasticforcing_enums = xNew<int>(num_fields);
      for(int i=0;i<num_fields;i++){
         stochasticforcing_enums[i] = StringToEnumx(fields[i]);
         xDelete<char>(fields[i]);
      }
      xDelete<char*>(fields);
      parameters->AddObject(new IntVecParam(StochasticForcingFieldsEnum,stochasticforcing_enums,num_fields));
      xDelete<int>(stochasticforcing_enums);
      parameters->AddObject(iomodel->CopyConstantObject("md.stochasticforcing.randomflag",StochasticForcingRandomflagEnum));
      iomodel->FetchData(&transparam,&M,&N,"md.stochasticforcing.dimensions");
      parameters->AddObject(new IntVecParam(StochasticForcingDimensionsEnum,transparam,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.stochasticforcing.timecovariance");
      parameters->AddObject(new DoubleVecParam(StochasticForcingTimeCovarianceEnum,transparam,N));
      xDelete<IssmDouble>(transparam);
      iomodel->FetchData(&transparam,&M,&N,"md.stochasticforcing.covariance");
      parameters->AddObject(new DoubleMatParam(StochasticForcingCovarianceEnum,transparam,M,N));
      xDelete<IssmDouble>(transparam);
   }

	/*Deal with mass flux segments: {{{*/
	iomodel->FetchData(&qmu_mass_flux_present,"md.qmu.mass_flux_segments_present");
	iomodel->FetchData(&autodiff_mass_flux_present,"md.autodiff.mass_flux_segments_present");

	if(qmu_mass_flux_present || autodiff_mass_flux_present)mass_flux_present=true;
	else mass_flux_present=false;
	parameters->AddObject(new BoolParam(MassFluxSegmentsPresentEnum,mass_flux_present));

	if(mass_flux_present){

		/*Fetch the mass flux segments necessary to compute the mass fluxes.  Build a DoubleMatArrayParam object out of them: */
		iomodel->FetchData(&array,&mdims_array,&ndims_array,&mass_flux_num_profiles,"md.qmu.mass_flux_segments");
		if(mass_flux_num_profiles==0)_error_("mass_flux_num_profiles is 0, when MassFlux computations were requested!");

		/*Go through segments, and extract those that belong to this cpu: */
		for(i=0;i<mass_flux_num_profiles;i++){
			temp_matrix=array[i];
			temp_m=mdims_array[i];
			temp_n=ndims_array[i];
			_assert_(temp_n==5);

			m=0;
			for(j=0;j<temp_m;j++){
				if (  iomodel->my_elements[reCast<int>(*(temp_matrix+5*j+4))-1] )m++;
			}
			if(m){
				matrix=xNewZeroInit<IssmDouble>(5*m);
				count=0;
				for(j=0;j<temp_m;j++){
					if (iomodel->my_elements[reCast<int>(*(temp_matrix+5*j+4))-1]){
						for(k=0;k<5;k++)*(matrix+5*count+k)=*(temp_matrix+5*j+k);
						count++;
					}
				}
			}
			else{
				matrix=NULL;
			}

			/*Assign: */
			array[i]=matrix;
			mdims_array[i]=m;
			ndims_array[i]=5;

			/*Free temporary matrix: */
			xDelete<IssmDouble>(temp_matrix);
		}

		/*Ok, we have an array of segments, different on every cpu. Create a DoubleMatArrayParam object with it: */
		parameters->AddObject(new DoubleMatArrayParam(MassFluxSegmentsEnum,array,mass_flux_num_profiles,mdims_array,ndims_array));

		/*Free data: */
		for(i=0;i<mass_flux_num_profiles;i++){
			IssmDouble* matrix=array[i];
			xDelete<IssmDouble>(matrix);
		}
		xDelete<int>(mdims_array);
		xDelete<int>(ndims_array);
		xDelete<IssmDouble*>(array);
	}
	/*}}}*/

	/*Before returning, create parameters in case we are running Qmu or control types runs: */
	CreateParametersControl(parameters,iomodel,solution_type);

	#ifdef _HAVE_DAKOTA_
	CreateParametersDakota(parameters,iomodel,rootpath);
	#endif

	/*Now, deal with toolkits options, which need to be put into the parameters dataset: */
	ParseToolkitsOptionsx(parameters,toolkitsoptionsfid);

	#ifdef _HAVE_AD_
	if(VerboseMProcessor()) _printf0_("   starting autodiff parameters \n");
	CreateParametersAutodiff(parameters,iomodel);
	if(VerboseMProcessor()) _printf0_("   ending autodiff parameters \n");
	#endif
}
