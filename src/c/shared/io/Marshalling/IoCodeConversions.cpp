#include <math.h>
#include <float.h>
#include <cstring>

#include "../../Enum/EnumDefinitions.h"
#include "../../MemOps/MemOps.h"
#include "../../Exceptions/exceptions.h"

void FieldAndEnumFromCode(int* out_enum,char** pfield,const char* string_in){/*{{{*/

	/*output*/
	char* fieldname = NULL;
	int   input_enum = -1;

	if(strcmp(string_in,"Thickness")==0 || strcmp(string_in,"md.geometry.thickness")==0){
		const char* field = "md.geometry.thickness";
		input_enum        = ThicknessEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"MaterialsRheologyBbar")==0){
		const char* field = "md.materials.rheology_B";
		input_enum        = MaterialsRheologyBbarEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"MaterialsRheologyB")==0){
		const char* field = "md.materials.rheology_B";
		input_enum        = MaterialsRheologyBEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"MaterialsRheologyN")==0){
		const char* field = "md.materials.rheology_n";
		input_enum        = MaterialsRheologyNEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbMassBalance")==0){
		const char* field = "md.smb.mass_balance";
		input_enum        = SmbMassBalanceEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbAccumulation")==0){
		const char* field = "md.smb.accumulation";
		input_enum        = SmbAccumulationEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbMelt")==0){
		const char* field = "md.smb.melt";
		input_enum        = SmbMeltEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbRefreeze")==0){
		const char* field = "md.smb.refreeze";
		input_enum        = SmbRefreezeEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbRunoff")==0){
		const char* field = "md.smb.runoff";
		input_enum        = SmbRunoffEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbEvaporation")==0){
		const char* field = "md.smb.evaporation";
		input_enum        = SmbEvaporationEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbTa")==0){
		const char* field = "md.smb.Ta";
		input_enum        = SmbTaEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbV")==0){
		const char* field = "md.smb.V";
		input_enum        = SmbVEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbDswrf")==0){
		const char* field = "md.smb.dswrf";
		input_enum        = SmbDswrfEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbDlwrf")==0){
		const char* field = "md.smb.dlwrf";
		input_enum        = SmbDlwrfEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbP")==0){
		const char* field = "md.smb.P";
		input_enum        = SmbPEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbEAir")==0){
		const char* field = "md.smb.eAir";
		input_enum        = SmbEAirEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbPAir")==0){
		const char* field = "md.smb.pAir";
		input_enum        = SmbPAirEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbVz")==0){
		const char* field = "md.smb.Vz";
		input_enum        = SmbVzEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbTz")==0){
		const char* field = "md.smb.Tz";
		input_enum        = SmbTzEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"SmbC")==0){
		const char* field = "md.smb.C";
		input_enum        = SmbCEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsFloatingiceMeltingRate")==0){
		const char* field = "md.basalforcings.floatingice_melting_rate";
		input_enum        = BasalforcingsFloatingiceMeltingRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsGeothermalflux")==0){
		const char* field = "md.basalforcings.geothermalflux";
		input_enum        = BasalforcingsGeothermalfluxEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"FrictionCoefficient")==0 || strcmp(string_in,"md.friction.coefficient")==0){
		const char* field = "md.friction.coefficient";
		input_enum        = FrictionCoefficientEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"FrictionC")==0 || strcmp(string_in,"md.friction.C")==0){
		const char* field = "md.friction.C";
		input_enum        = FrictionCEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"FrictionEffectivePressure")==0){
		const char* field = "md.friction.effective_pressure";
		input_enum        = FrictionEffectivePressureEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"Vx")==0){
		 const char* field = "md.initialization.vx";
		 input_enum        = VxEnum;
		 fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	 }
	else if(strcmp(string_in,"Vy")==0){
		 const char* field = "md.initialization.vy";
		 input_enum        = VyEnum;
		 fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	 }
	else if(strcmp(string_in,"BalancethicknessThickeningRate")==0){
		 const char* field = "md.balancethickness.thickening_rate";
		 input_enum        = BalancethicknessThickeningRateEnum;
		 fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BalancethicknessSpcthickness")==0){
		 const char* field = "md.balancethickness.spcthickness";
		 input_enum        = BalancethicknessSpcthicknessEnum;
		 fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"CalvingMeltingrate")==0){
		const char* field = "md.calving.meltingrate";
		input_enum        = CalvingMeltingrateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"CalvingStressThresholdGroundedice")==0){
		const char* field = "md.calving.stress_threshold_groundedice";
		input_enum        = CalvingStressThresholdGroundediceEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"CalvingADStressThresholdGroundedice")==0){
		const char* field = "md.calving.stress_threshold_groundedice";
		input_enum        = CalvingADStressThresholdGroundediceEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"DamageDbar")==0){
		const char* field = "md.damage.D";
		input_enum        = DamageDbarEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"IceLoad")==0){
		const char* field = "md.masstransport.spcthickness";
		input_enum        = MasstransportSpcthicknessEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"NGiaRate")==0){
		const char* field = "md.gia.Ngia";
		input_enum        = NGiaRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"DslGlobalAverageThermostericSeaLevel")==0){
		const char* field = "md.dsl.global_average_thermosteric_sea_level";
		input_enum        = OceantransportSpcstrEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"DslSeaWaterPressureAtSeaFloor")==0){
		const char* field = "md.dsl.sea_water_pressure_at_sea_floor";
		input_enum        = OceantransportSpcbottompressureEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"DslSeaSurfaceHeightAboveGeoid")==0){
		const char* field = "md.dsl.sea_surface_height_above_geoid";
		input_enum        = OceantransportSpcdslEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"TwsLoad")==0){
		const char* field = "md.hydrology.spcwatercolumn";
		input_enum        = HydrologyTwsSpcEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"UGiaRate")==0){
		const char* field = "md.gia.Ugia";
		input_enum        = UGiaRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"MaskIceLevelset")==0){
		const char* field = "md.mask.ice_levelset";
		input_enum        = MaskIceLevelsetEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsPerturbationMeltingRate")==0){
		const char* field = "md.basalforcings.perturbation_melting_rate";
		input_enum        = BasalforcingsPerturbationMeltingRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsMeltrateFactor")==0){
		const char* field = "md.basalforcings.meltrate_factor";
		input_enum        = BasalforcingsMeltrateFactorEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsSpatialDeepwaterMeltingRate")==0){
		const char* field = "md.basalforcings.deepwater_melting_rate";
		input_enum        = BasalforcingsSpatialDeepwaterMeltingRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"BasalforcingsDeepwaterMeltingRate")==0){
		const char* field = "md.basalforcings.deepwater_melting_rate";
		input_enum        = BasalforcingsDeepwaterMeltingRateEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else if(strcmp(string_in,"Bed")==0){
		const char* field = "md.geometry.bed";
		input_enum        = BedEnum;
		fieldname=xNew<char>((strlen(field)+1)); xMemCpy<char>(fieldname,field,(strlen(field)+1));
	}
	else{
		_error_("Field \""<<string_in<<"\" not supported yet");
	}

	/*Assign output pointers*/
	*out_enum = input_enum;
	*pfield   = fieldname;
	return;
}/*}}}*/
int IoCodeToEnumSMB(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return SMBforcingEnum;
		case 2: return SMBcomponentsEnum;
		case 3: return SMBmeltcomponentsEnum;
		case 4: return SMBpddEnum;
		case 5: return SMBd18opddEnum;
		case 6: return SMBgradientsEnum;
		case 7: return SMBhenningEnum;
		case 8: return SMBgembEnum;
		case 9: return SMBgradientselaEnum;
		case 10: return SMBpddSicopolisEnum;
		case 11: return SMBgradientscomponentsEnum;
		case 12: return SMBsemicEnum;	 
		case 13: return SMBarmaEnum;
		case 14: return SMBdebrisEvattEnum;
		case 15: return SMBpddGCMEnum;
		default: _error_("Marshalled SMB code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumBasal(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return FloatingMeltRateEnum;
		case 2: return LinearFloatingMeltRateEnum;
		case 3: return MismipFloatingMeltRateEnum;
		case 4: return MantlePlumeGeothermalFluxEnum;
		case 5: return BasalforcingsPicoEnum;
		case 6: return SpatialLinearFloatingMeltRateEnum;
		case 7: return BasalforcingsIsmip6Enum;
		case 8: return BeckmannGoosseFloatingMeltRateEnum;
		case 9: return LinearFloatingMeltRatearmaEnum;
		default: _error_("Marshalled Basal Forcings code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumCalving(int enum_in){/*{{{*/
	switch(enum_in){
		case 1:  return DefaultCalvingEnum;
		case 2:  return CalvingVonmisesEnum;
		case 3:  return CalvingLevermannEnum;
		case 4:  return CalvingMinthicknessEnum;
		case 5:  return CalvingHabEnum;
		case 6:  return CalvingCrevasseDepthEnum;
		case 7:  return CalvingDev2Enum;
		case 8:  return CalvingTestEnum;
		case 9:  return CalvingParameterizationEnum;
		case 10: return CalvingPollardEnum;
		case 11: return CalvingVonmisesADEnum;
		case 12:  return CalvingCalvingMIPEnum;
		default: _error_("Marshalled Calving law code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumFrontalforcings(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return FrontalForcingsDefaultEnum;
		case 2: return FrontalForcingsRignotEnum;
		case 3: return FrontalForcingsRignotarmaEnum;
		default: _error_("Marshalled Frontalforcings code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumHydrology(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return HydrologydcEnum;
		case 2: return HydrologyshreveEnum;
		case 3: return HydrologyshaktiEnum;
		case 4: return HydrologypismEnum;
		case 5: return HydrologyGlaDSEnum;
		case 6: return HydrologyTwsEnum;
		case 7: return HydrologyarmapwEnum;
		default: _error_("Marshalled hydrology code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumMaterials(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return MatdamageiceEnum;
		case 2: return MatestarEnum;
		case 3: return MaticeEnum;
		case 4: return MatenhancediceEnum;
		case 5: return MaterialsEnum; //This should not happen anymore??
		case 6: return MatlithoEnum;
		case 7: return MathydroEnum;
		default: _error_("Marshalled materials code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumNature(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return MatdamageiceEnum;
		case 2: return MatestarEnum;
		case 3: return MaticeEnum;
		case 4: return MatenhancediceEnum;
		case 5: return MaterialsEnum;
		case 6: return MatlithoEnum;
		case 7: return MathydroEnum;
		default: _error_("Marshalled materials nature code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumTimestepping(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return FixedTimesteppingEnum;
		case 2: return AdaptiveTimesteppingEnum;
		default: _error_("Marshalled materials code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumAmr(int enum_in){/*{{{*/
	switch(enum_in){
		case 1: return AmrBamgEnum;
		case 2: return AmrNeopzEnum;
		default: _error_("Marshalled AMR code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/
int IoCodeToEnumGrd(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return NoneEnum;
		case 1: return ElasticEnum;
		case 2: return IvinsEnum;
		default: _error_("Marshalled GRD code \""<<enum_in<<"\" not supported yet");
	}
}/*}}}*/

int IoCodeToEnumVertexEquation(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return NoneApproximationEnum;
		case 1: return SIAApproximationEnum;
		case 2: return SSAApproximationEnum;
		case 3: return L1L2ApproximationEnum;
		case 4: return MOLHOApproximationEnum;
		case 5: return HOApproximationEnum;
		case 6: return FSApproximationEnum;
		case 7: return SSAHOApproximationEnum;
		case 8: return HOFSApproximationEnum;
		case 9: return SSAFSApproximationEnum;
		default: _error_("Marshalled vertex equation code \""<<enum_in<<"\" not supported yet.");
	}
}/*}}}*/
int IoCodeToEnumElementEquation(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return NoneApproximationEnum;
		case 1: return SIAApproximationEnum;
		case 2: return SSAApproximationEnum;
		case 3: return L1L2ApproximationEnum;
		case 4: return MOLHOApproximationEnum;
		case 5: return HOApproximationEnum;
		case 6: return FSApproximationEnum;
		case 7: return SSAHOApproximationEnum;
		case 8: return SSAFSApproximationEnum;
		case 9: return HOFSApproximationEnum;
		default: _error_("Marshalled element equation code \""<<enum_in<<"\" not supported yet.");
	}

}/*}}}*/

int IoRiftfillToEnum(int enum_in){/*{{{*/
	switch(enum_in){
		case 0: return AirEnum;
		case 1: return IceEnum;
		case 2: return MelangeEnum;
		case 3: return WaterEnum;
		default: _error_("Marshalled Riftfill enum \""<<enum_in<<"\" not supported yet.");
	}
}/*}}}*/
