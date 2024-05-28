/*!\file: levelset_core.cpp
 * \brief: levelset-module to update the ice domain
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../solutionsequences/solutionsequences.h"
#include "../modules/modules.h"

void movingfront_core(FemModel* femmodel){

	/*Start profiler*/
	femmodel->profiler->Start(MOVINGFRONTCORE);

	/* intermediaries */
	bool save_results,isstressbalance,ismasstransport,isthermal,isenthalpy,islevelset,ismovingfront,killicebergs;
	int  domaintype, reinit_frequency,step;
	Analysis  *analysis=NULL;
	IssmDouble maxVel;

	/* recover parameters */
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isstressbalance,TransientIsstressbalanceEnum);
	femmodel->parameters->FindParam(&ismasstransport,TransientIsmasstransportEnum);
	femmodel->parameters->FindParam(&isthermal,TransientIsthermalEnum);
	femmodel->parameters->FindParam(&ismovingfront,TransientIsmovingfrontEnum);
	femmodel->parameters->FindParam(&reinit_frequency,LevelsetReinitFrequencyEnum);
	femmodel->parameters->FindParam(&killicebergs,LevelsetKillIcebergsEnum);
	femmodel->parameters->FindParam(&step,StepEnum);
	if(isthermal && domaintype==Domain3DEnum) femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);

	if(!ismovingfront) return;

	/* Many calving parameterizations and the level set equations require depth
	 * average velocities so do this calculation once for all here */
	if(domaintype!=Domain2DhorizontalEnum){
		femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
		if(domaintype==Domain3DEnum){
			femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
		}
	}

	/* smoothen slope of lsf for computation of normal on ice domain*/
	levelsetfunctionslope_core(femmodel);

	/* compute the maximal velocity over the whole domain */
	if(isstressbalance){
		femmodel->MaxVelx(&maxVel);
		femmodel->parameters->SetParam(maxVel, CalvingVelMaxEnum);
	}

	/* start the work from here */
	if(VerboseSolution()) _printf0_("   computing calving and undercutting\n");
	Calvingx(femmodel);
	FrontalForcingsx(femmodel);
	if(VerboseSolution()) _printf0_("   computing new ice front position\n");

	/* determine variables for extrapolation */
	std::vector<int>  extrapol_vars;
	if(isstressbalance){
		extrapol_vars.push_back(VxEnum);
		extrapol_vars.push_back(VyEnum);
		if(domaintype==Domain3DEnum) extrapol_vars.push_back(VzEnum);
	}
	if(ismasstransport) extrapol_vars.push_back(ThicknessEnum);
	if(isthermal && domaintype==Domain3DEnum){
		if(isenthalpy) extrapol_vars.push_back(EnthalpyEnum);
		else           extrapol_vars.push_back(TemperatureEnum);
	}

	/* extrapolate */
	analysis = new ExtrapolationAnalysis();
	for(int iv=0;iv<extrapol_vars.size();iv++){
		femmodel->parameters->SetParam(extrapol_vars[iv],ExtrapolationVariableEnum); 
		analysis->Core(femmodel);
	}
	delete analysis;	

	/* Need to do it again after extrapolation! */
	if(domaintype!=Domain2DhorizontalEnum){
		femmodel->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
		femmodel->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(femmodel);
		if(domaintype==Domain3DEnum){
			femmodel->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
			femmodel->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
			depthaverage_core(femmodel);
		}
	}

	/* Calculate the frontal velocity for levelset function */
	MovingFrontalVelx(femmodel);

	/* solve level set equation */
	LevelsetAnalysis lsanalysis;
	lsanalysis.Core(femmodel);
	lsanalysis.PostProcess(femmodel);

	/*Kill ice berg to avoid free body motion*/
	if(killicebergs){
		int killberg = 0;
		if(VerboseSolution()) _printf0_("   looking for icebergs to kill\n");
		killberg = KillIcebergsx(femmodel);
		/*wait for all cores*/
		int totalkill;
		ISSM_MPI_Reduce(&killberg,&totalkill,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(&totalkill,1,ISSM_MPI_INT,0,IssmComm::GetComm());

		if (totalkill > 0) {
			if(VerboseSolution()) _printf0_("   reinitializing level set after killing " << totalkill << " icebergs\n");
			femmodel->ResetLevelset();
			ResetBoundaryConditions(femmodel,LevelsetAnalysisEnum);
		}
	}

	/*Reset levelset if needed*/
	if(reinit_frequency && (step%reinit_frequency==0)){
		if(VerboseSolution()) _printf0_("   reinitializing level set\n");
		femmodel->ResetLevelset();
		ResetBoundaryConditions(femmodel,LevelsetAnalysisEnum);
	}

	/* update vertices included for next calculation */
	GetMaskOfIceVerticesLSMx(femmodel);

	/*Save results*/
	if(save_results){
		int outputs[1] = {MaskIceLevelsetEnum};
		femmodel->RequestedOutputsx(&femmodel->results,&outputs[0],1);
	}

	/*End profiler*/
	femmodel->profiler->Stop(MOVINGFRONTCORE);
}
