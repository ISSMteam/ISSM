/*!\file: transient_3d_core.cpp
 * \brief: core of the transient_3d solution
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <float.h>
#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"
#include "../toolkits/codipack/CoDiPackGlobal.h"

/*Prototypes*/
void transient_step(FemModel* femmodel);

void transient_core(FemModel* femmodel){/*{{{*/

	/*parameters: */
	IssmDouble finaltime,dt,yts;
	bool       iscontrol,isautodiff;
	int        timestepping;
	int        output_frequency,checkpoint_frequency;
	int        amr_frequency;
	char     **requested_outputs = NULL;

	/*intermediary: */
	int        step;
	IssmDouble time;

	/*first, figure out if there was a check point, if so, do a reset of the FemModel* femmodel structure. */
	femmodel->parameters->FindParam(&checkpoint_frequency,SettingsCheckpointFrequencyEnum);
	if(checkpoint_frequency) femmodel->Restart();

	/*then recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
	femmodel->parameters->FindParam(&output_frequency,SettingsOutputFrequencyEnum);
	femmodel->parameters->FindParam(&timestepping,TimesteppingTypeEnum);
	femmodel->parameters->FindParam(&amr_frequency,TransientAmrFrequencyEnum);
	femmodel->parameters->FindParam(&iscontrol,InversionIscontrolEnum);
	femmodel->parameters->FindParam(&isautodiff,AutodiffIsautodiffEnum);

	/*call modules that are not dependent on time stepping:*/
	transient_precore(femmodel);

	while(time < finaltime - (yts*DBL_EPSILON)){ //make sure we run up to finaltime.

		/*Time Increment*/
		switch(timestepping){
			case AdaptiveTimesteppingEnum:
				femmodel->TimeAdaptx(&dt);

				break;
			case FixedTimesteppingEnum:
				femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
				break;
			default:
				_error_("Time stepping \""<<EnumToStringx(timestepping)<<"\" not supported yet");
		}

		/*Do not exceed final time*/
		if(time+dt>finaltime) dt=finaltime-time;
		femmodel->parameters->SetParam(dt,TimesteppingTimeStepEnum);

		/*Set new step number and time step size*/
		step+=1;
		time+=dt;
		femmodel->parameters->SetParam(time,TimeEnum);
		femmodel->parameters->SetParam(step,StepEnum);

		if(VerboseSolution()){
			//_printf0_("iteration " << step << "/" << ceil((finaltime-time)/dt)+step << \
			//			"  time [yr]: " <<std::fixed<<setprecision(2)<< time/yts << " (time step: " << dt/yts << ")\n");
			_printf0_("\e[92miteration " << step << "/" << ceil((finaltime-time)/dt)+step << \
						"  time [yr]: " <<std::fixed<<setprecision(2)<< time/yts << "\e[m (time step: " << dt/yts << ")\n");
		}
		bool save_results=false;
		if(step%output_frequency==0 || (time >= finaltime - (yts*DBL_EPSILON)) || step==1) save_results=true;
		femmodel->parameters->SetParam(save_results,SaveResultsEnum);

		/*Run transient step!*/
		transient_step(femmodel);

		/*unload results*/
		if(save_results){
			if(VerboseSolution()) _printf0_("   saving temporary results\n");
			OutputResultsx(femmodel);
		}

		if(checkpoint_frequency && (step%checkpoint_frequency==0)){
			if(VerboseSolution()) _printf0_("   checkpointing model \n");
			femmodel->CheckPoint();
		}

		/*Adaptive mesh refinement*/
		if(amr_frequency){

#if !defined(_HAVE_AD_)
			if(save_results) femmodel->WriteMeshInResults();
			if(step%amr_frequency==0 && time<finaltime){
				if(VerboseSolution()) _printf0_("   refining mesh\n");
				femmodel->ReMesh();//Do not refine the last step
			}

#else
			_error_("AMR not suppored with AD");
#endif
		}

		if(iscontrol && isautodiff){
			/*Go through our dependent variables, and compute the response:*/
			DataSet* dependent_objects=((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
			for(Object* & object:dependent_objects->objects){
				DependentObject* dep=(DependentObject*)object;
				dep->RecordResponsex(femmodel);
			}
		}
	}

	if(!iscontrol || !isautodiff) femmodel->RequestedDependentsx();

	/*finalize:*/
	transient_postcore(femmodel);

}/*}}}*/
void transient_step(FemModel* femmodel){/*{{{*/

	/*parameters: */
	bool isstressbalance,ismasstransport,ismmemasstransport,isage,isoceantransport,issmb,isthermal,isgroundingline,isesa,issampling;
	bool isslc,ismovingfront,isdamageevolution,ishydrology,isstochasticforcing,save_results;
	bool isdebris;
	int  step,sb_coupling_frequency,isoceancoupling;
	int  domaintype,numoutputs;

	/*then recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&sb_coupling_frequency,SettingsSbCouplingFrequencyEnum);
	femmodel->parameters->FindParam(&isstressbalance,TransientIsstressbalanceEnum);
	femmodel->parameters->FindParam(&ismasstransport,TransientIsmasstransportEnum);
	femmodel->parameters->FindParam(&ismmemasstransport,TransientIsmmemasstransportEnum);
	femmodel->parameters->FindParam(&isage,TransientIsageEnum);
	femmodel->parameters->FindParam(&isoceantransport,TransientIsoceantransportEnum);
	femmodel->parameters->FindParam(&issmb,TransientIssmbEnum);
	femmodel->parameters->FindParam(&isthermal,TransientIsthermalEnum);
	femmodel->parameters->FindParam(&isesa,TransientIsesaEnum);
	femmodel->parameters->FindParam(&isslc,TransientIsslcEnum);
	femmodel->parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);
	femmodel->parameters->FindParam(&ismovingfront,TransientIsmovingfrontEnum);
	femmodel->parameters->FindParam(&isoceancoupling,TransientIsoceancouplingEnum);
	femmodel->parameters->FindParam(&isdamageevolution,TransientIsdamageevolutionEnum);
	femmodel->parameters->FindParam(&ishydrology,TransientIshydrologyEnum);
	femmodel->parameters->FindParam(&isdebris,TransientIsdebrisEnum);
	femmodel->parameters->FindParam(&issampling,TransientIssamplingEnum);
	femmodel->parameters->FindParam(&numoutputs,TransientNumRequestedOutputsEnum);
	femmodel->parameters->FindParam(&isstochasticforcing,StochasticForcingIsStochasticForcingEnum);

	if(isstochasticforcing) StochasticForcingx(femmodel);

	if(isthermal && domaintype==Domain3DEnum){
		if(issmb){
			bool isenthalpy;
			int  smb_model;
			femmodel->parameters->FindParam(&isenthalpy,ThermalIsenthalpyEnum);
			femmodel->parameters->FindParam(&smb_model,SmbEnum);
			if(isenthalpy){
				if(smb_model==SMBpddEnum || smb_model==SMBd18opddEnum || smb_model==SMBpddSicopolisEnum){
					femmodel->SetCurrentConfiguration(EnthalpyAnalysisEnum);
					ResetBoundaryConditions(femmodel,EnthalpyAnalysisEnum);
				}
			}
			else{
				if(smb_model==SMBpddEnum || smb_model==SMBd18opddEnum || smb_model==SMBpddSicopolisEnum){
					femmodel->SetCurrentConfiguration(ThermalAnalysisEnum);
					ResetBoundaryConditions(femmodel,ThermalAnalysisEnum);
				}
			}
		}
		thermal_core(femmodel);
	}

	/* Using Hydrology dc  coupled we need to compute smb in the hydrology inner time loop*/
	if(issmb) smb_core(femmodel);

	if(ishydrology) hydrology_core(femmodel);

	if(isstressbalance && (step%sb_coupling_frequency==0 || step==1)) stressbalance_core(femmodel);

	if(isdamageevolution) damage_core(femmodel);

	if(ismovingfront)	movingfront_core(femmodel);

	if(isdebris) debris_core(femmodel);

#if defined(_HAVE_OCEAN_)
	if(isoceancoupling) {
		/*First calculate thickness change without melt (dynamic thinning) to send to ocean
		 * then receive ocean melt 
		 * then go back to the previous geometry to continue the transient with the melt received*/
		InputUpdateFromConstantx(femmodel,0.,BasalforcingsFloatingiceMeltingRateEnum,P1Enum);
		masstransport_core(femmodel);
		OceanExchangeDatax(femmodel,false);
		InputDuplicatex(femmodel,ThicknessOldEnum,ThicknessEnum);
		InputDuplicatex(femmodel,BaseOldEnum,BaseEnum);
		InputDuplicatex(femmodel,SurfaceOldEnum,SurfaceEnum);
	}
#endif

	/* from here on, prepare geometry for next time step*/

	if(ismasstransport){
		bmb_core(femmodel);
		masstransport_core(femmodel);
	}
	if(ismmemasstransport){
		mmemasstransport_core(femmodel);
	}

	if(isoceantransport) oceantransport_core(femmodel);

	if(isgroundingline) groundingline_core(femmodel);

	/*Update mesh vertices now that we have changed the geometry*/
	if(ismasstransport || isgroundingline) femmodel->UpdateVertexPositionsx();

	if(isesa) esa_core(femmodel);

	/*Sea level change: */
	if(isslc){
#ifdef _HAVE_SEALEVELCHANGE_
		sealevelchange_core(femmodel);
#else
		_error_("Compiled with SeaLevelChange capability");
#endif
	}

	/*Sampling: */
	if(issampling) sampling_core(femmodel);

	/*Any requested output that needs to be saved?*/
	if(numoutputs){
		char **requested_outputs = NULL;
		femmodel->parameters->FindParam(&requested_outputs,&numoutputs,TransientRequestedOutputsEnum);

		if(VerboseSolution()) _printf0_("   computing transient requested outputs\n");
		femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs,save_results);

		/*Free resources:*/
		for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);
	}
}/*}}}*/
void transient_precore(FemModel* femmodel){/*{{{*/

	bool       isslc,isuq;
	int        amr_frequency,amr_restart,isoceancoupling;

	femmodel->parameters->FindParam(&isoceancoupling,TransientIsoceancouplingEnum);
	femmodel->parameters->FindParam(&amr_frequency,TransientAmrFrequencyEnum);
	femmodel->parameters->FindParam(&isslc,TransientIsslcEnum);
	femmodel->parameters->FindParam(&isuq,QmuIsdakotaEnum);

#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
	if(amr_frequency){
		femmodel->parameters->FindParam(&amr_restart,AmrRestartEnum);
		if(amr_restart) femmodel->ReMesh();
	}
#endif

#if defined(_HAVE_OCEAN_ )
	if(isoceancoupling) OceanExchangeDatax(femmodel,true);
#endif

#if defined(_HAVE_SEALEVELCHANGE_)
	if(isslc) sealevelchange_initialgeometry(femmodel);
#endif

	//Resolve Mmes prior to running in transient: 
	if(!isuq)UpdateMmesx(femmodel);

}/*}}}*/
void transient_postcore(FemModel* femmodel){/*{{{*/

	bool       isslc;
	femmodel->parameters->FindParam(&isslc,TransientIsslcEnum);

	#if defined(_HAVE_SEALEVELCHANGE_)
	if(isslc) sealevelchange_finalize(femmodel);
	#endif
	
}/*}}}*/

#ifdef _HAVE_CODIPACK_
double transient_ad(FemModel* femmodel, double* G, double* Jlist){/*{{{*/

	/*parameters: */
	IssmDouble finaltime,dt,yts,time;
	int        isoceancoupling;
	int        step,timestepping;
	int        checkpoint_frequency,num_responses;
	int		 *control_enum;

	/*Get rank*/
	int my_rank = IssmComm::GetRank();

	/*then recover parameters common to all solutions*/
	femmodel->parameters->FindParam(&step,StepEnum);
	femmodel->parameters->FindParam(&time,TimeEnum);
	femmodel->parameters->FindParam(&finaltime,TimesteppingFinalTimeEnum);
	femmodel->parameters->FindParam(&yts,ConstantsYtsEnum);
	femmodel->parameters->FindParam(&timestepping,TimesteppingTypeEnum);
	femmodel->parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	femmodel->parameters->FindParam(&checkpoint_frequency,SettingsCheckpointFrequencyEnum); _assert_(checkpoint_frequency>0);
	femmodel->parameters->FindParam(&control_enum,NULL,InversionControlParametersEnum);

	std::vector<IssmDouble> time_all;
	std::vector<IssmDouble> dt_all;
	std::vector<int>        checkpoint_steps;
	int                     Ysize = 0;
	CoDi_global            codi_y_data = {};
	CountDoublesFunctor   *hdl_countdoubles = NULL;
	RegisterInputFunctor  *hdl_regin        = NULL;
	RegisterOutputFunctor *hdl_regout       = NULL;

	while(time < finaltime - (yts*DBL_EPSILON)){ //make sure we run up to finaltime.

		/*Time Increment*/
		switch(timestepping){
			case AdaptiveTimesteppingEnum:
				femmodel->TimeAdaptx(&dt);
				if(time+dt>finaltime) dt=finaltime-time;
				femmodel->parameters->SetParam(dt,TimesteppingTimeStepEnum);
				break;
			case FixedTimesteppingEnum:
				femmodel->parameters->FindParam(&dt,TimesteppingTimeStepEnum);
				break;
			default:
				_error_("Time stepping \""<<EnumToStringx(timestepping)<<"\" not supported yet");
		}
		step+=1;
		time+=dt;
		femmodel->parameters->SetParam(time,TimeEnum);
		femmodel->parameters->SetParam(step,StepEnum);
		femmodel->parameters->SetParam(false,SaveResultsEnum);
		time_all.push_back(time);
		dt_all.push_back(dt);

		if(VerboseSolution()){
			_printf0_("iteration " << step << "/" << ceil((finaltime-time)/dt)+step << \
					"  time [yr]: " <<std::fixed<<setprecision(2)<< time/yts << " (time step: " << dt/yts << ")\n");
		}

		/*Store Model State at the beginning of the step*/
		if(step%checkpoint_frequency==0 || step==1){
			if(VerboseSolution()) _printf0_("   checkpointing model (step: "<<step<<")\n");
			femmodel->CheckPointAD(step);
			checkpoint_steps.push_back(step);
		}

		/*Run transient step!*/
		transient_step(femmodel);

		/*Go through our dependent variables, and compute the response:*/
		DataSet* dependent_objects=((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
		for(Object* & object:dependent_objects->objects){
			DependentObject* dep=(DependentObject*)object;
			dep->RecordResponsex(femmodel);
		}

		if(VerboseSolution()) _printf0_("   counting number of active variables\n");
		hdl_countdoubles = new CountDoublesFunctor();
		femmodel->Marshall(hdl_countdoubles);
		if(hdl_countdoubles->DoubleCount()>Ysize) Ysize= hdl_countdoubles->DoubleCount();
		delete hdl_countdoubles;
	}

	int finalstep = step;
	if(VerboseSolution()) _printf0_("   done with initial complete transient\n");

	/*__________________________________________________________________________________*/

	/*Get X (control)*/
	IssmDouble *X = NULL; int Xsize;
	GetVectorFromControlInputsx(&X,&Xsize,femmodel->elements,femmodel->nodes,femmodel->vertices,femmodel->loads,femmodel->materials,femmodel->parameters,"value");

	/*Initialize model state adjoint (Yb)*/
	double *Yb  = xNewZeroInit<double>(Ysize);

	/*Get final Ysize*/
	hdl_countdoubles = new CountDoublesFunctor();
	femmodel->Marshall(hdl_countdoubles);
	int Ysize_i= hdl_countdoubles->DoubleCount();
	delete hdl_countdoubles;

	/*Start tracing*/
	codi_global.start();

	/*Reverse dependent (f)*/
	hdl_regin = new RegisterInputFunctor(&codi_y_data);
	femmodel->Marshall(hdl_regin);
	delete hdl_regin;
	if(my_rank==0) for(int i=0; i < Xsize; i++) codi_global.registerInput(X[i]);
	SetControlInputsFromVectorx(femmodel,X);

	IssmDouble J     = 0.;
	int        count = 0;
	DataSet* dependent_objects=((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
	for(Object* & object:dependent_objects->objects){
		DependentObject* dep=(DependentObject*)object;
		IssmDouble       output_value = dep->GetValue();

		J += output_value;

		/*Keep track of output for printing*/
		Jlist[count] = output_value.getValue();
		count++;
	}
	Jlist[count] = J.getValue();
	_assert_(count == num_responses);

	codi_global.registerOutput(J);

	codi_global.stop();

	if(VerboseAutodiff())_printf0_("   CoDiPack fos_reverse\n");
	if(my_rank==0) codi_global.setGradient(0, 1.0);
	codi_global.evaluate();

	/*Initialize Xb and Yb*/
	double *Xb  = xNewZeroInit<double>(Xsize);
	codi_global.updateFullGradient(Xb, Xsize);
	codi_y_data.getFullGradient(Yb, Ysize);

	/*reverse loop for transient step (G)*/
	for(vector<int>::reverse_iterator iter = checkpoint_steps.rbegin(); iter != checkpoint_steps.rend(); iter++){

		/*Restore model from this step*/
		int reverse_step = *iter;
		femmodel->RestartAD(reverse_step);

		codi_y_data.clear();
		codi_global.start();

		/*Get new Ysize*/
		hdl_countdoubles = new CountDoublesFunctor();
		femmodel->Marshall(hdl_countdoubles);
		int Ysize_i= hdl_countdoubles->DoubleCount();
		delete hdl_countdoubles;

		/*We need to store the CoDiPack identifier here, since y is overwritten.*/
		hdl_regin = new RegisterInputFunctor(&codi_y_data);
		femmodel->Marshall(hdl_regin);
		delete hdl_regin;

		/*Tell codipack that X is the independent*/
		for(int i=0; i<Xsize; i++) codi_global.registerInput(X[i]);
		SetControlInputsFromVectorx(femmodel,X);

		/*Get New state*/
		for(int ii=0;ii<checkpoint_frequency;ii++){
			int        thisstep = reverse_step+ii;
			IssmDouble thistime = time_all[reverse_step+ii-1];
			IssmDouble thisdt   = dt_all[reverse_step+ii-1];
			femmodel->parameters->SetParam(thistime,TimeEnum);
			femmodel->parameters->SetParam(thisstep,StepEnum);
			femmodel->parameters->SetParam(thisdt,TimesteppingTimeStepEnum);

			if(VerboseSolution()){
				_printf0_("step "<<thisstep<<" ("<<ii+1<<"/"<<checkpoint_frequency<<") time [yr]: "\
						<<std::fixed<<std::setprecision(2)<<thistime/yts<< " (time step: " << thisdt/yts << ")\n");
			}

			transient_step(femmodel);

			/*Go through our dependent variables, and compute the response:*/
			DataSet* dependent_objects=((DataSetParam*)femmodel->parameters->FindParamObject(AutodiffDependentObjectsEnum))->value;
			for(Object* & object:dependent_objects->objects){
				DependentObject* dep=(DependentObject*)object;
				dep->RecordResponsex(femmodel);
			}

			/*First and last segment need special treatment*/
			if(thisstep==finalstep) break;
			if(reverse_step==1 && ii==checkpoint_frequency-2) break;
		}

		/*Register output*/
		hdl_regout = new RegisterOutputFunctor(&codi_y_data);
		femmodel->Marshall(hdl_regout);
		delete hdl_regout;

		/*stop tracing*/
		codi_global.stop();

		/*Reverse transient step (G)*/
		/* Using y_b here to seed the next reverse iteration there y_b is always overwritten*/
		codi_y_data.setFullGradient(Yb, Ysize);

		if(VerboseSolution()) _printf0_("computing gradient...\n");
		codi_global.evaluate();

		/* here we access the gradient data via the stored identifiers.*/
		codi_global.updateFullGradient(Xb, Xsize);
		codi_y_data.getFullGradient(Yb, Ysize);	// Yb is overwritten here.
	}

	/*Clear tape*/
	codi_global.clear();

	/*Broadcast gradient to other ranks (make sure to sum all gradients)*/
	ISSM_MPI_Allreduce(Xb,G,Xsize,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	#ifdef _ISSM_DEBUG_
	for(int i=0; i<Xsize; i++){
		if(xIsNan(Xb[i])) _error_("Found NaN in gradient at position "<<i);
		if(xIsInf(Xb[i])) _error_("Found Inf in gradient at position "<<i);
	}
	#endif

	/*Cleanup and return misfit*/
	xDelete<IssmDouble>(X);
	xDelete<double>(Xb);
	xDelete<double>(Yb);
	xDelete<int>(control_enum);
	return J.getValue();
}/*}}}*/
#endif
