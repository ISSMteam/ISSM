#include "./LoveAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void LoveAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
}/*}}}*/
void LoveAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
}/*}}}*/
void LoveAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
}/*}}}*/
int  LoveAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	_error_("not needed!");
}/*}}}*/
void LoveAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

}/*}}}*/
void LoveAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/
	IssmDouble* frequencies = NULL;
	IssmDouble* hypergeom_z = NULL;
	IssmDouble* hypergeom_table1 = NULL;
	IssmDouble* hypergeom_table2 = NULL;
	int         nfreq,nz, nalpha;
	iomodel->FetchData(&nfreq,"md.love.nfreq");
	iomodel->FetchData(&nz,"md.love.hypergeom_nz");
	iomodel->FetchData(&nalpha,"md.love.hypergeom_nalpha");
	iomodel->FetchData(&frequencies,NULL,NULL,"md.love.frequencies");
	iomodel->FetchData(&hypergeom_z,NULL,NULL,"md.love.hypergeom_z");
	iomodel->FetchData(&hypergeom_table1,NULL,NULL,"md.love.hypergeom_table1");
	iomodel->FetchData(&hypergeom_table2,NULL,NULL,"md.love.hypergeom_table2");

	parameters->AddObject(new DoubleVecParam(LoveFrequenciesEnum,frequencies,nfreq));
	parameters->AddObject(new DoubleVecParam(LoveHypergeomZEnum,hypergeom_z,nz));
	parameters->AddObject(new DoubleMatParam(LoveHypergeomTable1Enum,hypergeom_table1,nz,nalpha));
	parameters->AddObject(new DoubleMatParam(LoveHypergeomTable2Enum,hypergeom_table2,nz,nalpha));

	xDelete<IssmDouble>(frequencies);
	xDelete<IssmDouble>(hypergeom_z);
	xDelete<IssmDouble>(hypergeom_table1);
	xDelete<IssmDouble>(hypergeom_table2);

	parameters->AddObject(iomodel->CopyConstantObject("md.love.nfreq",LoveNfreqEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.sh_nmax",LoveShNmaxEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.sh_nmin",LoveShNminEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.g0",LoveG0Enum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.r0",LoveR0Enum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.mu0",LoveMu0Enum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.Gravitational_Constant",LoveGravitationalConstantEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.chandler_wobble",LoveChandlerWobbleEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.allow_layer_deletion",LoveAllowLayerDeletionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.underflow_tol",LoveUnderflowTolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.pw_threshold",LovePostWidderThresholdEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.min_integration_steps",LoveMinIntegrationStepsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.max_integration_dr",LoveMaxIntegrationdrEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.integration_scheme",LoveIntegrationSchemeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.istemporal",LoveIsTemporalEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.n_temporal_iterations",LoveNTemporalIterationsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.love_kernels",LoveKernelsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.forcing_type",LoveForcingTypeEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.inner_core_boundary",LoveInnerCoreBoundaryEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.core_mantle_boundary",LoveCoreMantleBoundaryEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.complex_computation",LoveComplexComputationEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.quad_precision",LoveQuadPrecisionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.debug",LoveDebugEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.equatorialmoi",RotationalEquatorialMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.polarmoi",RotationalPolarMoiEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.rotational.angularvelocity",RotationalAngularVelocityEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.lovenumbers.tk2secular",TidalLoveK2SecularEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.hypergeom_nz",LoveHypergeomNZEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.love.hypergeom_nalpha",LoveHypergeomNAlphaEnum));
}/*}}}*/
void LoveAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/

/*Finite Element Analysis*/
void           LoveAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_(" not needed!");
}/*}}}*/
void           LoveAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_(" not needed!");
}/*}}}*/
ElementVector* LoveAnalysis::CreateDVector(Element* element){/*{{{*/
	_error_(" not needed!");
}/*}}}*/
ElementMatrix* LoveAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_(" not needed!");
}/*}}}*/
ElementMatrix* LoveAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_(" not needed!");
}/*}}}*/
ElementVector* LoveAnalysis::CreatePVector(Element* element){/*{{{*/
	_error_("not supported");
}/*}}}*/
void           LoveAnalysis::GetB(IssmDouble* B,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	_error_("not supported");
}/*}}}*/
void           LoveAnalysis::GetBprime(IssmDouble* Bprime,Element* element,int dim,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	_error_("not supported");
}/*}}}*/
void           LoveAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           LoveAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           LoveAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	_error_("not supported");
}/*}}}*/
