/*!\file: CreateParametersControl.cpp
 * \brief driver for creating parameters dataset, for control analysis.
 */ 

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../ModelProcessorx.h"

void CreateParametersControl(Parameters* parameters,IoModel* iomodel,int solution_type){

	bool        control_analysis;
	int         inversiontype;
	int         nsteps;
	int         num_controls;
	int         num_costfunc;
	char**      controls      = NULL;
	int        *maxiter       = NULL;
	char**      cm_responses  = NULL;
	IssmDouble *cm_jump       = NULL;
	IssmDouble *optscal       = NULL;
	IssmDouble *control_scaling_factors = NULL;

	/*retrieve some parameters: */
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	iomodel->FindConstant(&inversiontype,"md.inversion.type");

	if(control_analysis){

		switch(inversiontype){
			case 0:/*Brent Search*/
			case 1:/*TAO*/
			case 2:/*M1QN3*/
			case 3:/*Validation*/
				  {
				/*How many controls and how many responses?*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.num_control_parameters",InversionNumControlParametersEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.num_cost_functions",InversionNumCostFunctionsEnum));

				/*recover controls and convert to Enums*/
				iomodel->FindConstant(&controls,&num_controls,"md.inversion.control_parameters");
				if(num_controls<1) _error_("no controls found");
				int* control_enums=xNew<int>(num_controls);
				for(int i=0;i<num_controls;i++){
					control_enums[i]=StringToEnumx(controls[i]);
					xDelete<char>(controls[i]);
				}
				xDelete<char*>(controls);
				parameters->AddObject(new IntVecParam(InversionControlParametersEnum,control_enums,num_controls));

				iomodel->FindConstant(&cm_responses,&num_costfunc,"md.inversion.cost_functions");
				if(num_costfunc<1) _error_ ("no cost functions found");
				int* costfunc_enums=xNew<int>(num_costfunc);
				for(int i=0;i<num_costfunc;i++){
					costfunc_enums[i]=StringToEnumx(cm_responses[i]);
					xDelete<char>(cm_responses[i]);
				}
				xDelete<char*>(cm_responses);
				parameters->AddObject(new IntVecParam(InversionCostFunctionsEnum,costfunc_enums,num_costfunc));

				xDelete<int>(control_enums);
				xDelete<int>(costfunc_enums);

				break;
				  }
			case 4:/*AD M1QN3*/
				  {
			 /*Intermediaries*/
			int     num_independent_objects,M;
			char**  names = NULL;

			/*this is done somewhere else*/
			parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.num_independent_objects",InversionNumControlParametersEnum));
			parameters->AddObject(iomodel->CopyConstantObject("md.autodiff.num_dependent_objects",InversionNumCostFunctionsEnum));

			/*Step 1: create controls (independents)*/
			iomodel->FetchData(&num_independent_objects,"md.autodiff.num_independent_objects");  _assert_(num_independent_objects>0);
			iomodel->FetchMultipleData(&names,&M,"md.autodiff.independent_name");                _assert_(M==num_independent_objects);
			iomodel->FetchMultipleData(&control_scaling_factors,&M,"md.autodiff.independent_scaling_factor"); _assert_(M==num_independent_objects);
			int* ind_enums=xNew<int>(num_independent_objects);
			for(int i=0;i<num_independent_objects;i++){
				ind_enums[i]=StringToEnumx(names[i]);
				xDelete<char>(names[i]);
			}
			xDelete<char*>(names);
			parameters->AddObject(new DoubleVecParam(InversionControlScalingFactorsEnum,control_scaling_factors,num_independent_objects));
			parameters->AddObject(new IntVecParam(InversionControlParametersEnum,ind_enums,num_independent_objects));
			xDelete<int>(ind_enums);	

			/*Step 2: create cost functions (dependent)*/
			iomodel->FindConstant(&cm_responses,&num_costfunc,"md.autodiff.dependent_object_names");
			if(num_costfunc<1) _error_ ("no cost functions found");
			int* costfunc_enums=xNew<int>(num_costfunc);
			for(int i=0;i<num_costfunc;i++){
				costfunc_enums[i]=StringToEnumx(cm_responses[i]);
				xDelete<char>(cm_responses[i]);
			}
			xDelete<char*>(cm_responses);
			parameters->AddObject(new IntVecParam(InversionCostFunctionsEnum,costfunc_enums,num_costfunc));
			xDelete<int>(costfunc_enums);

			break;
			}
			default:
				_error_("not supported");
		}

		/*Inversion type specifics*/
		switch(inversiontype){
			case 0:/*Brent Search*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.incomplete_adjoint",InversionIncompleteAdjointEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.nsteps",InversionNstepsEnum));
				iomodel->FetchData(&cm_jump,&nsteps,NULL,"md.inversion.step_threshold");
				iomodel->FetchData(&optscal,NULL,NULL,"md.inversion.gradient_scaling");
				iomodel->FetchData(&maxiter,NULL,NULL,"md.inversion.maxiter_per_step");
				parameters->AddObject(new DoubleMatParam(InversionGradientScalingEnum,optscal,nsteps,num_controls));
				parameters->AddObject(new DoubleVecParam(InversionStepThresholdEnum,cm_jump,nsteps));
				parameters->AddObject(new IntVecParam(InversionMaxiterPerStepEnum,maxiter,nsteps));
				break;
			case 1:/*TAO*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.incomplete_adjoint",InversionIncompleteAdjointEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.gatol",InversionGatolEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.grtol",InversionGrtolEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.gttol",InversionGttolEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxsteps",InversionMaxstepsEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxiter",InversionMaxiterEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.algorithm",InversionAlgorithmEnum));
				break;
			case 2:/*M1QN3*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.incomplete_adjoint",InversionIncompleteAdjointEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.dxmin",InversionDxminEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.dfmin_frac",InversionDfminFracEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.gttol",InversionGttolEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxsteps",InversionMaxstepsEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxiter",InversionMaxiterEnum));
				iomodel->FetchData(&control_scaling_factors,NULL,NULL,"md.inversion.control_scaling_factors");
				parameters->AddObject(new DoubleVecParam(InversionControlScalingFactorsEnum,control_scaling_factors,num_controls));
				break;
			case 3:/*Validation*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.incomplete_adjoint",InversionIncompleteAdjointEnum));
				iomodel->FetchData(&control_scaling_factors,NULL,NULL,"md.inversion.control_scaling_factors");
				parameters->AddObject(new DoubleVecParam(InversionControlScalingFactorsEnum,control_scaling_factors,num_controls));
				break;
			case 4:/*AD M1QN3*/
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.dxmin",InversionDxminEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.dfmin_frac",InversionDfminFracEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.gttol",InversionGttolEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxsteps",InversionMaxstepsEnum));
				parameters->AddObject(iomodel->CopyConstantObject("md.inversion.maxiter",InversionMaxiterEnum));
			break;
			default:
				_error_("not supported");
		}

		xDelete<int>(maxiter);
		xDelete<IssmDouble>(control_scaling_factors);
		iomodel->DeleteData(cm_jump,"md.inversion.step_threshold");
		iomodel->DeleteData(optscal,"md.inversion.gradient_scaling");
	}
}
