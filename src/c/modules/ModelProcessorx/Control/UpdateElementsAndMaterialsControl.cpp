/*
 * UpdateElementsAndMaterialsControl:
 */

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../ModelProcessorx.h"

void	UpdateElementsAndMaterialsControl(Elements* elements,Parameters* parameters,Inputs* inputs,Materials* materials, IoModel* iomodel){
	/*Intermediary*/
	bool       control_analysis;
	int        M,N;
	int        control,cost_function,domaintype;
	int        num_controls,num_cost_functions;
	IssmDouble yts,scale;
	Element     *element          = NULL;
	Material    *material         = NULL;
	int         *control_enums    = NULL;
	char       **controls         = NULL;
	char       **cost_functions   = NULL;
	IssmDouble  *independent      = NULL;
	IssmDouble  *independents_min = NULL;
	IssmDouble  *independents_max = NULL;
	IssmDouble  *weights          = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	if(!control_analysis) return;

	/*Fetch parameters: */
	bool isautodiff;
	iomodel->FindConstant(&isautodiff,"md.autodiff.isautodiff");
	if(isautodiff){
		UpdateElementsAndMaterialsControlAD(elements,parameters,inputs,materials,iomodel);
		return;
	}

	/*Process controls and convert from string to enums*/
	iomodel->FindConstant(&num_controls,"md.inversion.num_control_parameters");
	iomodel->FindConstant(&controls,&num_controls,"md.inversion.control_parameters");
	if(num_controls<1) _error_("no controls found");
	control_enums=xNew<int>(num_controls);
	for(int i=0;i<num_controls;i++){
		control_enums[i]=StringToEnumx(controls[i]);
	}

	/*Process cost functions and convert from string to enums*/
	iomodel->FindConstant(&num_cost_functions,"md.inversion.num_cost_functions");
	iomodel->FindConstant(&cost_functions,&num_cost_functions,"md.inversion.cost_functions");
	if(num_cost_functions<1) _error_("No cost functions found");
	int* cost_function_enums=xNew<int>(num_cost_functions);
	for(int i=0;i<num_cost_functions;++i){
		cost_function_enums[i]=StringToEnumx(cost_functions[i]);
	}

	/*Fetch Observations and add to inputs*/
	iomodel->FindConstant(&domaintype,"md.mesh.domain_type");
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FetchData(&weights,&M,&N,"md.inversion.cost_functions_coefficients");

	/*Transpose weights for simplicity!*/
	if(M*N && N>1){
		IssmDouble* weights_transp = xNew<IssmDouble>(M*N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) weights_transp[j*M+i] = weights[i*N+j];
		xDelete<IssmDouble>(weights);
		weights = weights_transp;
	}

	if(M!=iomodel->numberofvertices && N!=num_cost_functions) _error_("not supported");
	for(int i=0;i<num_cost_functions;i++){
		cost_function=cost_function_enums[i];
		if(     cost_function==ThicknessAbsMisfitEnum) iomodel->FetchDataToInput(inputs,elements,"md.inversion.thickness_obs",InversionThicknessObsEnum);
		else if(cost_function==SurfaceAbsMisfitEnum)   iomodel->FetchDataToInput(inputs,elements,"md.inversion.surface_obs",InversionSurfaceObsEnum);
		else if(cost_function==RheologyBInitialguessMisfitEnum) iomodel->FetchDataToInput(inputs,elements,"md.materials.rheology_B",RheologyBInitialguessEnum);
		else if(cost_function==SurfaceAbsVelMisfitEnum
			  || cost_function==SurfaceRelVelMisfitEnum
			  || cost_function==SurfaceLogVelMisfitEnum
			  || cost_function==SurfaceLogVxVyMisfitEnum
			  || cost_function==SurfaceAverageVelMisfitEnum){
			iomodel->FetchDataToInput(inputs,elements,"md.inversion.vx_obs",InversionVxObsEnum);
			if(domaintype!=Domain2DverticalEnum) iomodel->FetchDataToInput(inputs,elements,"md.inversion.vy_obs",InversionVyObsEnum); 
		}
		for(Object* & object : elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			element->DatasetInputAdd(InversionCostFunctionsCoefficientsEnum,&weights[i*iomodel->numberofvertices],inputs,iomodel,M,1,1,cost_function,cost_function);
		}
	}
	parameters->AddObject(new IntParam(ControlInputSizeMEnum,iomodel->numberofvertices));
	xDelete<IssmDouble>(weights);

	/*Get controls*/
	iomodel->FetchData(&independents_min,&M,&N,"md.inversion.min_parameters");
	if(M!=iomodel->numberofvertices && N!=num_controls) _error_("not supported");
	iomodel->FetchData(&independents_max,&M,&N,"md.inversion.max_parameters");
	if(M!=iomodel->numberofvertices && N!=num_controls) _error_("not supported");

	/*Transpose weights for simplicity!*/
	if(M*N && N>1){
		IssmDouble* independents_min_transp = xNew<IssmDouble>(M*N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) independents_min_transp[j*M+i] = independents_min[i*N+j];
		xDelete<IssmDouble>(independents_min);
		independents_min = independents_min_transp;

		IssmDouble* independents_max_transp = xNew<IssmDouble>(M*N);
		for(int i=0;i<M;i++) for(int j=0;j<N;j++) independents_max_transp[j*M+i] = independents_max[i*N+j];
		xDelete<IssmDouble>(independents_max);
		independents_max = independents_max_transp;
	}

	int* M_all = xNew<int>(num_controls);
	int* N_all = xNew<int>(num_controls);
	int* Interp_all = xNew<int>(num_controls);

	int offset = 0;
	for(int i=0;i<num_controls;i++){
		control = control_enums[i];
		if(!IsInputEnum(control)) _error_("Only inputs can be parameters except if you use AD");
		scale   = 1.;

		switch(control){
			/*List of supported controls*/
			case BalancethicknessThickeningRateEnum:      iomodel->FetchData(&independent,&M,&N,"md.balancethickness.thickening_rate");scale = 1./yts; break; 
			case BalancethicknessSpcthicknessEnum:        iomodel->FetchData(&independent,&M,&N,"md.balancethickness.spcthickness");                   break; 
			case VxEnum:                                  iomodel->FetchData(&independent,&M,&N,"md.initialization.vx");scale = 1./yts;                break; 
			case VyEnum:                                  iomodel->FetchData(&independent,&M,&N,"md.initialization.vy");scale = 1./yts;                break; 
			case ThicknessEnum:                           iomodel->FetchData(&independent,&M,&N,"md.geometry.thickness");                              break; 
			case FrictionCoefficientEnum:                 iomodel->FetchData(&independent,&M,&N,"md.friction.coefficient");                            break; 
			case FrictionCEnum:									 iomodel->FetchData(&independent,&M,&N,"md.friction.C");				                            break; 
			case FrictionAsEnum:                          iomodel->FetchData(&independent,&M,&N,"md.friction.As");                                     break; 
			case BalancethicknessApparentMassbalanceEnum: iomodel->FetchData(&independent,&M,&N,"md.balancethickness.apparent_massbalance");           break; 
			case BalancethicknessOmegaEnum:               iomodel->FetchData(&independent,&M,&N,"md.balancethickness.omega");                          break; 
			case MaterialsRheologyBEnum:                  iomodel->FetchData(&independent,&M,&N,"md.materials.rheology_B");                            break; 
			/*Special cases*/
			case MaterialsRheologyBbarEnum:               iomodel->FetchData(&independent,&M,&N,"md.materials.rheology_B");                            break; 
			case DamageDbarEnum:                          iomodel->FetchData(&independent,&M,&N,"md.damage.D");                                        break; 
			default:
				_error_("Control " << EnumToStringx(control) << " not implemented yet");
		}

		/*Transient independents not supported outside of AD*/
		if(N!=1) _error_("Transient controls not supported yet");
		N_all[i] = N;
		M_all[i] = M;

		if(M==iomodel->numberofvertices){
			Interp_all[i] = P1Enum;
		}
		else if(M==iomodel->numberofelements){
			Interp_all[i] = P0Enum;
		}
		else{
			_error_("Control size not supported");
		}

		/*Special case if 3d*/
		if(iomodel->domaintype==Domain3DEnum){
			if(control==MaterialsRheologyBbarEnum) control=MaterialsRheologyBEnum;
			if(control==DamageDbarEnum)            control=DamageDEnum;
		}

		for(Object* & object : elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			element->ControlInputCreate(independent,&independents_min[offset],&independents_max[offset],inputs,iomodel,M,N,scale,control,i+1);
		}
		xDelete<IssmDouble>(independent);

		offset += M*N;
	}
	parameters->AddObject(new IntVecParam(ControlInputSizeMEnum,M_all,num_controls));
	parameters->AddObject(new IntVecParam(ControlInputSizeNEnum,N_all,num_controls));
	parameters->AddObject(new IntVecParam(ControlInputInterpolationEnum,Interp_all,num_controls));
	xDelete<int>(M_all);
	xDelete<int>(N_all);
	xDelete<int>(Interp_all);
	xDelete<IssmDouble>(independents_min);
	xDelete<IssmDouble>(independents_max);

	/*Free data: */
	for(int i=0;i<num_controls;i++){
		switch(control_enums[i]){
			/*List of supported controls*/
			case BalancethicknessThickeningRateEnum:      iomodel->DeleteData(1,"md.balancethickness.thickening_rate"); break;
			case BalancethicknessSpcthicknessEnum:        iomodel->DeleteData(1,"md.balancethickness.spcthickness"); break;
			case VxEnum:                                  iomodel->DeleteData(1,"md.initialization.vx"); break;
			case VyEnum:                                  iomodel->DeleteData(1,"md.initialization.vy"); break;
			case ThicknessEnum:                           iomodel->DeleteData(1,"md.geometry.thickness"); break;
			case FrictionCoefficientEnum:                 iomodel->DeleteData(1,"md.friction.coefficient"); break;
			case FrictionCEnum:			                   iomodel->DeleteData(1,"md.friction.C"); break;
			case FrictionAsEnum:                          iomodel->DeleteData(1,"md.friction.As"); break;
			case BalancethicknessApparentMassbalanceEnum: iomodel->DeleteData(1,"md.balancethickness.apparent_massbalance"); break;
			case BalancethicknessOmegaEnum:               iomodel->DeleteData(1,"md.balancethickness.omega"); break;
			case MaterialsRheologyBEnum:                  iomodel->DeleteData(1,"md.materials.rheology_B"); break;
			/*Special cases*/
			case MaterialsRheologyBbarEnum: iomodel->DeleteData(1,"md.materials.rheology_B"); break;
			case DamageDbarEnum:            iomodel->DeleteData(1,"md.damage.D");            break;
			default:
				_error_("Control " << EnumToStringx(control_enums[i]) << " not implemented yet");
		}
	}

	xDelete<int>(control_enums);
	xDelete<int>(cost_function_enums);
	for(int i=0;i<num_cost_functions;i++) xDelete<char>(cost_functions[i]);
	xDelete<char*>(cost_functions);
	for(int i=0;i<num_controls;i++) xDelete<char>(controls[i]);
	xDelete<char*>(controls);
}
void UpdateElementsAndMaterialsControlAD(Elements* elements,Parameters* parameters,Inputs* inputs,Materials* materials, IoModel* iomodel){

	#if defined(_HAVE_AD_)
	/*Intermediaries*/
	int          num_independent_objects,M,N;
	char       **names                = NULL;
	int         *types                = NULL;
	int         *control_sizes        = NULL;
	IssmDouble  *independent          = NULL;
	IssmDouble **independents_fullmin = NULL;
	IssmDouble **independents_fullmax = NULL;
	bool         control_analysis     = false;

	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");

	/*Now, return if no control*/
	if(!control_analysis) return;

	/*Step1: create controls (independents)*/
	iomodel->FetchData(&num_independent_objects,"md.autodiff.num_independent_objects"); _assert_(num_independent_objects>0); 
	iomodel->FetchMultipleData(&names,&M,"md.autodiff.independent_name"); _assert_(M==num_independent_objects);

	int* M_all = NULL;
	int* N_all = NULL;
	int* Interp_all = xNew<int>(num_independent_objects);

	/*create independent objects, and at the same time, fetch the corresponding independent variables, 
	 *and declare them as such in ADOLC: */
	iomodel->FetchMultipleData(&independents_fullmin,&M_all,&N_all,&M,"md.autodiff.independent_min_parameters"); _assert_(M==num_independent_objects);
	iomodel->FetchMultipleData(&independents_fullmax,NULL  ,NULL  ,&M,"md.autodiff.independent_max_parameters"); _assert_(M==num_independent_objects);
	iomodel->FetchMultipleData(&control_sizes,&M,"md.autodiff.independent_control_size");                        _assert_(M==num_independent_objects);

	for(int i=0;i<num_independent_objects;i++){

		/*Get field name and input Enum from independent name*/
		char* iofieldname  = NULL;
		int   input_enum;
		IssmDouble* independents_min = NULL;
		IssmDouble*	independents_max = NULL;

		/*Fetch required data*/
		FieldAndEnumFromCode(&input_enum,&iofieldname,names[i]);
		iomodel->FetchData(&independent,&M,&N,iofieldname);
		_assert_(independent && N==control_sizes[i]);
		xDelete<char>(iofieldname);

		independents_min = NULL; independents_min = xNew<IssmDouble>(M*N);
		independents_max = NULL; independents_max = xNew<IssmDouble>(M*N);
		for(int m=0;m<M;m++){
			for(int n=0;n<N;n++){
				independents_min[N*m+n]=independents_fullmin[i][N*m+n];
				independents_max[N*m+n]=independents_fullmax[i][N*m+n];
			}
		}

		if(IsInputEnum(input_enum)){

			/*remove last row if time series*/
			if(N!=1) M_all[i]=M-1;

			if(M_all[i]==iomodel->numberofvertices){
				Interp_all[i] = P1Enum;
			}
			else if(M_all[i]==iomodel->numberofelements){
				Interp_all[i] = P0Enum;
			}
			else{
				_error_("Control size not supported");
			}

			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->ControlInputCreate(independent,independents_min,independents_max,inputs,iomodel,M,N,1.,input_enum,i+1);
			}
		}
		else if(IsParamEnum(input_enum)){
			//_error_("not supported yet");
			Interp_all[i] = DummyEnum; //Placeholder
			parameters->AddObject(new ControlParam(independent,independents_min,independents_max,input_enum,M_all[i],N_all[i]));

			if(M!=1){
				_assert_(M==2); //TransientParam
				M_all[i]=M-1;
			}
		}
		xDelete<IssmDouble>(independent);
		xDelete<IssmDouble>(independents_min);
		xDelete<IssmDouble>(independents_max);
	}
	parameters->AddObject(new IntVecParam(ControlInputSizeNEnum,N_all,num_independent_objects));
	parameters->AddObject(new IntVecParam(ControlInputSizeMEnum,M_all,num_independent_objects));
	parameters->AddObject(new IntVecParam(ControlInputInterpolationEnum,Interp_all,num_independent_objects));

	/*cleanup*/
	for(int i=0;i<num_independent_objects;i++){
		xDelete<char>(names[i]);
		xDelete<IssmDouble>(independents_fullmin[i]);
		xDelete<IssmDouble>(independents_fullmax[i]);
	}
	xDelete<char*>(names);
	xDelete<int>(types);
	xDelete<int>(M_all);
	xDelete<int>(N_all);
	xDelete<int>(Interp_all);
	xDelete<IssmDouble*>(independents_fullmin);
	xDelete<IssmDouble*>(independents_fullmax);
	xDelete<int>(control_sizes);

	return;
#else 
	_error_("AD not compiled");
#endif
}
