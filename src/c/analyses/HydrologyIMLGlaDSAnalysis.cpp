#include <float.h> /* defines DBL_EPSILON*/
#include "./HydrologyIMLGlaDSAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include <mpi.h> 

/*Model processing*/
void HydrologyIMLGlaDSAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/

	/*retrieve some parameters: */
	int hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	if(hydrology_model!=HydrologyIMLGlaDSEnum) return;

	IoModelToConstraintsx(constraints,iomodel,"md.hydrology.spcphi",HydrologyIMLGlaDSAnalysisEnum,P1Enum);

	bool islakes;
	iomodel->FindConstant(&islakes,"md.hydrology.islakes");
	if(islakes){
		/*Fetch lake datalake*/
		IssmDouble* lake_data = NULL;
		IssmDouble* spcphi    = NULL;
		IssmDouble* phi       = NULL;
		iomodel->FetchData(&lake_data,NULL,NULL,"md.hydrology.lake_mask");
		iomodel->FetchData(&spcphi,NULL,NULL,"md.hydrology.spcphi");
		iomodel->FetchData(&phi,NULL,NULL,"md.initialization.hydraulic_potential");

		/*Create spcs from x,y,z, as well as the spc values on those spcs: */
		int count = constraints->Size();
		for(int i=0;i<iomodel->numberofvertices;i++){
			if(iomodel->my_vertices[i]){
				if(xIsNan(spcphi[i]) && lake_data[i]>0.){
					constraints->AddObject(new SpcDynamic(count+1,i+1,0,0., HydrologyIMLGlaDSAnalysisEnum));
					count++;
				}
			}
		}
	
		/*Cleanup*/
		xDelete<IssmDouble>(lake_data);
		xDelete<IssmDouble>(spcphi);
		xDelete<IssmDouble>(phi);
	}

}/*}}}*/
void HydrologyIMLGlaDSAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want IML-GlaDS?*/
	if(hydrology_model!=HydrologyIMLGlaDSEnum) return;

	/*Add channels?*/
	int Ka,La;
	int Kd,Ld;
	bool ischannels;
	IssmDouble* channelarea;
	IssmDouble* channeldischarge;
	iomodel->FindConstant(&ischannels,"md.hydrology.ischannels");
	iomodel->FetchData(&channelarea,&Ka,&La,"md.initialization.channelarea");
	iomodel->FetchData(&channeldischarge,&Kd,&Ld,"md.initialization.channel_discharge");
	if(ischannels){
		/*Get faces (edges in 2d)*/
		CreateFaces(iomodel);
		for(int i=0;i<iomodel->numberoffaces;i++){

			/*Get left and right elements*/
			int element=iomodel->faces[4*i+2]-1; //faces are [node1 node2 elem1 elem2]

			/*Now, if this element is not in the partition*/
			if(!iomodel->my_elements[element]) continue;

			/*Add channelarea and discharge from initialization if exists*/
			if(Ka==0 && Kd==0){
				loads->AddObject(new Channel(i+1,0.,0.,i,iomodel));
			}
			if(Ka==0 && Kd==iomodel->numberoffaces){
				loads->AddObject(new Channel(i+1,0.,channeldischarge[i],i,iomodel));
			}
			if(Ka==iomodel->numberoffaces && Kd==0){
				loads->AddObject(new Channel(i+1,channelarea[i],0.,i,iomodel));
			}
			if(Ka==iomodel->numberoffaces && Kd==iomodel->numberoffaces){
				loads->AddObject(new Channel(i+1,channelarea[i],channeldischarge[i],i,iomodel));
			}
			else{
				if (Ka!=0 && Ka!=iomodel->numberoffaces) _error_("Unknown dimension for channel area initialization.");
				if (Kd!=0 && Kd!=iomodel->numberoffaces) _error_("Unknown dimension for channel discharge initialization.");
			}
			iomodel->DeleteData(1,"md.initialization.channelarea");
            iomodel->DeleteData(1,"md.initialization.channel_discharge");
		}

	}
	xDelete<IssmDouble>(channelarea);
    xDelete<IssmDouble>(channeldischarge);


	/*Create discrete loads for Moulins*/
	CreateSingleNodeToElementConnectivity(iomodel);
	if(iomodel->domaintype!=Domain2DhorizontalEnum && iomodel->domaintype!=Domain3DsurfaceEnum) iomodel->FetchData(1,"md.mesh.vertexonbase");
	for(int i=0;i<iomodel->numberofvertices;i++){
		if (iomodel->domaintype!=Domain3DEnum){
			/*keep only this partition's nodes:*/
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}
		}
		else if(reCast<int>(iomodel->Data("md.mesh.vertexonbase")[i])){
			if(iomodel->my_vertices[i]){
				loads->AddObject(new Moulin(i+1,i,iomodel));
			}	
		}
	}
	iomodel->DeleteData(1,"md.mesh.vertexonbase");
	/*Deal with Neumann BC*/
	int M,N;
	int *segments = NULL;
	iomodel->FetchData(&segments,&M,&N,"md.mesh.segments");

	/*Check that the size seem right*/
	_assert_(N==3); _assert_(M>=3);
	for(int i=0;i<M;i++){
		if(iomodel->my_elements[segments[i*3+2]-1]){
		loads->AddObject(new Neumannflux(i+1,i,iomodel,segments));
		}
	}
	xDelete<int>(segments);
	
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	/*Fetch parameters: */
	int  hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want IML-GlaDS?*/
	if(hydrology_model!=HydrologyIMLGlaDSEnum) return;

	if(iomodel->domaintype==Domain3DEnum) iomodel->FetchData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
	::CreateNodes(nodes,iomodel,HydrologyIMLGlaDSAnalysisEnum,P1Enum);
	iomodel->DeleteData(2,"md.mesh.vertexonbase","md.mesh.vertexonsurface");
}/*}}}*/
int  HydrologyIMLGlaDSAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Fetch data needed: */
	int    hydrology_model;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
	int    meltflag;	
	iomodel->FindConstant(&meltflag,"md.hydrology.melt_flag");
	bool islakes;
	iomodel->FindConstant(&islakes,"md.hydrology.islakes");

	/*Now, do we really want IML-GlaDS?*/
	if(hydrology_model!=HydrologyIMLGlaDSEnum) return;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	iomodel->FetchDataToInput(inputs,elements,"md.geometry.thickness",ThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.base",BaseEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.geometry.bed",BedEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sealevel",SealevelEnum,0);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.geothermalflux",BasalforcingsGeothermalfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.basalforcings.groundedice_melting_rate",BasalforcingsGroundediceMeltingRateEnum);
	if(meltflag==2){
		iomodel->FetchDataToInput(inputs,elements,"md.smb.runoff",SmbRunoffEnum);
	}
	if(iomodel->domaintype!=Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonbase",MeshVertexonbaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonsurface",MeshVertexonsurfaceEnum);
	}
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ocean_levelset",MaskOceanLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonboundary",MeshVertexonboundaryEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.bump_height",HydrologyBumpHeightEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.sheet_conductivity",HydrologySheetConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.channel_conductivity",HydrologyChannelConductivityEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.englacial_void_ratio",HydrologyEnglacialVoidRatioEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.neumannflux",HydrologyNeumannfluxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.moulin_input",HydrologyMoulinInputEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.watercolumn",HydrologySheetThicknessEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.hydraulic_potential",HydraulicPotentialEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.hydrology.rheology_B_base",HydrologyRheologyBBaseEnum);
	if(islakes){
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.lake_mask",HydrologyLakeMaskEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.characteristic_outlet_length",HydrologyLakeOutletLengthEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.max_lake_area",HydrologyMaxLakeAreaEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.lake_channelQr",HydrologyLakeChannelQrEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.lake_outletQr",HydrologyLakeOutletQrEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.lake_height",HydrologyLakeHeightEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.hydrology.lake_Qin",HydrologyLakeQinEnum);
		/*Necessary to initialise lake model...*/
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.sheet_discharge",HydrologySheetDischargeEnum);
		}
	
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyEnum);
	if(iomodel->domaintype==Domain2DhorizontalEnum){
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vx",VxBaseEnum);
		iomodel->FetchDataToInput(inputs,elements,"md.initialization.vy",VyBaseEnum);
	}

	/*Friction*/
	FrictionUpdateInputs(elements, inputs, iomodel);

}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	/*retrieve some parameters: */
	int    hydrology_model;
	int    numoutputs;
	char** requestedoutputs = NULL;
	iomodel->FindConstant(&hydrology_model,"md.hydrology.model");

	/*Now, do we really want IML-GlaDS?*/
	if(hydrology_model!=HydrologyIMLGlaDSEnum) return;

	parameters->AddObject(new IntParam(HydrologyModelEnum,hydrology_model));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.pressure_melt_coefficient",HydrologyPressureMeltCoefficientEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.cavity_spacing",HydrologyCavitySpacingEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.ischannels",HydrologyIschannelsEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.melt_flag",HydrologyMeltFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.maxiter",HydrologyMaxiterEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.restol",HydrologyRestolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.reltol",HydrologyReltolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.abstol",HydrologyAbstolEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_sheet_width",HydrologyChannelSheetWidthEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_alpha",HydrologyChannelAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.channel_beta",HydrologyChannelBetaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_alpha",HydrologySheetAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.sheet_beta",HydrologySheetBetaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.omega",HydrologyOmegaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.istransition",HydrologyIsTransitionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.isincludesheetthickness",HydrologyIsIncludeSheetThicknessEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.creep_open_flag",HydrologyCreepOpenFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.islakes",HydrologyLakeFlagEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.hydrology.num_lakes",HydrologyNumLakesEnum));
	
	/*Friction*/
	FrictionUpdateParameters(parameters, iomodel);

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.hydrology.requested_outputs");
	parameters->AddObject(new IntParam(HydrologyNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(HydrologyRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.hydrology.requested_outputs");
}/*}}}*/

/*Finite Element Analysis*/
void           HydrologyIMLGlaDSAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           HydrologyIMLGlaDSAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* HydrologyIMLGlaDSAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* HydrologyIMLGlaDSAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
	_error_("Not implemented");
}/*}}}*/
ElementMatrix* HydrologyIMLGlaDSAnalysis::CreateKMatrix(Element* element){/*{{{*/
	/*Use the  GlaDS implementation*/
	return HydrologyGlaDSAnalysis::CreateKMatrix(element);
}/*}}}*/
ElementVector* HydrologyIMLGlaDSAnalysis::CreatePVector(Element* element){/*{{{*/
	/*Use the  GlaDS implementation*/
	return HydrologyGlaDSAnalysis::CreatePVector(element);
}/*}}}*/
void           HydrologyIMLGlaDSAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	/*if lakes are activated, update phi at the lake outlet boundary*/
	bool islakes;
	element->FindParam(&islakes,HydrologyLakeFlagEnum);
	/*update the element phi*/
	element->GetSolutionFromInputsOneDof(solution,HydraulicPotentialEnum);
	/*update the element phi at the lake outlet boundary*/
	if(!islakes) return;
	else {
		UpdateLakeOutletPhi(element);
	}

}/*}}}*/
void           HydrologyIMLGlaDSAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           HydrologyIMLGlaDSAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	
	/*if lakes are activated, update phi at the lake outlet boundary*/
	bool islakes;
	element->FindParam(&islakes,HydrologyLakeFlagEnum);
	/*update the element phi*/
	element->InputUpdateFromSolutionOneDof(solution,HydraulicPotentialEnum);
	/*update the element phi at the lake outlet boundary*/
	if(islakes){
		UpdateLakeOutletPhi(element);
	}
	/*Use the GlaDS implementation to update outputs*/
	HydrologyGlaDSAnalysis::UpdateOutputs(element);
}/*}}}*/
void           HydrologyIMLGlaDSAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	
	/*Use the GlaDS implementation to update constraints*/
	return HydrologyGlaDSAnalysis::UpdateConstraints(femmodel);
}/*}}}*/

/*IML-GlaDS specifics*/
void HydrologyIMLGlaDSAnalysis::UpdateLakeOutletPhi(FemModel* femmodel){/*{{{*/

	for(Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		UpdateLakeOutletPhi(element);
	}
}/*}}}*/

void HydrologyIMLGlaDSAnalysis::UpdateLakeOutletPhi(Element* element){/*{{{*/

	if(!element->IsAnyLake()){
		return;
	}
	else{
		/*Intermediaries*/
		IssmDouble zb, lh;
		IssmDouble lakeID;
		IssmDouble phi;

		/*Get number of vertices for this element*/
		int numvertices = element->GetNumberOfVertices();

		/*Initialise phi at the lake outlet*/
		IssmDouble *phiLO = xNew<IssmDouble>(numvertices);

		/*Retrieve all parameters*/
		IssmDouble rho_ice = element->FindParam(MaterialsRhoIceEnum);
		IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
		IssmDouble g = element->FindParam(ConstantsGEnum);
		Input *zb_input = element->GetInput(BaseEnum); 					_assert_(zb_input);
		Input *lh_input = element->GetInput(HydrologyLakeHeightEnum); 	_assert_(lh_input);
		Input *phi_input = element->GetInput(HydraulicPotentialEnum); 	_assert_(phi_input);
		Input *lakeID_input = element->GetInput(HydrologyLakeMaskEnum); _assert_(lakeID_input);
	
		/*Get values at gauss points*/
		Gauss *gauss = element->NewGauss();
		for (int in = 0; in < numvertices; in++){
			gauss->GaussVertex(in);

			/*Get input values at gauss points*/
			zb_input->GetInputValue(&zb, gauss);
			lh_input->GetInputValue(&lh, gauss);
			phi_input->GetInputValue(&phi, gauss);
			lakeID_input->GetInputValue(&lakeID, gauss);

			/*Update the lake outlet phi*/
			if (lakeID > 0.){
				phiLO[in] = rho_water * g * zb + rho_water * g * lh;
			}
			else{
				phiLO[in] = phi;
			}
		}
		element->AddInput(HydraulicPotentialEnum, phiLO, P1Enum);
		/*clean up*/
		xDelete<IssmDouble>(phiLO);
		delete gauss;
	}
	
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::SetChannelCrossSectionOld(FemModel* femmodel){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::SetChannelCrossSectionOld(femmodel);
}/*}}}*/

void HydrologyIMLGlaDSAnalysis::SetLakeOutletDischargeOld(FemModel* femmodel){/*{{{*/
	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	bool islakes;
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!islakes){
		return;
	}
	/*Initialise QrcOld vector*/
	int nv = femmodel->vertices->NumberOfVertices();
	Vector<IssmDouble>* QrcOld = new Vector<IssmDouble>(nv);
	
	if (!ischannels){
		/*Populate empty (because no channels) Qrc vector, sheet discharge
		will then be added later */
		QrcOld->Set(0.0);

	}
	else {
		/*Loop over all channels, and have them add their discharge*/
		for(int i=0; i<femmodel->loads->Size(); i++){
			if(femmodel->loads->GetEnum(i)==ChannelEnum){
				Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
				channel->AddDischargeToVector(QrcOld);
			}
		}
	}
	/*Assemble vector*/
	QrcOld->Assemble();

	/*Update inputs*/
	InputUpdateFromVectorx(femmodel, QrcOld, HydrologyLakeChannelQrOldEnum, VertexSIdEnum);
	delete QrcOld;

}/*}}}*/

void HydrologyIMLGlaDSAnalysis::UpdateLakeDepth(FemModel* femmodel){/*{{{*/

	bool islakes;
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!islakes) return;

	int numlakes;
	femmodel->parameters->FindParam(&numlakes,HydrologyNumLakesEnum);

	/*Initialise arrays for alllakes, skip index 0 (no lake)*/
	IssmDouble* dt             = xNewZeroInit<IssmDouble>(numlakes+1);
    IssmDouble* local_qr       = xNewZeroInit<IssmDouble>(numlakes+1);
    IssmDouble* total_qr       = xNewZeroInit<IssmDouble>(numlakes+1);
    IssmDouble* lake_Qin       = xNewZeroInit<IssmDouble>(numlakes+1);
	IssmDouble* lake_area_old  = xNewZeroInit<IssmDouble>(numlakes+1);
    IssmDouble* lake_height_old = xNewZeroInit<IssmDouble>(numlakes+1);
    IssmDouble* lake_height     = xNewZeroInit<IssmDouble>(numlakes+1);

	/* First, loop over all elements to gather data for each lake */
	for (Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		if(element->IsAllFloating() || !element->IsIceInElement()){
			continue;
		}
		
		int numvertices = element->GetNumberOfVertices();
		IssmDouble ele_dt = element->FindParam(TimesteppingTimeStepEnum);
		IssmDouble lc	  = element->FindParam(HydrologyChannelSheetWidthEnum);

		Input* le_input    	= element->GetInput(HydrologyLakeOutletLengthEnum); _assert_(le_input);
		Input* qin_input   	= element->GetInput(HydrologyLakeQinEnum); 			_assert_(qin_input);
		Input* qrc_input   	= element->GetInput(HydrologyLakeChannelQrEnum); 	_assert_(qrc_input);
		Input* qs_input    	= element->GetInput(HydrologySheetDischargeEnum); 	_assert_(qs_input);
		Input* lamax_input 	= element->GetInput(HydrologyMaxLakeAreaEnum); 		_assert_(lamax_input); 
		Input* lh_old_input = element->GetInput(HydrologyLakeHeightOldEnum); 	_assert_(lh_old_input);
		Input* lakeID_input = element->GetInput(HydrologyLakeMaskEnum); 		_assert_(lakeID_input);
		Input* phi_input    = element->GetInput(HydraulicPotentialEnum); 		_assert_(phi_input);
	
		/*Params for lake area scaling (not yet implemented)*/
        Input* la_old_input = NULL;
		/*Get values at gauss points*/
		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			IssmDouble le, qin, qrc, qs, lamax, lh_old;

			/*Read lakeID safely by using an intermediate Issmdouble variable.
			This prevents corruption from casting pointer types.*/
			IssmDouble lake_id_double;
			lakeID_input->GetInputValue(&lake_id_double,gauss);
			int lakeID = (int)lake_id_double;

			le_input->GetInputValue(&le,gauss);
			qin_input->GetInputValue(&qin,gauss);
			qrc_input->GetInputValue(&qrc,gauss);
			qs_input->GetInputValue(&qs,gauss);
			lamax_input->GetInputValue(&lamax,gauss);
			lh_old_input->GetInputValue(&lh_old,gauss);

			if(lakeID > 0){
				/*Aggregate Qr contribution, dividing by element connectivity to avoid double-counting*/
				/* See eq.4 in Hepburn et al., 2026 */
				IssmDouble connectivity = (IssmDouble)element->VertexConnectivity(iv);
				local_qr[lakeID] += qrc/connectivity;
				if (qrc > 0.){
					local_qr[lakeID] += qs*(le-lc)/connectivity; /*Subtract lc to avoid double counting the sheet contribution to the channel*/
				}
				else if (qrc < 0.){
					local_qr[lakeID] += -qs*(le-lc)/connectivity;
				}
				else {
					/* When qrc == 0 (no channels), determine flow direction from hydraulic potential difference */
                    /* Get hydraulic potential at current lake vertex */
					IssmDouble phi_lake;
					phi_input->GetInputValue(&phi_lake,gauss);

					/* Calculate the average hydraulic potential of non-lake vertices within this element*/
					IssmDouble phi_avg = 0.0;
					int non_lake_count = 0;
					Gauss* temp_gauss = element->NewGauss();
					for(int jv = 0; jv < numvertices; jv++){
						temp_gauss->GaussVertex(jv);
						IssmDouble temp_lake_id;
						lakeID_input->GetInputValue(&temp_lake_id, temp_gauss);
						if((int)temp_lake_id == 0){ //non-lake vertex
							IssmDouble phi_temp;
							phi_input->GetInputValue(&phi_temp, temp_gauss);
							phi_avg += phi_temp;
							non_lake_count++;
						}
					}
					delete temp_gauss;

					if(non_lake_count > 0){
						phi_avg /= non_lake_count;
						/*If surrounding potential > lake potential, flow is into the lake (negative)*/
						/*If surrounding potential < lake potential, flow is out of the lake (positive)*/
						if(phi_avg > phi_lake){
							local_qr[lakeID] += -qs*(le-lc)/connectivity; //Flow into the lake
						}
						else{
							local_qr[lakeID] += qs*(le-lc)/connectivity; //Flow out of the lake
						}
					} else {
						/* All vertices are lake vertices, assume outflow*/
						local_qr[lakeID] += qs*(le-lc)/connectivity;
					}

				}

				lake_area_old[lakeID] = lamax;
				lake_Qin[lakeID] = qin;
				lake_height_old[lakeID] = lh_old;
				dt[lakeID] = ele_dt;
			}
		}
		delete gauss;
	}

	/* Sum the local_qr contributions from all processors into total_qr*/
	ISSM_MPI_Allreduce(local_qr, total_qr, 				numlakes+1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, IssmComm::GetComm());

	/*MPI_MAX ensures the non-zero valu from the processor that owns a lake element propogates to all others*/
	ISSM_MPI_Allreduce(MPI_IN_PLACE, lake_area_old, 	numlakes+1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());	
	ISSM_MPI_Allreduce(MPI_IN_PLACE, lake_Qin, 			numlakes+1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());
	ISSM_MPI_Allreduce(MPI_IN_PLACE, lake_height_old, 	numlakes+1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());
	ISSM_MPI_Allreduce(MPI_IN_PLACE, dt, 				numlakes+1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, IssmComm::GetComm());

	/* Now claculate a new lake height for all lakes using the globally aggregated data*/
	for (int lake = 1; lake<=numlakes; lake++){
		
		IssmDouble A 		= 1.0 / lake_area_old[lake];
		IssmDouble B 		= lake_Qin[lake] - total_qr[lake];

		IssmDouble alpha 	= 0.0;
		IssmDouble beta 	= A * B;

		lake_height[lake] = ODE1(alpha, beta, lake_height_old[lake], dt[lake],2);
		
		if(lake_height[lake] < DBL_EPSILON) lake_height[lake] = DBL_EPSILON;
	}

	/* Write the calculated results back to the elements*/
	for (Object* & object : femmodel->elements->objects){
		Element* element=xDynamicCast<Element*>(object);

		if(element->IsAllFloating() || !element->IsIceInElement()){
			continue;
		}

		int numvertices = element->GetNumberOfVertices();
		IssmDouble* lh_new = xNew<IssmDouble>(numvertices);
		IssmDouble* Qr_out = xNew<IssmDouble>(numvertices);
		IssmDouble* la_new = xNew<IssmDouble>(numvertices); //redundant, but will eventually store variable lake area when implemented

		Input* lakeID_input = element->GetInput(HydrologyLakeMaskEnum); _assert_(lakeID_input);
		Input* lamax_input = element->GetInput(HydrologyMaxLakeAreaEnum); _assert_(lamax_input);

		Gauss* gauss=element->NewGauss();
		for(int iv=0;iv<numvertices;iv++){
			gauss->GaussVertex(iv);
			/*Read lakeID safely here too*/
			IssmDouble lake_id_double;
			lakeID_input->GetInputValue(&lake_id_double,gauss);
			int lakeID = (int)lake_id_double;

			if(lakeID > 0){
				/*use the pre-calculate values from our arrays*/
				lh_new[iv] = lake_height[lakeID];
				Qr_out[iv] = total_qr[lakeID]; // FIXME: this writes the total QR to all lake vertices, not just the outlet. Need to fix this later.

				/*Get max lake area*/
				IssmDouble lamax;
				lamax_input->GetInputValue(&lamax,gauss);
				la_new[iv] = lamax; 
			}
			else{
				lh_new[iv] = 0.0;
				Qr_out[iv] = 0.0;
				la_new[iv] = 0.0;
			}
		}
		delete gauss;

		/*save results*/
		element->AddInput(HydrologyLakeAreaEnum, la_new, P1Enum);
		element->AddInput(HydrologyLakeHeightEnum, lh_new, P1Enum);
		element->AddInput(HydrologyLakeOutletQrEnum, Qr_out, P1Enum);

		/* clean up */
		xDelete<IssmDouble>(lh_new);
		xDelete<IssmDouble>(Qr_out);
		xDelete<IssmDouble>(la_new);
	}

	/* clean up */
	xDelete<IssmDouble>(dt);
	xDelete<IssmDouble>(local_qr);
	xDelete<IssmDouble>(total_qr);
	xDelete<IssmDouble>(lake_Qin);
	xDelete<IssmDouble>(lake_area_old);
	xDelete<IssmDouble>(lake_height_old);
	xDelete<IssmDouble>(lake_height);

}/*}}}*/

void HydrologyIMLGlaDSAnalysis::UpdateSheetThickness(FemModel* femmodel){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::UpdateSheetThickness(femmodel);
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateSheetThickness(Element* element){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::UpdateSheetThickness(element);
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateEffectivePressure(FemModel* femmodel){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::UpdateEffectivePressure(femmodel);
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateEffectivePressure(Element* element){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::UpdateEffectivePressure(element);
}/*}}}*/
void HydrologyIMLGlaDSAnalysis::UpdateChannelCrossSection(FemModel* femmodel){/*{{{*/
	/*Use the  GlaDS implementation*/
	HydrologyGlaDSAnalysis::UpdateChannelCrossSection(femmodel);
}/*}}}*/

void HydrologyIMLGlaDSAnalysis::UpdateLakeOutletDischarge(FemModel* femmodel){/*{{{*/
	bool ischannels;
	femmodel->parameters->FindParam(&ischannels,HydrologyIschannelsEnum);
	bool islakes;
	femmodel->parameters->FindParam(&islakes,HydrologyLakeFlagEnum);
	if(!islakes){
		return;
	}
	int nv = femmodel->vertices->NumberOfVertices();
	Vector<IssmDouble>* Qrc = new Vector<IssmDouble>(nv);
	if (!ischannels){
		/*Fill empty Qrc vector, sheet discharge will be added*/
		Qrc->Set(0.0);
		/*Assemble vector*/
		Qrc->Assemble();
	}
	else {
		/*populate qrc vector*/
		for(int i=0; i<femmodel->loads->Size(); i++){
			if(femmodel->loads->GetEnum(i)==ChannelEnum){
				Channel* channel=(Channel*)femmodel->loads->GetObjectByOffset(i);
				channel->AddDischargeToVector(Qrc);
			}
		}
		/*Assemble vector*/
		Qrc->Assemble();
	}
	InputUpdateFromVectorx(femmodel, Qrc, HydrologyLakeChannelQrEnum, VertexSIdEnum);
	delete Qrc;	
}/*}}}*/
