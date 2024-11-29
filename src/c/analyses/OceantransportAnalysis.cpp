#include "./OceantransportAnalysis.h"
#include <math.h>
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../classes/Inputs/TransientInput.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void OceantransportAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void OceantransportAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void OceantransportAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,OceantransportAnalysisEnum,P1Enum);
}/*}}}*/
int  OceantransportAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 3;
}/*}}}*/
void OceantransportAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int  nature=0;
	bool isdakota=0;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Plug inputs into element:*/
	iomodel->FetchDataToInput(inputs,elements,"md.dsl.sea_water_pressure_at_sea_floor", OceantransportSpcbottompressureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.dsl.sea_surface_height_above_geoid",  OceantransportSpcdslEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.dsl.global_average_thermosteric_sea_level",OceantransportSpcstrEnum);

	/*Resolve mmes if we not running Dakota. Otherwise, Dakota will provide a modelid, which will be use to resolve Mme, 
	 *but it will be done in InputUpdateFromDakota:*/
	iomodel->FetchData(&isdakota,"md.qmu.isdakota");
	if (inputs->GetInputObjectEnum(OceantransportSpcbottompressureEnum)==DatasetInputEnum && isdakota){
		int modelid;

		/*retrieve model id: */
		iomodel->FetchData(&modelid,"md.dsl.modelid");

		/*replace dataset of forcings with only one, the modelid'th:*/
		MmeToInputFromIdx(inputs,elements,modelid,OceantransportSpcbottompressureEnum, P1Enum);
		MmeToInputFromIdx(inputs,elements,modelid,OceantransportSpcdslEnum, P1Enum);
		MmeToInputFromIdx(inputs,elements,modelid,OceantransportSpcstrEnum, P0Enum);
	}

	/*Initialize sea level cumulated sea level loads :*/
	iomodel->ConstantToInput(inputs,elements,0.,AccumulatedDeltaBottomPressureEnum,P1Enum);
	iomodel->ConstantToInput(inputs,elements,0.,OldAccumulatedDeltaBottomPressureEnum,P1Enum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.bottompressure",BottomPressureEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.dsl",DslEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.str",StrEnum);

}/*}}}*/
void OceantransportAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int dslmodel=0;
	int     numoutputs;
	char**  requestedoutputs = NULL;

	/*Deal with dsl multi-model ensembles: {{{*/
	iomodel->FetchData(&dslmodel,"md.dsl.model");
	if(dslmodel==2){
		IssmDouble modelid; 
		int nummodels;

		/*create double param, not int param, because Dakota will be updating it as 
		 * a double potentially: */
		iomodel->FetchData(&modelid,"md.dsl.modelid");
		parameters->AddObject(new DoubleParam(DslModelidEnum,modelid));
		parameters->AddObject(iomodel->CopyConstantObject("md.dsl.nummodels",DslNummodelsEnum));
		iomodel->FetchData(&nummodels,"md.dsl.nummodels");

		/*quick checks: */
		if(nummodels<=0)_error_("dslmme object in  md.dsl field should contain at least 1 ensemble model!");
		if(modelid<=0 || modelid>nummodels)_error_("modelid field in dslmme object of md.dsl field should be between 1 and the number of ensemble runs!");
	} /*}}}*/
	/*Requested outputs {{{*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.solidearth.requested_outputs");
	if(numoutputs)parameters->AddObject(new StringArrayParam(SealevelchangeRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.solidearth.requested_outputs");
	/*}}}*/

}/*}}}*/

/*Finite Element Analysis*/
void           OceantransportAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           OceantransportAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* OceantransportAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* OceantransportAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* OceantransportAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* OceantransportAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           OceantransportAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/

	/*retrieve bottom pressure, dsl and str from the spcs in our element:*/

	IssmDouble bp,dsl,str;
	int       *doflist = NULL;

	/*Fetch number of nodes and initialize values*/
	int         numnodes = element->GetNumberOfNodes();
	int         numdof   = numnodes*3;
	IssmDouble* values   = xNew<IssmDouble>(numdof);

	/*Get dof list and inputs */
	element->GetDofList(&doflist,NoneApproximationEnum,GsetEnum);
	Input* bp_input=element->GetInput(OceantransportSpcbottompressureEnum); _assert_(bp_input);
	Input* dsl_input=element->GetInput(OceantransportSpcdslEnum); _assert_(dsl_input);
	Input* str_input=element->GetInput(OceantransportSpcstrEnum); _assert_(str_input);

	/*Ok, we have the velocities in inputs, fill in solution */
	Gauss* gauss=element->NewGauss();
	for(int i=0;i<numnodes;i++){
		gauss->GaussVertex(i);
		bp_input->GetInputValue(&bp,gauss);
		dsl_input->GetInputValue(&dsl,gauss);
		str_input->GetInputValue(&str,gauss);
		values[i*3+0]=bp;
		values[i*3+1]=dsl;
		values[i*3+2]=str;
	}

	/*Add value to global vector*/
	solution->SetValues(numdof,doflist,values,INS_VAL);

	/*Free resources:*/
	delete gauss;
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);
}/*}}}*/
void           OceantransportAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           OceantransportAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

	int         i,domaintype;
	int*        doflist=NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();
	int numdof   = numnodes*3;

	/*Fetch dof list and allocate solution vectors*/
	element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numdof);
	IssmDouble* bp        = xNew<IssmDouble>(numnodes);
	IssmDouble* dsl        = xNew<IssmDouble>(numnodes);
	IssmDouble* str        = xNew<IssmDouble>(numnodes);
	IssmDouble  strmean;

	/*Use the dof list to index into the solution vector: */
	for(i=0;i<numdof;i++) values[i]=solution[doflist[i]];

	/*Retrieve bp,dsl and str:*/
	strmean=0;
	for(i=0;i<numnodes;i++){
		bp[i]=values[i*3+0];
		dsl[i]=values[i*3+1];
		str[i]=values[i*3+2];
		strmean+=str[i]/numnodes;

		/*Check solution*/
		if(xIsNan<IssmDouble>(bp[i])) _error_("NaN found in bottom pressure solution vector");
		if(xIsInf<IssmDouble>(bp[i])) _error_("Inf found in bottom pressure  solution vector");
		if(xIsNan<IssmDouble>(dsl[i])) _error_("NaN found in dsl solution vector");
		if(xIsInf<IssmDouble>(dsl[i])) _error_("Inf found in dsl solution vector");
		if(xIsNan<IssmDouble>(str[i])) _error_("NaN found in str solution vector");
		if(xIsInf<IssmDouble>(str[i])) _error_("Inf found in str solution vector");
	}

	/*Add bp, dsl and str as inputs to the tria element: */
	element->AddInput(BottomPressureEnum,bp,P1Enum);
	element->AddInput(DslEnum,dsl,P1Enum);
	element->AddInput(StrEnum,&strmean,P0Enum); 

	/*Free resources:*/
	xDelete<IssmDouble>(bp);
	xDelete<IssmDouble>(str);
	xDelete<IssmDouble>(dsl);
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);

}/*}}}*/
void           OceantransportAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
