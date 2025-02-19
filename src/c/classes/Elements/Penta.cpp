/*!\file Penta.cpp
 * \brief: implementation of the Penta object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "../Inputs/PentaInput.h"
#include "../Inputs/ControlInput.h"
#include "../Inputs/TransientInput.h"
#include "../Inputs/DatasetInput.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES   6
#define NUMVERTICES2D 3

/*Constructors/destructor/copy*/
Penta::~Penta(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Penta::Penta(int penta_id, int penta_sid,int penta_lid,IoModel* iomodel,int nummodels)/*{{{*/
	:ElementHook(nummodels,penta_id,NUMVERTICES,iomodel){

	int penta_elements_ids[2];

	/*Checks in debugging mode*/
	_assert_(iomodel->Data("md.mesh.upperelements"));
	_assert_(iomodel->Data("md.mesh.lowerelements"));

	/*id: */
	this->id  = penta_id;
	this->sid = penta_sid;
	this->lid = penta_lid;

	/*surface and base*/
	this->isonsurface = false;
	this->isonbase    = false;

	/*Build neighbors list*/
	if (xIsNan<IssmDouble>(iomodel->Data("md.mesh.upperelements")[penta_sid]) || iomodel->Data("md.mesh.upperelements")[penta_sid]==-1.) penta_elements_ids[1]=this->id; //upper penta is the same penta
	else                                    penta_elements_ids[1]=reCast<int,IssmDouble>((iomodel->Data("md.mesh.upperelements")[penta_sid]));
	if (xIsNan<IssmDouble>(iomodel->Data("md.mesh.lowerelements")[penta_sid]) || iomodel->Data("md.mesh.lowerelements")[penta_sid]==-1.) penta_elements_ids[0]=this->id; //lower penta is the same penta
	else                                    penta_elements_ids[0]=reCast<int,IssmDouble>((iomodel->Data("md.mesh.lowerelements")[penta_sid]));
	this->InitHookNeighbors(penta_elements_ids);

	//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
	this->parameters=NULL;

	/*initialize pointers:*/
	this->nodes             = NULL;
	this->vertices          = NULL;
	this->material          = NULL;
	this->verticalneighbors = NULL;

	/*Only allocate pointer*/
	this->element_type_list=xNew<int>(nummodels);

	/*surface and base*/
	_assert_(iomodel->Data("md.mesh.vertexonsurface"));
	_assert_(iomodel->Data("md.mesh.vertexonbase"));
	this->isonsurface = false;
	this->isonbase    = false;
	IssmDouble sum = 0.;
	for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonsurface")[reCast<int>(iomodel->elements[(penta_id-1)*NUMVERTICES+i])-1];
	_assert_(sum>=0 && sum<4);
	if(sum>2.5) this->isonsurface = true;
	sum = 0.;
	for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonbase")[reCast<int>(iomodel->elements[(penta_id-1)*NUMVERTICES+i])-1];
	_assert_(sum>=0 && sum<4);
	if(sum>2.5) this->isonbase = true;
}
/*}}}*/
Object* Penta::copy() {/*{{{*/

	int i;
	Penta* penta=NULL;

	penta=new Penta();

	//deal with PentaRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		penta->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++) {
			if (this->element_type_list[i]) penta->element_type_list[i]=this->element_type_list[i];
			else penta->element_type_list[i] = 0;
		}
	}
	else penta->element_type_list = NULL;
	penta->element_type=this->element_type;
	penta->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		penta->hnodes=xNew<Hook*>(penta->numanalyses);
		for(i=0;i<penta->numanalyses;i++){
			if (this->hnodes[i]) penta->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else penta->hnodes[i] = NULL;
		}
	}
	else penta->hnodes = NULL;

	penta->hvertices = (Hook*)this->hvertices->copy();
	penta->hmaterial = (Hook*)this->hmaterial->copy();
	if (this->hneighbors) penta->hneighbors = (Hook*)(this->hneighbors->copy());
	else penta->hneighbors = NULL;

	/*deal with Tria fields: */
	penta->id  = this->id;
	penta->sid = this->sid;
	penta->lid = this->lid;
	penta->isonbase  = this->isonbase;
	penta->isonsurface  = this->isonsurface;

	/*point parameters: */
	penta->parameters=this->parameters;

	/*recover objects: */
	if (this->nodes) {
		unsigned int num_nodes = 6;
		penta->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
		for(i=0;i<num_nodes;i++) if(this->nodes[i]) penta->nodes[i]=this->nodes[i]; else penta->nodes[i] = NULL;
	}
	else penta->nodes = NULL;

	penta->vertices = (Vertex**)this->hvertices->deliverp();
	penta->material = (Material*)this->hmaterial->delivers();
	penta->verticalneighbors = (Penta**)this->hneighbors->deliverp();

	return penta;

}
/*}}}*/
void Penta::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = PentaEnum;
   marshallhandle->call(object_enum);
	marshallhandle->call(this->isonsurface);
	marshallhandle->call(this->isonbase);

	/*Call parent classes: */
	ElementHook::Marshall(marshallhandle);
	Element::MarshallElement2(marshallhandle,this->numanalyses);
	PentaRef::Marshall(marshallhandle);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();
	verticalneighbors = (Penta**)this->hneighbors->deliverp();
}/*}}}*/

/*Other*/
void       Penta::AddBasalInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	_assert_(this->inputs);
	if(!IsOnBase()) return;
	else{
		if(interpolation_enum==P1Enum || interpolation_enum==P1DGEnum){
			IssmDouble extrudedvalues[NUMVERTICES];
			for(int i=0;i<NUMVERTICES2D;i++){
				extrudedvalues[i]=values[i];
				extrudedvalues[i+NUMVERTICES2D]=values[i];
			}
			Penta* penta=this;
			for(;;){
				penta->AddInput(input_enum,&extrudedvalues[0],interpolation_enum);
				if (penta->IsOnSurface()) break;
				penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
			}
		}
		else if(interpolation_enum==P0Enum){
			Penta* penta=this;
			for(;;){
				penta->AddInput(input_enum,&values[0],interpolation_enum);
				if (penta->IsOnSurface()) break;
				penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
			}
		}
		else _error_("Interpolation "<<EnumToStringx(interpolation_enum)<<" not implemented yet");
	}

}
/*}}}*/
void       Penta::AddInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/**/
	int vertexlids[NUMVERTICES];

	/*Call inputs method*/
	_assert_(this->inputs);
	switch(interpolation_enum){
		case P1Enum:
			for(int i=0;i<NUMVERTICES;i++) vertexlids[i]=this->vertices[i]->lid;
			inputs->SetPentaInput(input_enum,interpolation_enum,NUMVERTICES,vertexlids,values);
			break;
		case P1DGEnum:
			inputs->SetPentaInput(input_enum,interpolation_enum,this->lid,NUMVERTICES,values);
			break;
		default:
			inputs->SetPentaInput(input_enum,interpolation_enum,this->lid,this->GetNumberOfNodes(interpolation_enum),values);
	}

}
/*}}}*/
void       Penta::AddControlInput(int input_enum,Inputs* inputs,IoModel* iomodel,IssmDouble* values,IssmDouble* values_min,IssmDouble* values_max, int interpolation_enum,int id){/*{{{*/

	/*Intermediaries*/
	int vertexlids[NUMVERTICES];

	_assert_(iomodel->elements);
	for(int i=0;i<NUMVERTICES;i++){
		int vertexid =reCast<int>(iomodel->elements[NUMVERTICES*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
		vertexlids[i]=iomodel->my_vertices_lids[vertexid-1];
	}

	/*Create Control Input*/
	inputs->SetControlInput(input_enum,PentaInputEnum,interpolation_enum,id);
	ControlInput* control_input = inputs->GetControlInput(input_enum); _assert_(input_enum);

	/*Call inputs method*/
	switch(interpolation_enum){
		case P0Enum:
			control_input->SetControl(interpolation_enum,1,&this->lid,values,values_min,values_max);
			break;
		case P1Enum:
			control_input->SetControl(interpolation_enum,NUMVERTICES,&vertexlids[0],values,values_min,values_max);
			break;
		default:
			_error_("Cannot add \""<<EnumToStringx(input_enum)<<"\" interpolation "<<EnumToStringx(interpolation_enum)<<" not supported");
	}

}
/*}}}*/
void       Penta::DatasetInputCreate(IssmDouble* array,int M,int N,int* individual_enums,int num_inputs,Inputs* inputs,IoModel* iomodel,int input_enum){/*{{{*/

	/*Intermediaries*/
	int        vertexsids[NUMVERTICES];
	int        vertexlids[NUMVERTICES];
	IssmDouble nodeinputs[NUMVERTICES];

	/*Some sanity checks*/
	if(num_inputs<1)                 _error_("Cannot create a DatasetInput of size <1");
	if(M!=iomodel->numberofvertices) _error_("Input size not supported yet");
	if(N!=num_inputs)                _error_("Sizes are not consistent");

	/*Get indices*/
	_assert_(iomodel->elements);
	for(int i=0;i<NUMVERTICES;i++){
		vertexsids[i] = reCast<int>(iomodel->elements[NUMVERTICES*this->Sid()+i])-1;
		vertexlids[i] = iomodel->my_vertices_lids[vertexsids[i]];
	}

	/*Create inputs and add to DataSetInput*/
	for(int i=0;i<num_inputs;i++){
		for(int j=0;j<NUMVERTICES;j++) nodeinputs[j]=array[vertexsids[j]*N+i];
		inputs->SetPentaDatasetInput(input_enum,individual_enums[i],P1Enum,NUMVERTICES,vertexlids,nodeinputs);
	}
}
/*}}}*/
void       Penta::BasalNodeIndices(int* pnumindices,int** pindices,int finiteelement){/*{{{*/

	PentaRef::BasalNodeIndices(pnumindices,pindices,finiteelement);

}
/*}}}*/
void       Penta::AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part){/*{{{*/
	_error_("Not supported yet!");
}
/*}}}*/
void       Penta::CalvingRateVonmises(){/*{{{*/

	if(!this->IsOnBase()) return;
	this->ComputeSigmaVM();

	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  vx,vy;
	IssmDouble  sigma_vm,sigma_max,sigma_max_floating,sigma_max_grounded;
	IssmDouble  groundedice,bed,sealevel;

	/*Depth average velocity for stress calculation*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input = this->GetInput(VxAverageEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyAverageEnum); _assert_(vy_input);
	Input* gr_input = this->GetInput(MaskOceanLevelsetEnum); _assert_(gr_input);
	Input* bs_input = this->GetInput(BaseEnum);                    _assert_(bs_input);
	Input* smax_fl_input = this->GetInput(CalvingStressThresholdFloatingiceEnum); _assert_(smax_fl_input);
	Input* smax_gr_input = this->GetInput(CalvingStressThresholdGroundediceEnum); _assert_(smax_gr_input);
	Input* sl_input  = this->GetInput(SealevelEnum); _assert_(sl_input);
	Input* sigma_vm_input = this->GetInput(SigmaVMEnum); _assert_(sigma_vm_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<3;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components and thickness*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		sigma_vm_input->GetInputValue(&sigma_vm,&gauss);
		gr_input->GetInputValue(&groundedice,&gauss);
		bs_input->GetInputValue(&bed,&gauss);
		smax_fl_input->GetInputValue(&sigma_max_floating,&gauss);
		smax_gr_input->GetInputValue(&sigma_max_grounded,&gauss);
		sl_input->GetInputValue(&sealevel,&gauss);

		/*Tensile stress threshold*/
		if(groundedice<0)
		 sigma_max = sigma_max_floating;
		else
		 sigma_max = sigma_max_grounded;

		/*Assign values*/
		if(bed>sealevel){
			calvingrate[iv] = 0.;
		}
		else{
			calvingrate[iv] = sqrt(vx*vx+vy*vy)*sigma_vm/sigma_max;
		}
	}

	/*Add input*/
	this->AddBasalInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->CalvingRateToVector();

	/*Extrude*/
	this->InputExtrude(CalvingCalvingrateEnum,-1);
	this->InputExtrude(CalvingratexEnum,-1);
	this->InputExtrude(CalvingrateyEnum,-1);
}
/*}}}*/
void       Penta::CalvingRateLevermann(){/*{{{*/

	IssmDouble  strainparallel;
	IssmDouble  propcoeff;
	IssmDouble  strainperpendicular;
	IssmDouble  calvingrate[NUMVERTICES];

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input=this->GetInput(VxEnum);																		_assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);																		_assert_(vy_input);
	Input* strainparallel_input=this->GetInput(StrainRateparallelEnum);								_assert_(strainparallel_input);
	Input* strainperpendicular_input=this->GetInput(StrainRateperpendicularEnum);              _assert_(strainperpendicular_input);
	Input* levermanncoeff_input=this->GetInput(CalvinglevermannCoeffEnum);                     _assert_(levermanncoeff_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		strainparallel_input->GetInputValue(&strainparallel,&gauss);
		strainperpendicular_input->GetInputValue(&strainperpendicular,&gauss);
		levermanncoeff_input->GetInputValue(&propcoeff,&gauss);

		/*Calving rate proportionnal to the positive product of the strain rate along the ice flow direction and the strain rate perpendicular to the ice flow */
		calvingrate[iv]=propcoeff*strainparallel*strainperpendicular;
		if(calvingrate[iv]<0){
			calvingrate[iv]=0;
		}
	}

	/*Add input*/
	this->AddBasalInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->CalvingRateToVector();
}/*}}}*/
void       Penta::CalvingFluxLevelset(){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)){
		IssmDouble flux_per_area=0;
		this->AddInput(CalvingFluxLevelsetEnum,&flux_per_area,P0Enum);
	}
	else{
		int               index1,index2;
		const IssmPDouble epsilon = 1.e-15;
		IssmDouble        s1,s2;
		IssmDouble        gl[NUMVERTICES];
		IssmDouble        xyz_front[2][3];

		IssmDouble  xyz_list[NUMVERTICES][3];
      ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

		/*Recover parameters and values*/
		Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

		/*Be sure that values are not zero*/
		if(gl[0]==0.) gl[0]=gl[0]+epsilon;
		if(gl[1]==0.) gl[1]=gl[1]+epsilon;
		if(gl[2]==0.) gl[2]=gl[2]+epsilon;

		int pt1 = 0;
		int pt2 = 1;
		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

			/*Portion of the segments*/
			s1=gl[2]/(gl[2]-gl[1]);
			s2=gl[2]/(gl[2]-gl[0]);
			if(gl[2]<0.){
				pt1 = 1; pt2 = 0;
			}
			xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
			xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
			xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
			xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
			xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
			xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

			/*Portion of the segments*/
			s1=gl[0]/(gl[0]-gl[1]);
			s2=gl[0]/(gl[0]-gl[2]);
			if(gl[0]<0.){
				pt1 = 1; pt2 = 0;
			}

			xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
			xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
			xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
			xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
			xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
			xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

			/*Portion of the segments*/
			s1=gl[1]/(gl[1]-gl[0]);
			s2=gl[1]/(gl[1]-gl[2]);
			if(gl[1]<0.){
				pt1 = 1; pt2 = 0;
			}

			xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
			xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
			xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
			xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
			xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
			xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
		}
		else{
			_error_("case not possible");
		}

		/*Some checks in debugging mode*/
		_assert_(s1>=0 && s1<=1.);
		_assert_(s2>=0 && s2<=1.);

		/*Get normal vector*/
		IssmDouble normal[3];
		this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
		normal[0] = -normal[0];
		normal[1] = -normal[1];

		/*Get inputs*/
		IssmDouble flux = 0.;
		IssmDouble area = 0.;
		IssmDouble calvingratex,calvingratey,thickness,Jdet,flux_per_area;
		IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
		Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* calvingratex_input=NULL;
		Input* calvingratey_input=NULL;
		calvingratex_input=this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
		calvingratey_input=this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);

		/*Start looping on Gaussian points*/
		Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
		while(gauss->next()){
			thickness_input->GetInputValue(&thickness,gauss);
			calvingratex_input->GetInputValue(&calvingratex,gauss);
			calvingratey_input->GetInputValue(&calvingratey,gauss);
			this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

			flux += rho_ice*Jdet*gauss->weight*thickness*(calvingratex*normal[0] + calvingratey*normal[1]);
			area += Jdet*gauss->weight*thickness;

			flux_per_area=flux/area;
		}

		this->AddInput(CalvingFluxLevelsetEnum,&flux_per_area,P0Enum);

		/*Clean up and return*/
		delete gauss;
	}
}
/*}}}*/
void       Penta::CalvingMeltingFluxLevelset(){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)){
		IssmDouble flux_per_area=0;
		this->AddInput(CalvingMeltingFluxLevelsetEnum,&flux_per_area,P0Enum);
	}
	else{
		int               index1,index2;
		const IssmPDouble epsilon = 1.e-15;
		IssmDouble        s1,s2;
		IssmDouble        gl[NUMVERTICES];
		IssmDouble        xyz_front[2][3];

		IssmDouble  xyz_list[NUMVERTICES][3];
      ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

		/*Recover parameters and values*/
		Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

		/*Be sure that values are not zero*/
		if(gl[0]==0.) gl[0]=gl[0]+epsilon;
		if(gl[1]==0.) gl[1]=gl[1]+epsilon;
		if(gl[2]==0.) gl[2]=gl[2]+epsilon;

		int pt1 = 0;
		int pt2 = 1;
		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

			/*Portion of the segments*/
			s1=gl[2]/(gl[2]-gl[1]);
			s2=gl[2]/(gl[2]-gl[0]);
			if(gl[2]<0.){
				pt1 = 1; pt2 = 0;
			}
			xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
			xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
			xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
			xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
			xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
			xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

			/*Portion of the segments*/
			s1=gl[0]/(gl[0]-gl[1]);
			s2=gl[0]/(gl[0]-gl[2]);
			if(gl[0]<0.){
				pt1 = 1; pt2 = 0;
			}

			xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
			xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
			xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
			xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
			xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
			xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

			/*Portion of the segments*/
			s1=gl[1]/(gl[1]-gl[0]);
			s2=gl[1]/(gl[1]-gl[2]);
			if(gl[1]<0.){
				pt1 = 1; pt2 = 0;
			}

			xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
			xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
			xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
			xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
			xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
			xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
		}
		else{
			_error_("case not possible");
		}

		/*Some checks in debugging mode*/
		_assert_(s1>=0 && s1<=1.);
		_assert_(s2>=0 && s2<=1.);

		/*Get normal vector*/
		IssmDouble normal[3];
		this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
		normal[0] = -normal[0];
		normal[1] = -normal[1];

		/*Get inputs*/
		IssmDouble flux = 0.;
		IssmDouble area = 0.;
		IssmDouble calvingratex,calvingratey,vx,vy,vel,meltingrate,meltingratex,meltingratey,thickness,Jdet,flux_per_area;
		IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
		Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* calvingratex_input=NULL;
		Input* calvingratey_input=NULL;
		Input* vx_input=NULL;
		Input* vy_input=NULL;
		Input* meltingrate_input=NULL;
		calvingratex_input=this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
		calvingratey_input=this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
		meltingrate_input=this->GetInput(CalvingMeltingrateEnum); _assert_(meltingrate_input);

		/*Start looping on Gaussian points*/
		Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
		while(gauss->next()){
			thickness_input->GetInputValue(&thickness,gauss);
			calvingratex_input->GetInputValue(&calvingratex,gauss);
			calvingratey_input->GetInputValue(&calvingratey,gauss);
			vx_input->GetInputValue(&vx,gauss);
			vy_input->GetInputValue(&vy,gauss);
			vel=vx*vx+vy*vy;
			meltingrate_input->GetInputValue(&meltingrate,gauss);
			meltingratex=meltingrate*vx/(sqrt(vel)+1.e-14);
			meltingratey=meltingrate*vy/(sqrt(vel)+1.e-14);
			this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

			flux += rho_ice*Jdet*gauss->weight*thickness*((calvingratex+meltingratex)*normal[0] + (calvingratey+meltingratey)*normal[1]);
			area += Jdet*gauss->weight*thickness;

			flux_per_area=flux/area;
		}

		this->AddInput(CalvingMeltingFluxLevelsetEnum,&flux_per_area,P0Enum);

		/*Clean up and return*/
		delete gauss;
	}
}
/*}}}*/
void       Penta::CalvingRateCalvingMIP(){/*{{{*/

	if(!this->IsOnBase()) return;

	IssmDouble  calvingrate[NUMVERTICES];
	int			experiment = 1;  /* exp:1 by default */
	int         dim, domaintype;
	IssmDouble	vx, vy, vel, c, wrate;
	IssmDouble  time, groundedice, yts;

	/*Get problem dimension and whether there is moving front or not*/
	this->FindParam(&domaintype,DomainTypeEnum);
	this->FindParam(&time,TimeEnum);
	this->FindParam(&yts,ConstantsYtsEnum);

	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}
	if(dim==1) _error_("not implemented in 1D...");

	/*Depth average velocity for stress calculation*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Retrieve all inputs and parameters we will need*/
	Input *vx_input      = this->GetInput(VxAverageEnum);                                _assert_(vx_input);
	Input *vy_input      = this->GetInput(VyAverageEnum);                                _assert_(vy_input);
	Input *wrate_input   = this->GetInput(CalvingAblationrateEnum);               _assert_(wrate_input); 
	Input* gr_input      = this->GetInput(MaskOceanLevelsetEnum);						_assert_(gr_input);

	/* Use which experiment: use existing Enum */
	this->FindParam(&experiment, CalvingUseParamEnum);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<3;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components */
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=sqrt(vx*vx+vy*vy)+1.e-14;

		/* no calving for grounded ice in EXP4 */
		gr_input->GetInputValue(&groundedice,&gauss);

		switch (experiment) { 
			case 1:
			case 3:
				/* Exp 1 and 3: set c=v-wrate, wrate=0, so that w=0 */
				wrate = 0.0;
				break;
			case 2:
				/* Exp 2: set c=v-wrate*/
				wrate = -300*sin(2.0*M_PI*time/yts/1000)/yts;  // m/a -> m/s
				break;
			case 4:
				/* Exp 4: set c=v-wrate(given), for the first 500 years, then c=0 for the second 500 years*/
				if((groundedice<0) && (time<=500.0*yts)) {
				//	wrate_input->GetInputValue(&wrate,&gauss);
					wrate = -750*sin(2.0*M_PI*time/yts/1000)/yts;  // m/a -> m/s
				}
				else {
					/* no calving on the grounded ice*/
					wrate = vel;
				}
				break;
			default:
				_error_("The experiment is not supported yet!");
		}

		calvingrate[iv] = vel - wrate;
	}
	/*Add input*/
	this->AddBasalInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->CalvingRateToVector();

	/*Extrude*/
	this->InputExtrude(CalvingCalvingrateEnum,-1);
	this->InputExtrude(CalvingratexEnum,-1);
	this->InputExtrude(CalvingrateyEnum,-1);
}
/*}}}*/
void       Penta::ComputeBasalStress(void){/*{{{*/

	_error_("not implemented (needs to be redone)");
	int         i,j;
	int         dofv[3]={0,1,2};
	int         dofp[1]={3};
	int         analysis_type,approximation;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_tria[3][3];
	IssmDouble  rho_ice,gravity,FSreconditioning;
	IssmDouble  pressure,viscosity,Jdet2d;
	IssmDouble  bed_normal[3];
	IssmDouble  basalforce[3] = {0.};
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble  stresstensor[6]={0.0};
	IssmDouble  sigma_xx,sigma_yy,sigma_zz;
	IssmDouble  sigma_xy,sigma_xz,sigma_yz;
	IssmDouble  surface=0,value=0;
	GaussPenta* gauss;

	/*retrive parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	this->Element::GetInputValue(&approximation,ApproximationEnum);

	/*Check analysis_types*/
	if (analysis_type!=StressbalanceAnalysisEnum) _error_("Not supported yet!");
	if (approximation!=FSApproximationEnum) _error_("Not supported yet!");

	/*retrieve some parameters: */
	this->parameters->FindParam(&FSreconditioning,StressbalanceFSreconditioningEnum);

	if(!IsOnBase()){
		//put zero
		//sigma_b->SetValue(id-1,0.0,INS_VAL);
		return;
	}

	/*recovre material parameters: */
	rho_ice=FindParam(MaterialsRhoIceEnum);
	gravity=FindParam(ConstantsGEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	for(i=0;i<3;i++) for(j=0;j<3;j++) xyz_list_tria[i][j]=xyz_list[i][j];

	/*Retrieve all inputs we will be needing: */
	Input* pressure_input=this->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=this->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);             _assert_(vy_input);
	Input* vz_input=this->GetInput(VzEnum);             _assert_(vz_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussPenta(0,1,2,2);
	while(gauss->next()){

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],gauss,vx_input,vy_input,vz_input);
		pressure_input->GetInputValue(&pressure,gauss);

		/*Compute Stress*/
		sigma_xx=2*viscosity*epsilon[0]-pressure*FSreconditioning; // sigma = nu eps - pressure
		sigma_yy=2*viscosity*epsilon[1]-pressure*FSreconditioning;
		sigma_zz=2*viscosity*epsilon[2]-pressure*FSreconditioning;
		sigma_xy=2*viscosity*epsilon[3];
		sigma_xz=2*viscosity*epsilon[4];
		sigma_yz=2*viscosity*epsilon[5];

		/*Get normal vector to the bed */
		NormalBase(&bed_normal[0],&xyz_list_tria[0][0]);

		/*basalforce*/
		basalforce[0] += sigma_xx*bed_normal[0] + sigma_xy*bed_normal[1] + sigma_xz*bed_normal[2];
		basalforce[1] += sigma_xy*bed_normal[0] + sigma_yy*bed_normal[1] + sigma_yz*bed_normal[2];
		basalforce[2] += sigma_xz*bed_normal[0] + sigma_yz*bed_normal[1] + sigma_zz*bed_normal[2];

		GetTriaJacobianDeterminant(&Jdet2d, &xyz_list_tria[0][0],gauss);
		value+=sigma_zz*Jdet2d*gauss->weight;
		surface+=Jdet2d*gauss->weight;
	}
	value=value/surface;

	/*Add value to output*/
	//sigma_b->SetValue(id-1,value,INS_VAL);
}
/*}}}*/
void       Penta::ComputeDeviatoricStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  viscosity;
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,ezz,exy,exz,eyz];*/
	IssmDouble  tau_xx[NUMVERTICES];
	IssmDouble	tau_yy[NUMVERTICES];
	IssmDouble	tau_zz[NUMVERTICES];
	IssmDouble  tau_xy[NUMVERTICES];
	IssmDouble	tau_xz[NUMVERTICES];
	IssmDouble	tau_yz[NUMVERTICES];
	IssmDouble	tau_eff[NUMVERTICES];

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input=this->GetInput(VzEnum); _assert_(vz_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);

		/*Compute Stress*/
		tau_xx[iv]=2*viscosity*epsilon[0]; // tau = nu eps
		tau_yy[iv]=2*viscosity*epsilon[1];
		tau_zz[iv]=2*viscosity*epsilon[2];
		tau_xy[iv]=2*viscosity*epsilon[3];
		tau_xz[iv]=2*viscosity*epsilon[4];
		tau_yz[iv]=2*viscosity*epsilon[5];

		tau_eff[iv] = tau_xx[iv]*tau_xx[iv] + tau_yy[iv]*tau_yy[iv] + tau_zz[iv]*tau_zz[iv] +
		  2*tau_xy[iv]*tau_xy[iv] + 2*tau_xz[iv]*tau_xz[iv] + 2*tau_yz[iv]*tau_yz[iv];

		tau_eff[iv] = sqrt(tau_eff[iv]/2.);
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(DeviatoricStressxxEnum,&tau_xx[0],P1DGEnum);
	this->AddInput(DeviatoricStressxyEnum,&tau_xy[0],P1DGEnum);
	this->AddInput(DeviatoricStressxzEnum,&tau_xz[0],P1DGEnum);
	this->AddInput(DeviatoricStressyyEnum,&tau_yy[0],P1DGEnum);
	this->AddInput(DeviatoricStressyzEnum,&tau_yz[0],P1DGEnum);
	this->AddInput(DeviatoricStresszzEnum,&tau_zz[0],P1DGEnum);
	this->AddInput(DeviatoricStresseffectiveEnum,&tau_eff[0],P1DGEnum);
}
/*}}}*/
void       Penta::ComputeStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  pressure,viscosity;
	IssmDouble  epsilon[6]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  sigma_xx[NUMVERTICES];
	IssmDouble	sigma_yy[NUMVERTICES];
	IssmDouble	sigma_zz[NUMVERTICES];
	IssmDouble  sigma_xy[NUMVERTICES];
	IssmDouble	sigma_xz[NUMVERTICES];
	IssmDouble	sigma_yz[NUMVERTICES];

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* pressure_input=this->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=this->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);             _assert_(vy_input);
	Input* vz_input=this->GetInput(VzEnum);             _assert_(vz_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		pressure_input->GetInputValue(&pressure,&gauss);

		/*Compute Stress*/
		sigma_xx[iv]=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
		sigma_yy[iv]=2*viscosity*epsilon[1]-pressure;
		sigma_zz[iv]=2*viscosity*epsilon[2]-pressure;
		sigma_xy[iv]=2*viscosity*epsilon[3];
		sigma_xz[iv]=2*viscosity*epsilon[4];
		sigma_yz[iv]=2*viscosity*epsilon[5];
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(StressTensorxxEnum,&sigma_xx[0],P1DGEnum);
	this->AddInput(StressTensorxyEnum,&sigma_xy[0],P1DGEnum);
	this->AddInput(StressTensorxzEnum,&sigma_xz[0],P1DGEnum);
	this->AddInput(StressTensoryyEnum,&sigma_yy[0],P1DGEnum);
	this->AddInput(StressTensoryzEnum,&sigma_yz[0],P1DGEnum);
	this->AddInput(StressTensorzzEnum,&sigma_zz[0],P1DGEnum);
}
/*}}}*/
void       Penta::ComputeSigmaVM(){/*{{{*/

	if(!this->IsOnBase()) return;

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  lambda1,lambda2,ex,ey,vx,vy,vel;
	IssmDouble  B,n;
	IssmDouble  sigma_vm[NUMVERTICES];

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Depth average B for stress calculation*/
	this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input = this->GetInput(VxAverageEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyAverageEnum); _assert_(vy_input);
	Input* B_input  = this->GetInput(MaterialsRheologyBbarEnum);   _assert_(B_input);
	Input* n_input  = this->GetInput(MaterialsRheologyNEnum);   _assert_(n_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<3;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components and thickness*/
		B_input->GetInputValue(&B,&gauss);
		n_input->GetInputValue(&n,&gauss);
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=sqrt(vx*vx+vy*vy)+1.e-14;

		/*Compute strain rate and viscosity: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);

		/*Get Eigen values*/
		Matrix2x2Eigen(&lambda1,&lambda2,&ex,&ey,epsilon[0],epsilon[2],epsilon[1]);
		_assert_(!xIsNan<IssmDouble>(lambda1));
		_assert_(!xIsNan<IssmDouble>(lambda2));

		/*Process Eigen values (only account for extension)*/
		lambda1 = max(lambda1,0.);
		lambda2 = max(lambda2,0.);

		/*Calculate sigma_vm*/
		IssmDouble epse_2    = 1./2. *(lambda1*lambda1 + lambda2*lambda2);
		sigma_vm[iv] = sqrt(3.) * B * pow(epse_2,1./(2.*n));
	}

	/*Add input*/
	this->AddBasalInput(SigmaVMEnum,&sigma_vm[0],P1DGEnum);
	this->InputExtrude(SigmaVMEnum,-1);
}
/*}}}*/
void       Penta::Configure(Elements* elementsin, Loads* loadsin, Nodes* nodesin,Vertices* verticesin, Materials* materialsin, Parameters* parametersin,Inputs* inputsin){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	this->element_type=this->element_type_list[analysis_counter];

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective
	 * datasets, using internal ids and offsets hidden in hooks: */
	if (this->hnodes[analysis_counter]) this->hnodes[analysis_counter]->configure(nodesin);
	this->hvertices->configure(verticesin);
	this->hmaterial->configure(materialsin);
	this->hneighbors->configure(elementsin);

	/*Now, go pick up the objects inside the hooks: */
	if (this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices          = (Vertex**)this->hvertices->deliverp();
	this->material          = (Material*)this->hmaterial->delivers();
	this->verticalneighbors = (Penta**)this->hneighbors->deliverp();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
	this->inputs=inputsin;
}
/*}}}*/
void       Penta::ControlInputSetGradient(IssmDouble* gradient,int control_enum,int control_index,int offset,int M,int N,int interp){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         lidlist[NUMVERTICES];
	int         idlist[NUMVERTICES];

	if(control_enum==MaterialsRheologyBbarEnum) control_enum = MaterialsRheologyBEnum;
	if(control_enum==DamageDbarEnum)            control_enum = DamageDEnum;

	ElementInput* input=this->inputs->GetControlInputData(control_enum,"gradient");   _assert_(input);
	this->GetVerticesLidList(&lidlist[0]);
	GradientIndexing(&idlist[0],control_index);

	/*Get values on vertices*/
	if(input->ObjectEnum()==PentaInputEnum && input->GetInputInterpolationType()==P1Enum){
		_assert_(N==1);
		for(int i=0;i<NUMVERTICES;i++){
			values[i] = gradient[idlist[i]];
		}
		input->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
	}
	else if(input->ObjectEnum()==PentaInputEnum && input->GetInputInterpolationType()==P0Enum){
		_assert_(N==1);
		input->SetInput(P0Enum,this->lid,gradient[idlist[0]]);
	}
	else if(input->ObjectEnum()==TransientInputEnum){
		for(int n=0;n<N;n++){
			_error_("not implemented");
			//Input* new_input = new PentaInput(control_enum,gradient,P1Enum);
			//controlinput->SetInput(new_input,n);
			//controlinput->Configure(parameters);
		}
	}
	else _error_("Type not supported");

}
/*}}}*/
void       Penta::ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum,int control_interp){/*{{{*/

	int         sidlist[NUMVERTICES];
	int         lidlist[NUMVERTICES];
	int         connectivity[NUMVERTICES];
	IssmPDouble values[NUMVERTICES];
	IssmPDouble gradients[NUMVERTICES];
	IssmDouble  value,gradient;

	/*Get relevant inputs*/
	if(control_enum==MaterialsRheologyBbarEnum) control_enum = MaterialsRheologyBEnum;
	if(control_enum==DamageDbarEnum)            control_enum = DamageDEnum;
	ElementInput* control_value    = this->inputs->GetControlInputData(control_enum,"value");    _assert_(control_value);
	ElementInput* control_gradient = this->inputs->GetControlInputData(control_enum,"gradient"); _assert_(control_gradient);

	if(control_interp==P1Enum){
		_assert_(control_value->GetInputInterpolationType()==P1Enum);
		_assert_(control_gradient->GetInputInterpolationType()==P1Enum);

		this->GetVerticesConnectivityList(&connectivity[0]);
		this->GetVerticesSidList(&sidlist[0]);
		this->GetVerticesLidList(&lidlist[0]);

		control_value->Serve(NUMVERTICES,&lidlist[0]);
		control_gradient->Serve(NUMVERTICES,&lidlist[0]);

		GaussPenta gauss;
		for (int iv=0;iv<NUMVERTICES;iv++){
			gauss.GaussVertex(iv);

			control_value->GetInputValue(&value,&gauss);
			control_gradient->GetInputValue(&gradient,&gauss);

			values[iv]    = reCast<IssmPDouble>(value)/reCast<IssmPDouble>(connectivity[iv]);
			gradients[iv] = reCast<IssmPDouble>(gradient)/reCast<IssmPDouble>(connectivity[iv]);
		}

		vector_control->SetValues(NUMVERTICES,&sidlist[0],&values[0],ADD_VAL);
		vector_gradient->SetValues(NUMVERTICES,&sidlist[0],&gradients[0],ADD_VAL);
	}
	else if(control_interp==P0Enum){
		_assert_(control_value->GetInputInterpolationType()==P0Enum);
		_assert_(control_gradient->GetInputInterpolationType()==P0Enum);

		control_value->Serve(1,&this->lid);
		control_gradient->Serve(1,&this->lid);

		vector_control->SetValue(this->sid,reCast<IssmPDouble>(control_value->element_values[0]),ADD_VAL);
		vector_gradient->SetValue(this->sid,reCast<IssmPDouble>(control_gradient->element_values[0]),ADD_VAL);
	}
	else{
		_error_("not supported");
	}

}/*}}}*/
void       Penta::CreateDistanceInputFromSegmentlist(IssmDouble* distances,int distanceenum){/*{{{*/

	/*Get current field and vertex coordinates*/
	IssmDouble ls[NUMVERTICES],distance;
	Element::GetInputListOnVertices(&ls[0],distanceenum);

	/*Get distance from list of segments and reset ls*/
	for(int j=0;j<NUMVERTICES;j++){
		distance=distances[this->vertices[j]->Lid()];
		if(xIsNan<IssmDouble>(distance)) _error_("NaN found in vector");
		if(xIsInf<IssmDouble>(distance)) _error_("Inf found in vector");

		if(ls[j]>0){
			ls[j] = distance;
		}
		else{
			ls[j] = - distance;
		}
	}

	/*Update Levelset*/
	this->AddInput(distanceenum,&ls[0],P1Enum);
}
/*}}}*/
void       Penta::CreateInputTimeAverage(int transientinput_enum,int averagedinput_enum,IssmDouble start_time,IssmDouble end_time,int averaging_method){/*{{{*/
	_assert_(end_time>start_time);

	/*Get transient input time steps*/
	TransientInput* transient_input  = this->inputs->GetTransientInput(transientinput_enum);
	PentaInput* averaged_input = transient_input->GetPentaInput(start_time,end_time,averaging_method);
	Input* averaged_copy = averaged_input->copy();

	averaged_copy->ChangeEnum(averagedinput_enum);
	this->inputs->AddInput(averaged_copy);
}
/*}}}*/
void       Penta::ElementResponse(IssmDouble* presponse,int response_enum){/*{{{*/

	switch(response_enum){
		case MaterialsRheologyBbarEnum:
			*presponse=this->material->GetBbar(NULL);
			break;
		case DamageDbarEnum:
			*presponse=this->material->GetDbar(NULL);
			break;
		case VelEnum:
			{

				/*Get input:*/
				IssmDouble vel;
				Input* vel_input=this->GetInput(VelEnum); _assert_(vel_input);
				vel_input->GetInputAverage(&vel);

				/*Assign output pointers:*/
				*presponse=vel;
			}
			break;
		default:
			_error_("Response type " << EnumToStringx(response_enum) << " not supported yet!");
	}

}
/*}}}*/
void       Penta::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble xmin,ymin,zmin;
	IssmDouble xmax,ymax,zmax;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	xmin=xyz_list[0][0]; xmax=xyz_list[0][0];
	ymin=xyz_list[0][1]; ymax=xyz_list[0][1];
	zmin=xyz_list[0][2]; zmax=xyz_list[0][2];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<xmin) xmin=xyz_list[i][0];
		if(xyz_list[i][0]>xmax) xmax=xyz_list[i][0];
		if(xyz_list[i][1]<ymin) ymin=xyz_list[i][1];
		if(xyz_list[i][1]>ymax) ymax=xyz_list[i][1];
		if(xyz_list[i][2]<zmin) zmin=xyz_list[i][2];
		if(xyz_list[i][2]>zmax) zmax=xyz_list[i][2];
	}

	*hx=xmax-xmin;
	*hy=ymax-ymin;
	*hz=zmax-zmin;
}
/*}}}*/
int        Penta::FiniteElement(void){/*{{{*/
	return this->element_type;
}
/*}}}*/
IssmDouble Penta::FloatingArea(bool scaled){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,base_area,scalefactor,floatingarea;
	IssmDouble  xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || !IsOnBase())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	base_area= 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	floatingarea=(1-phi)*base_area;

	if(scaled==true){
		Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
		scalefactor_input->GetInputAverage(&scalefactor);
		floatingarea=floatingarea*scalefactor;
	}

	/*Clean up and return*/
	return floatingarea;
}
/*}}}*/
void       Penta::FSContactMigration(Vector<IssmDouble>* vertex_sigmann,Vector<IssmDouble>* vertex_waterpressure){/*{{{*/

	if(!IsOnBase()) return;

	int approximation;
	this->Element::GetInputValue(&approximation,ApproximationEnum);
	if(approximation==HOApproximationEnum || approximation==SSAApproximationEnum || approximation==SSAHOApproximationEnum || approximation==HOFSApproximationEnum){
		_error_("Cannot compute contact condition for non FS elements");
	}

	/*Intermediaries*/
	IssmDouble  bed_normal[3],base[NUMVERTICES],bed[NUMVERTICES],surface[NUMVERTICES],phi[NUMVERTICES];
	IssmDouble  water_pressure[NUMVERTICES],pressureice[NUMVERTICES],pressure[NUMVERTICES];
	IssmDouble  sigmaxx[NUMVERTICES],sigmayy[NUMVERTICES],sigmazz[NUMVERTICES],sigmaxy[NUMVERTICES];
	IssmDouble  sigmayz[NUMVERTICES],sigmaxz[NUMVERTICES],sigma_nn[NUMVERTICES];
	IssmDouble  viscosity,epsilon[NUMVERTICES];
	Element::GetInputListOnVertices(&base[0],BaseEnum);
	Element::GetInputListOnVertices(&bed[0],BedEnum);
	Element::GetInputListOnVertices(&surface[0],SurfaceEnum);
	Element::GetInputListOnVertices(&pressure[0],PressureEnum);
	Element::GetInputListOnVertices(&phi[0],MaskOceanLevelsetEnum);
	IssmDouble rho_ice   = FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble gravity   = FindParam(ConstantsGEnum);

	/* Get node coordinates and dof list: */
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = this->GetInput(VzEnum); _assert_(vz_input);

	/*1. Recover stresses at the base*/
	GaussPenta gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		this->material->ViscosityFS(&viscosity,3,&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		/*FIXME: this is for Hongju only*/
		//pressureice[iv]=gravity*rho_ice*(surface[iv]-base[iv]);
		//if (pressure[iv]/pressureice[iv]>1)	pressure[iv]=pressureice[iv];

		/*Compute Stress*/
		sigmaxx[iv]=2*viscosity*epsilon[0]-pressure[iv]; // sigma = nu eps - pressure
		sigmayy[iv]=2*viscosity*epsilon[1]-pressure[iv];
		sigmazz[iv]=2*viscosity*epsilon[2]-pressure[iv];
		sigmaxy[iv]=2*viscosity*epsilon[3];
		sigmaxz[iv]=2*viscosity*epsilon[4];
		sigmayz[iv]=2*viscosity*epsilon[5];
	}

	/*2. compute contact condition*/
	for(int i=0;i<NUMVERTICES;i++){

		/*If was grounded*/
		if (phi[i]>=0.){
			NormalBase(&bed_normal[0],&xyz_list[0][0]);
			sigma_nn[i]=-1*(sigmaxx[i]*bed_normal[0]*bed_normal[0] + sigmayy[i]*bed_normal[1]*bed_normal[1] + sigmazz[i]*bed_normal[2]*bed_normal[2]+2*sigmaxy[i]*bed_normal[0]*bed_normal[1]+2*sigmaxz[i]*bed_normal[0]*bed_normal[2]+2*sigmayz[i]*bed_normal[1]*bed_normal[2]);
			water_pressure[i]=-gravity*rho_water*base[i];
			vertex_sigmann->SetValue(vertices[i]->Pid(),sigma_nn[i],ADD_VAL);
			vertex_waterpressure->SetValue(vertices[i]->Pid(),water_pressure[i],ADD_VAL);
		}

		/*If was floating*/
		else{
			/*Tricky part:
			 * 1. if base is now touching, we put 1 for sigma_nn and leave water pressure at 0 so that the rest of the module will reground this vertex
			 * 2. if base is still above bed, water pressure is set as 1, sigma_nn is left as 0, so the GL module will keep it afloat*/
			if(base[i]<bed[i]) vertex_sigmann->SetValue(vertices[i]->Pid(),+1.,ADD_VAL);
			else vertex_waterpressure->SetValue(vertices[i]->Pid(),+1.,ADD_VAL);
		}
	}
}
/*}}}*/
void       Penta::GetAreaCoordinates(IssmDouble* area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints){/*{{{*/
	/*Computeportion of the element that is grounded*/

	int         i,j,k;
	IssmDouble  area_init,area_portion;
	IssmDouble  xyz_bis[3][3];

	area_init=fabs(xyz_list[1*3+0]*xyz_list[2*3+1] - xyz_list[1*3+1]*xyz_list[2*3+0] + xyz_list[0*3+0]*xyz_list[1*3+1] - xyz_list[0*3+1]*xyz_list[1*3+0] + xyz_list[2*3+0]*xyz_list[0*3+1] - xyz_list[2*3+1]*xyz_list[0*3+0])/2.;

	/*Initialize xyz_list with original xyz_list of triangle coordinates*/
	for(j=0;j<3;j++){
		for(k=0;k<3;k++){
			xyz_bis[j][k]=xyz_list[j*3+k];
		}
	}
	for(i=0;i<numpoints;i++){
		for(j=0;j<3;j++){
			for(k=0;k<3;k++){
				/*Change appropriate line*/
				xyz_bis[j][k]=xyz_zero[i*3+k];
			}

			/*Compute area fraction*/
			area_portion=fabs(xyz_bis[1][0]*xyz_bis[2][1] - xyz_bis[1][1]*xyz_bis[2][0] + xyz_bis[0][0]*xyz_bis[1][1] - xyz_bis[0][1]*xyz_bis[1][0] + xyz_bis[2][0]*xyz_bis[0][1] - xyz_bis[2][1]*xyz_bis[0][0])/2.;
			*(area_coordinates+3*i+j)=area_portion/area_init;

			/*Reinitialize xyz_list*/
			for(k=0;k<3;k++){
				/*Reinitialize xyz_list with original coordinates*/
				xyz_bis[j][k]=xyz_list[j*3+k];
			}
		}
	}
}
/*}}}*/
Element*   Penta::GetBasalElement(void){/*{{{*/

	/*Output*/
	Element* element=this->GetBasalPenta();
	return element;
}
/*}}}*/
Penta*     Penta::GetBasalPenta(void){/*{{{*/

	/*Output*/
	Penta* penta=NULL;

	/*Go through all pentas till the bed is reached*/
	penta=this;
	for(;;){
		/*Stop if we have reached the surface, else, take lower penta*/
		if (penta->IsOnBase()) break;

		/* get lower Penta*/
		penta=penta->GetLowerPenta();
		_assert_(penta->Id()!=this->id);
	}

	/*return output*/
	return penta;
}
/*}}}*/
int        Penta::GetElementType(){/*{{{*/

	/*return PentaRef field*/
	return this->element_type;
}
/*}}}*/
void       Penta::GetFractionGeometry2D(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl){/*{{{*/

  /*Compute portion of element that is grounded based on levelset at the 3 lower vertices of the Penta element*/
   bool               trapezeisnegative=true;
   int                point;
   const IssmPDouble  epsilon= 1.e-15;
   IssmDouble         f1,f2,phi;

   /*Weights*/
   IssmDouble loadweights_g[NUMVERTICES];
   IssmDouble total_weight = 0;

   _assert_(!xIsNan<IssmDouble>(gl[0]));
   _assert_(!xIsNan<IssmDouble>(gl[1]));
   _assert_(!xIsNan<IssmDouble>(gl[2]));

   /*Be sure that values are not zero*/
   if(gl[0]==0.) gl[0] = gl[0]+epsilon;
   if(gl[1]==0.) gl[1] = gl[1]+epsilon;
   if(gl[2]==0.) gl[2] = gl[2]+epsilon;

   /*Check that not all nodes are positive or negative: */
   if(gl[0]>0 && gl[1]>0 && gl[2]>0){
      point = 0;
      f1    = 1.;
      f2    = 1.;
   }
   else if(gl[0]<0 && gl[1]<0 && gl[2]<0){
      point = 0;
      f1    = 0.;
      f2    = 0.;
   }
	else{
		if(gl[0]*gl[1]*gl[2]<0) trapezeisnegative = false;

		/*Find the similar nodes*/
		if(gl[0]*gl[1]>0){ 
			point = 2;
			f1    = gl[2]/(gl[2]-gl[0]);
			f2    = gl[2]/(gl[2]-gl[1]);
		}
		else if(gl[1]*gl[2]>0){ 
			point = 0;
			f1    = gl[0]/(gl[0]-gl[1]);
			f2    = gl[0]/(gl[0]-gl[2]);
		}
		else if(gl[0]*gl[2]>0){ 
			point = 1;
			f1    = gl[1]/(gl[1]-gl[2]);
			f2    = gl[1]/(gl[1]-gl[0]);
		}
		else _error_("case not possible");
	}
	if(trapezeisnegative) phi = 1.-f1*f2;
	else                  phi = f1*f2;

	/*Compute weights*/
	Gauss* gauss = this->NewGauss(point,f1,f2,1-trapezeisnegative,2);

	total_weight = 0.0;
	for(int i=0;i<NUMVERTICES2D;i++)weights[i] = 0;
	while(gauss->next()){
		GetNodalFunctions(&loadweights_g[0],gauss,P1Enum);
		for(int i=0;i<NUMVERTICES2D;i++)weights[i] += loadweights_g[i]*gauss->weight;
		total_weight += gauss->weight;
	}
	delete gauss;

	/*Normalizing to phi such that weights provide coefficients for integration over subelement (for averaging:phi*weights)*/
   if(total_weight>0.) for(int i=0;i<NUMVERTICES2D;i++) weights[i] = weights[i]*phi/total_weight;
	else for(int i=0;i<NUMVERTICES2D;i++) weights[i] = 0.0;

	/*Assign output pointers*/
	*pphi               = phi;
	*ppoint1            = point;
	*pfraction1         = f1;
	*pfraction2         = f2;
	*ptrapezeisnegative = trapezeisnegative;
}
/*}}}*/
void       Penta::GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlyfloating, int distance_enum, IssmDouble intrusion_distance){/*{{{*/
	/*Compute portion of the element that is grounded*/
	bool               floating=true;
	int                point, melt_style;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         gl[NUMVERTICES];
	IssmDouble         f1,f2;

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],distance_enum);

	/*Determine where to apply sub-element melt using intrusion_distance*/
	for(int i=0; i<NUMVERTICES; i++){
		gl[i] -= intrusion_distance;
	}

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are grounded or floating*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All grounded
		point=0;
		f1=1.;
		f2=1.;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All floating
		point=0;
		f1=0.;
		f2=0.;
	}
	else{
		if(gl[0]*gl[1]*gl[2]<0) floating=false;

		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			point=2;
			f1=gl[2]/(gl[2]-gl[0]);
			f2=gl[2]/(gl[2]-gl[1]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			point=0;
			f1=gl[0]/(gl[0]-gl[1]);
			f2=gl[0]/(gl[0]-gl[2]);
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			point=1;
			f1=gl[1]/(gl[1]-gl[2]);
			f2=gl[1]/(gl[1]-gl[0]);
		}
		else _error_("case not possible");
	}
	*point1=point;
	*fraction1=f1;
	*fraction2=f2;
	*mainlyfloating=floating;
}
/*}}}*/
IssmDouble Penta::GetGroundedPortion(IssmDouble* xyz_list){/*{{{*/
	/*Computeportion of the element that is grounded*/

	bool               mainlyfloating = true;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         phi,s1,s2;
	IssmDouble         gl[NUMVERTICES];

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are grounded or floating*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All grounded
		phi=1;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All floating
		phi=0;
	}
	else{
		/*Figure out if two nodes are floating or grounded*/
		if(gl[0]*gl[1]*gl[2]>0) mainlyfloating=false;

		if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2
			s1=gl[2]/(gl[2]-gl[1]);
			s2=gl[2]/(gl[2]-gl[0]);
		}
		else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2
			s1=gl[0]/(gl[0]-gl[1]);
			s2=gl[0]/(gl[0]-gl[2]);
		}
		else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2
			s1=gl[1]/(gl[1]-gl[0]);
			s2=gl[1]/(gl[1]-gl[2]);
		}
		else _error_("case not possible");
		if(mainlyfloating){
			phi = (1-s1*s2);
		}
		else{
			phi = s1*s2;
		}
	}

	_assert_(phi<=1. && phi>=0.);
	return phi;
}
/*}}}*/
IssmDouble Penta::GetIcefrontArea(){/*{{{*/

	/*We need to be on base and cross the levelset*/
	if(!IsZeroLevelset(MaskIceLevelsetEnum)) return 0;
	if(!this->IsOnBase()) return 0;

	/*Spawn Tria element from the base of the Penta: */
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	IssmDouble frontarea = tria->GetIcefrontArea();
	delete tria->material; delete tria;

	return frontarea;
}/*}}}*/
void       Penta::GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum){/*{{{*/

	/* Intermediaries */
	const int dim=3;
	int i, dir,nrfrontnodes;
	IssmDouble  levelset[NUMVERTICES];

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&levelset[0],levelsetenum);

	int* indicesfront = xNew<int>(NUMVERTICES);
	/* Get basal nodes where there is no ice */
	nrfrontnodes=0;
	for(i=0;i<NUMVERTICES2D;i++){
		if(levelset[i]>=0.){
			indicesfront[nrfrontnodes]=i;
			nrfrontnodes++;
		}
	}
	_assert_(nrfrontnodes==2);

	/* arrange order of basal frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES2D+indicesfront[0]-indicesfront[1])%NUMVERTICES2D!=NUMVERTICES2D-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}

	IssmDouble* xyz_front = xNew<IssmDouble>(2*dim*nrfrontnodes);
	/* Return basal and top front nodes */
	for(i=0;i<nrfrontnodes;i++){
		for(dir=0;dir<dim;dir++){
			int ind1=i*dim+dir, ind2=(2*nrfrontnodes-1-i)*dim+dir; // vertex structure front segment: base0, base1, top1, top0
			xyz_front[ind1]=xyz_list[dim*indicesfront[i]+dir];
			xyz_front[ind2]=xyz_list[dim*(indicesfront[i]+NUMVERTICES2D)+dir];
		}
	}

	*pxyz_front=xyz_front;

	xDelete<int>(indicesfront);
}/*}}}*/
Input*    Penta::GetInput(int inputenum){/*{{{*/

	/*Get Input from dataset*/
	PentaInput* input = this->inputs->GetPentaInput(inputenum);
	if(!input) return input;

	/*Intermediaries*/
	int numindices;
	int indices[30]; /*Max numnodes*/

	/*Check interpolation*/
	int interpolation = input->GetInterpolation();
	if(interpolation==P1Enum){
		numindices = 6;
		for(int i=0;i<6;i++) indices[i] = vertices[i]->lid;
		input->Serve(numindices,&indices[0]);
	}
	else{
		input->Serve(this->lid,this->GetNumberOfNodes(interpolation));
	}

	/*Tell input it is NOT collapsed*/
	//input->SetServeCollapsed(0); FIXME: not needed?

	/*Return*/
	return input;
}/*}}}*/
Input*    Penta::GetInput(int inputenum,IssmDouble time){/*{{{*/

	/*Get Input from dataset*/
	PentaInput* input = this->inputs->GetPentaInput(inputenum,time);
	if(!input) return input;

	/*Intermediaries*/
	int numindices;
	int indices[30]; /*Max numnodes*/

	/*Check interpolation*/
	int interpolation = input->GetInterpolation();
	if(interpolation==P1Enum){
		numindices = 6;
		for(int i=0;i<6;i++) indices[i] = vertices[i]->lid;
		input->Serve(numindices,&indices[0]);
	}
	else{
		input->Serve(this->lid,this->GetNumberOfNodes(interpolation));
	}

	/*Tell input it is NOT collapsed*/
	//input->SetServeCollapsed(0); FIXME: not needed?

	/*Return*/
	return input;
}/*}}}*/
void       Penta::GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/* Start looping on the number of vertices: */
	if(input){
		GaussPenta gauss;
		for(int iv=0;iv<NUMVERTICES;iv++){
			gauss.GaussVertex(iv);
			input->GetInputValue(&pvalue[iv],&gauss);
		}
	}
	else{
		for(int iv=0;iv<NUMVERTICES;iv++) pvalue[iv] = default_value;
	}
}
/*}}}*/
void       Penta::GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/*What type of finite element are we dealing with?*/
	int fe       = this->FiniteElement();
	int numnodes = this->GetNumberOfNodes();

	/* Start looping on the number of vertices: */
	if(input){
		GaussPenta gauss;
		for(int iv=0;iv<numnodes;iv++){
			gauss.GaussNode(fe,iv);
			input->GetInputValue(&pvalue[iv],&gauss);
		}
	}
	else{
		for(int iv=0;iv<numnodes;iv++) pvalue[iv] = default_value;
	}
}
/*}}}*/
DatasetInput* Penta::GetDatasetInput(int inputenum){/*{{{*/

	DatasetInput* datasetinput = this->inputs->GetDatasetInput(inputenum);
	if(!datasetinput) return NULL;

	for(int i=0;i<datasetinput->GetNumIds();i++){

		PentaInput* input = datasetinput->GetPentaInputByOffset(i); _assert_(input);

		/*Intermediaries*/
		int numindices;
		int indices[30]; /*Max numnodes*/

		/*Check interpolation*/
		int interpolation = input->GetInterpolation();
		if(interpolation==P1Enum){
			numindices = 6;
			for(int i=0;i<6;i++) indices[i] = vertices[i]->lid;
			input->Serve(numindices,&indices[0]);
		}
		else{
			input->Serve(this->lid,this->GetNumberOfNodes(interpolation));
		}

	}

	return datasetinput;
}/*}}}*/
void       Penta::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=this->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	int index = this->GetNodeIndex(node);

	GaussPenta gauss;
	gauss.GaussNode(this->element_type,index);
	input->GetInputValue(pvalue,&gauss);
}
/*}}}*/
void       Penta::GetInputValue(IssmDouble* pvalue,Vertex* vertex,int enumtype){/*{{{*/

	Input* input=this->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	int index = this->GetVertexIndex(vertex);

	GaussPenta gauss;
	gauss.GaussVertex(index);
	input->GetInputValue(pvalue,&gauss);
}
/*}}}*/
Penta*     Penta::GetLowerPenta(void){/*{{{*/

	Penta* lower_penta=NULL;

	lower_penta=(Penta*)verticalneighbors[0]; //first one (0) under, second one (1) above

	return lower_penta;
}
/*}}}*/
void       Penta::GetLevelsetIntersectionBase(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level){/*{{{*/

	/* GetLevelsetIntersection computes:
	 * 1. indices of element, sorted in [iceverts, noiceverts] in counterclockwise fashion,
	 * 2. fraction of intersected triangle edges intersected by levelset, lying below level*/

	/*Intermediaries*/
	int i, numiceverts, numnoiceverts;
	int ind0, ind1, lastindex;
	int indices_ice[NUMVERTICES2D],indices_noice[NUMVERTICES2D];
	IssmDouble lsf[NUMVERTICES];
	int* indices = xNew<int>(NUMVERTICES2D);

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&lsf[0],levelset_enum);

	/* Determine distribution of ice over element.
	 * Exploit: ice/no-ice parts are connected, so find starting vertex of segment*/
	lastindex=0;
	for(i=0;i<NUMVERTICES2D;i++){ // go backwards along vertices, and check for sign change
		ind0=(NUMVERTICES2D-i)%NUMVERTICES2D;
		ind1=(ind0-1+NUMVERTICES2D)%NUMVERTICES2D;
		if((lsf[ind0]-level)*(lsf[ind1]-level)<=0.){ // levelset has been crossed, find last index belonging to segment
			if(lsf[ind1]==level) //if levelset intersects 2nd vertex, choose this vertex as last
				lastindex=ind1;
			else
				lastindex=ind0;
			break;
		}
	}

	numiceverts=0;
	numnoiceverts=0;
	for(i=0;i<NUMVERTICES2D;i++){
		ind0=(lastindex+i)%NUMVERTICES2D;
		if(lsf[i]<=level){
			indices_ice[numiceverts]=i;
			numiceverts++;
		}
		else{
			indices_noice[numnoiceverts]=i;
			numnoiceverts++;
		}
	}
	//merge indices
	for(i=0;i<numiceverts;i++){indices[i]=indices_ice[i];}
	for(i=0;i<numnoiceverts;i++){indices[numiceverts+i]=indices_noice[i];}

	switch (numiceverts){
		case 0: // no vertex has ice: element is ice free, no intersection
			for(i=0;i<2;i++)
				fraction[i]=0.;
			break;
		case 1: // one vertex has ice:
			for(i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[0]])/(lsf[indices[numiceverts+i]]-lsf[indices[0]]);
			}
			break;
		case 2: // two vertices have ice: fraction is computed from first ice vertex to last in CCW fashion
			for(i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[i]])/(lsf[indices[numiceverts]]-lsf[indices[i]]);
			}
			break;
		case NUMVERTICES2D: // all vertices have ice: return triangle area
			for(i=0;i<2;i++)
				fraction[i]=1.;
			break;
		default:
			_error_("Wrong number of ice vertices in Penta::GetLevelsetIntersectionBase!");
			break;
	}

	*pindices=indices;
	*pnumiceverts=numiceverts;
}
/*}}}*/
int        Penta::GetVertexIndex(Vertex* vertex){/*{{{*/
	_assert_(vertices);
	for(int i=0;i<NUMVERTICES;i++){
		if(vertex==vertices[i])
		 return i;
	}
	_error_("Vertex provided not found among element vertices");
}
/*}}}*/
int        Penta::GetNumberOfNodes(void){/*{{{*/
	return this->NumberofNodes(this->element_type);
}
/*}}}*/
int        Penta::GetNumberOfNodes(int enum_type){/*{{{*/
	return this->NumberofNodes(enum_type);
}
/*}}}*/
int        Penta::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
Penta*     Penta::GetUpperPenta(void){/*{{{*/

	Penta* upper_penta=NULL;

	upper_penta=(Penta*)verticalneighbors[1]; //first one (0) under, second one (1) above

	return upper_penta;
}
/*}}}*/
void       Penta::GetVectorFromControlInputs(Vector<IssmDouble>* vector,int control_enum,int control_index,int N,const char* data,int offset){/*{{{*/

	/*Get input*/
	if(control_enum==MaterialsRheologyBbarEnum) control_enum=MaterialsRheologyBEnum;
	ElementInput* input=this->inputs->GetControlInputData(control_enum,data);   _assert_(input);

	/*Lid list once for all*/
	int lidlist[NUMVERTICES];
	for(int i=0;i<NUMVERTICES;i++) lidlist[i] = vertices[i]->lid;

	/*Check what input we are dealing with*/
	switch(input->ObjectEnum()){
		case PentaInputEnum:
			  {
				IssmDouble values[NUMVERTICES];
				int        idlist[NUMVERTICES];

				PentaInput* triainput = xDynamicCast<PentaInput*>(input);

				/*Create list of indices and values for global vector*/
				GradientIndexing(&idlist[0],control_index);

				if(triainput->GetInputInterpolationType()==P1Enum){
					input->Serve(NUMVERTICES,&lidlist[0]);
					for(int i=0;i<NUMVERTICES;i++) values[i] = triainput->element_values[i];
					vector->SetValues(NUMVERTICES,idlist,values,INS_VAL);
				}
				else if(triainput->GetInputInterpolationType()==P0Enum){
					input->Serve(1,&this->lid);
					vector->SetValue(idlist[0],triainput->element_values[0],INS_VAL);
				}
				else{
					_error_("not supported yet");
				}
				break;
			  }

		case TransientInputEnum:
				{
					TransientInput* transientinput = xDynamicCast<TransientInput*>(input);
					int  N = transientinput->numtimesteps;
					int* M = NULL;
					parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
					int* idlist = xNew<int>(NUMVERTICES*N);
					IssmDouble* values = xNew<IssmDouble>(NUMVERTICES*N);
					for(int t=0;t<transientinput->numtimesteps;t++) {
						IssmDouble time = transientinput->GetTimeByOffset(t);
						_error_("not implemented, SEE TRIA!");
						//PentaInput* timeinput = xDynamicCast<PentaInput*>(transientinput->GetTimeInput(time));
						//if(timeinput->interpolation_type!=P1Enum) _error_("not supported yet");
						//input->Serve(NUMVERTICES,&lidlist[0]);
						///*Create list of indices and values for global vector*/
						//for(int i=0;i<NUMVERTICES;i++){
						//		idlist[N*i+t] = offset + this->vertices[i]->Sid()+t*M[control_index];
						//		values[N*i+t] = timeinput->values[i];
						//}
					}
					vector->SetValues(NUMVERTICES*transientinput->numtimesteps,idlist,values,INS_VAL);
					xDelete<int>(M);
					xDelete<int>(idlist);
					xDelete<IssmDouble>(values);
					break;
				}
		default: _error_("input "<<EnumToStringx(input->ObjectEnum())<<" not supported yet");
	}
}/*}}}*/
void       Penta::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES2D*3);
	::GetVerticesCoordinates(xyz_list,this->vertices,NUMVERTICES2D);

	/*Assign output pointer*/
	*pxyz_list = xyz_list;

}/*}}}*/
void       Penta::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	IssmDouble* xyz_list = xNew<IssmDouble>(NUMVERTICES2D*3);
	::GetVerticesCoordinates(xyz_list,&this->vertices[3],NUMVERTICES2D);

	/*Assign output pointer*/
	*pxyz_list = xyz_list;

}/*}}}*/
IssmDouble Penta::GroundedArea(bool scaled){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,base_area,scalefactor,groundedarea;
	IssmDouble  xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || !IsOnBase())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	base_area= 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	groundedarea=phi*base_area;

	if(scaled==true){
		Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
		scalefactor_input->GetInputAverage(&scalefactor);
		groundedarea=groundedarea*scalefactor;
	}
	/*Clean up and return*/
	return groundedarea;
}
/*}}}*/
IssmDouble Penta::GroundinglineMassFlux(bool scaled){/*{{{*/

	/*Make sure there is a grounding line here*/
	if(!IsOnBase()) return 0;
	if(!IsIceInElement()) return 0;
	if(!IsZeroLevelset(MaskOceanLevelsetEnum)) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	int pt1 = 0;
	int pt2 = 1;
	if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

		/*Portion of the segments*/
		s1=gl[2]/(gl[2]-gl[1]);
		s2=gl[2]/(gl[2]-gl[0]);
		if(gl[2]<0.){
			pt1 = 1; pt2 = 0;
		}
		xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
		xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
		xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
		xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
		xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
		xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
	}
	else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

		/*Portion of the segments*/
		s1=gl[0]/(gl[0]-gl[1]);
		s2=gl[0]/(gl[0]-gl[2]);
		if(gl[0]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
		xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
		xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
		xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
		xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
		xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
	}
	else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

		/*Portion of the segments*/
		s1=gl[1]/(gl[1]-gl[0]);
		s2=gl[1]/(gl[1]-gl[2]);
		if(gl[1]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
		xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
		xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
		xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
		xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
		xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
	}
	else{
		_error_("case not possible");
	}

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input *thickness_input = this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input *vx_input        = this->GetInput(VxAverageEnum); _assert_(vx_input);
	Input *vy_input        = this->GetInput(VyAverageEnum); _assert_(vy_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(vx*normal[0] + vy*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Penta::IcefrontMassFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum) || !IsOnBase()) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	int pt1 = 0;
	int pt2 = 1;
	if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

		/*Portion of the segments*/
		s1=gl[2]/(gl[2]-gl[1]);
		s2=gl[2]/(gl[2]-gl[0]);
		if(gl[2]<0.){
			pt1 = 1; pt2 = 0;
		}
		xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
		xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
		xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
		xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
		xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
		xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
	}
	else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

		/*Portion of the segments*/
		s1=gl[0]/(gl[0]-gl[1]);
		s2=gl[0]/(gl[0]-gl[2]);
		if(gl[0]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
		xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
		xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
		xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
		xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
		xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
	}
	else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

		/*Portion of the segments*/
		s1=gl[1]/(gl[1]-gl[0]);
		s2=gl[1]/(gl[1]-gl[2]);
		if(gl[1]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
		xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
		xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
		xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
		xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
		xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
	}
	else{
		_error_("case not possible");
	}

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
	vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(vx*normal[0] + vy*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Penta::IceVolume(bool scaled){/*{{{*/

	/*The volume of a troncated prism is base * 1/3 sum(length of edges)*/
	IssmDouble base,height,scalefactor;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble lsf[NUMVERTICES];

	if(!IsIceInElement()) return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
	/*Deal with partially ice-covered elements*/
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		bool        istrapneg;
		int         point;
		IssmDouble  f1,f2,phi;
		IssmDouble  heights[NUMVERTICES2D];
		IssmDouble  weights[NUMVERTICES2D];
		IssmDouble  lsf2d[NUMVERTICES2D];
		for(int i=0;i<NUMVERTICES2D;i++){
			heights[i] = xyz_list[i+NUMVERTICES2D][2]-xyz_list[i][2];
			lsf2d[i]   = lsf[i];
		}
		GetFractionGeometry2D(&weights[0],&phi,&point,&f1,&f2,&istrapneg,lsf2d);

		IssmDouble basetot;
		height = 0.0;
		for(int i=0;i<NUMVERTICES2D;i++) height += weights[i]/phi*heights[i];
		basetot = 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));
		base    = basetot*phi;	

		/*Account for scaling factor averaged over subelement 2D area*/
		if(scaled==true){
			IssmDouble scalefactor_vertices[NUMVERTICES];
			Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
			/*Compute loop only over lower vertices: i<NUMVERTICES2D*/
			scalefactor = 0.0;
			for(int i=0;i<NUMVERTICES2D;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
			base = base*scalefactor;
		}
	}

	else{ 
		/*First calculate the area of the base (cross section triangle)
		 * http://en.wikipedia.org/wiki/Pentangle
		 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
		base = 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

		if(scaled==true){ //scale for area projection correction
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);
			base=base*scalefactor;
		}

		/*Now get the average height*/
		height = 1./3.*((xyz_list[3][2]-xyz_list[0][2])+(xyz_list[4][2]-xyz_list[1][2])+(xyz_list[5][2]-xyz_list[2][2]));
	}

	/*Return: */
	return base*height;
}
/*}}}*/
IssmDouble Penta::IceVolumeAboveFloatation(bool scaled){/*{{{*/

	/*Volume above floatation: H + rho_water/rho_ice*bathymetry for nodes on the bed*/
	IssmDouble rho_ice,rho_water;
	IssmDouble base,bed,surface,bathymetry,scalefactor;
	IssmDouble xyz_list[NUMVERTICES][3];

	if(!IsIceInElement() || IsAllFloating() || !IsOnBase())return 0;

	rho_ice=FindParam(MaterialsRhoIceEnum);
	rho_water=FindParam(MaterialsRhoSeawaterEnum);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Pentangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2.*fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));
	if(scaled==true){
		Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
		scalefactor_input->GetInputAverage(&scalefactor);
		base=base*scalefactor;
	}

	/*Now get the average height above floatation*/
	Input* surface_input    = this->GetInput(SurfaceEnum);    _assert_(surface_input);
	Input* base_input        = this->GetInput(BaseEnum);        _assert_(base_input);
	Input* bed_input = this->GetInput(BedEnum); _assert_(bed_input);
	if(!bed_input) _error_("Could not find bed");
	surface_input->GetInputAverage(&surface);
	base_input->GetInputAverage(&bed);
	bed_input->GetInputAverage(&bathymetry);

	/*Return: */
	return base*(surface - bed + min( rho_water/rho_ice * bathymetry, 0.) );
}
/*}}}*/
void       Penta::InputDepthAverageAtBase(int original_enum,int average_enum){/*{{{*/

	IssmDouble  Jdet,value;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_line[2][3];
	IssmDouble  total[NUMVERTICES]       = {0.};
	int         lidlist[NUMVERTICES];
	IssmDouble  intz[NUMVERTICES]        = {0.};
	Input     *original_input           = NULL;
	Input     *depth_averaged_input     = NULL;

	/*Are we on the base? If not, return*/
	if(!IsOnBase()) return;

	/*Now follow all the upper element from the base to the surface to integrate the input*/
	Penta* penta = this;
	int    step  = 0;
	Gauss* gauss[3];
	for(int iv=0;iv<3;iv++) gauss[iv] = penta->NewGaussLine(iv,iv+3,3);

	for(;;){

		/*Step1: Get original input (to be depth-avegaged): */
		original_input=penta->GetInput(original_enum);
		if(!original_input) _error_("could not find input with enum " << EnumToStringx(original_enum));

		/*Step2: Create element thickness input*/
		::GetVerticesCoordinates(&xyz_list[0][0],penta->vertices,NUMVERTICES);
		for(int iv=0;iv<3;iv++){
			/*Get segment length*/
			for(int i=0;i<3;i++){
				xyz_list_line[0][i]=xyz_list[iv][i];
				xyz_list_line[1][i]=xyz_list[iv+3][i];
			}
			/*Integrate over edge*/
			gauss[iv]->Reset();
			while(gauss[iv]->next()){
				penta->JacobianDeterminantLine(&Jdet,&xyz_list_line[0][0],gauss[iv]);
				original_input->GetInputValue(&value,gauss[iv]);
				total[iv] += value*Jdet*gauss[iv]->weight;
				intz[iv]  += Jdet*gauss[iv]->weight;
			}
		}

		/*Stop if we have reached the surface, else, take upper penta*/
		if(penta->IsOnSurface()) break;

		/* get upper Penta*/
		penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
		step++;
	}
	for(int iv=0;iv<3;iv++) delete gauss[iv];

	/*Now we only need to divide the depth integrated input by the total thickness!*/
	for(int iv=0;iv<3;iv++){
		total[iv  ] = total[iv]/intz[iv];
		total[iv+3] = total[iv];
	}
	GetVerticesLidList(&lidlist[0]);
	switch(original_input->ObjectEnum()){
		case PentaInputEnum:
		case ControlInputEnum:
			this->inputs->SetPentaInput(average_enum,P1Enum,NUMVERTICES,lidlist,&total[0]);
			break;
		default:
			_error_("Interpolation " << EnumToStringx(original_input->ObjectEnum()) << " not supported yet");
	}
}
/*}}}*/
void       Penta::DatasetInputExtrude(int enum_type,int start){/*{{{*/

	_assert_(start==-1 || start==+1);
	_assert_(this->inputs);

	/*Are we on the the boundary we want to be?*/
	if(start==-1 && !IsOnBase())    return;
	if(start==+1 && !IsOnSurface()) return;

	/*Get original input*/
	DatasetInput* dinput = this->inputs->GetDatasetInput(enum_type);

	int lidlist[NUMVERTICES];
	this->GetVerticesLidList(&lidlist[0]);

	for(int id=0;id<dinput->GetNumIds();id++){

		PentaInput* pentainput = dinput->GetPentaInputByOffset(id);
		pentainput->Serve(NUMVERTICES,&lidlist[0]);

		if(pentainput->GetInterpolation()==P1Enum){

			/*Extrude values first*/
			IssmDouble extrudedvalues[NUMVERTICES];
			this->GetInputListOnVertices(&extrudedvalues[0],pentainput,0.);

			if(start==-1){
				for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i+NUMVERTICES2D]=extrudedvalues[i];
			}
			else{
				for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i]=extrudedvalues[i+NUMVERTICES2D];
			}

			/*Propagate to other Pentas*/
			Penta* penta=this;
			for(;;){

				/*Add input of the basal element to penta->inputs*/
				int vertexlids[NUMVERTICES];
				penta->GetVerticesLidList(&vertexlids[0]);
				pentainput->SetInput(P1Enum,NUMVERTICES,&vertexlids[0],&extrudedvalues[0]);

				/*Stop if we have reached the surface/base*/
				if(start==-1 && penta->IsOnSurface()) break;
				if(start==+1 && penta->IsOnBase())    break;

				/*get upper/lower Penta*/
				if(start==-1) penta=penta->GetUpperPenta();
				else          penta=penta->GetLowerPenta();
				_assert_(penta->Id()!=this->id);
			}
		}
		else{
			_error_("not implemented yet");
		}
	}
}
/*}}}*/
void       Penta::ControlInputExtrude(int enum_type,int start){/*{{{*/

	_assert_(start==-1 || start==+1);
	_assert_(this->inputs);

	/*Are we on the the boundary we want to be?*/
	if(start==-1 && !IsOnBase())    return;
	if(start==+1 && !IsOnSurface()) return;

	/*Get original input*/
	ElementInput* input  = this->inputs->GetControlInputData(enum_type,"value");
	if(input->ObjectEnum()!=PentaInputEnum) _error_("not supported yet");
	PentaInput* pentainput = xDynamicCast<PentaInput*>(input);
	ElementInput* input2 = this->inputs->GetControlInputData(enum_type,"savedvalues");
	if(input->ObjectEnum()!=PentaInputEnum) _error_("not supported yet");
	PentaInput* pentainput2= xDynamicCast<PentaInput*>(input2);
	/*FIXME: this should not be necessary*/
	ElementInput* input3 = this->inputs->GetControlInputData(enum_type,"gradient");
	if(input->ObjectEnum()!=PentaInputEnum) _error_("not supported yet");
	PentaInput* pentainput3= xDynamicCast<PentaInput*>(input3);

	int lidlist[NUMVERTICES];
	this->GetVerticesLidList(&lidlist[0]);
	pentainput->Serve(NUMVERTICES,&lidlist[0]);
	pentainput2->Serve(NUMVERTICES,&lidlist[0]);
	pentainput3->Serve(NUMVERTICES,&lidlist[0]);

	if(pentainput->GetInterpolation()==P1Enum){

		/*Extrude values first*/
		IssmDouble extrudedvalues[NUMVERTICES];
		IssmDouble extrudedvalues2[NUMVERTICES];
		IssmDouble extrudedvalues3[NUMVERTICES];

		this->GetInputListOnVertices(&extrudedvalues[0],pentainput,0.);
		this->GetInputListOnVertices(&extrudedvalues2[0],pentainput2,0.);
		this->GetInputListOnVertices(&extrudedvalues3[0],pentainput3,0.);

		if(start==-1){
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i+NUMVERTICES2D]=extrudedvalues[i];
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues2[i+NUMVERTICES2D]=extrudedvalues2[i];
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues3[i+NUMVERTICES2D]=extrudedvalues3[i]/2.; /*FIXME: this is just for NR*/
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues3[i]=extrudedvalues3[i]/2.; /*FIXME: this is just for NR*/
		}
		else{
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i]=extrudedvalues[i+NUMVERTICES2D];
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues2[i]=extrudedvalues2[i+NUMVERTICES2D];
		}

		/*Propagate to other Pentas*/
		Penta* penta=this;
		for(;;){

			if(penta->IsOnSurface() && start==-1){ /*FIXME: this is just for NR*/
				for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues3[i+NUMVERTICES2D]=0.;
			}

			/*Add input of the basal element to penta->inputs*/
			int vertexlids[NUMVERTICES];
			penta->GetVerticesLidList(&vertexlids[0]);
			pentainput->SetInput(P1Enum,NUMVERTICES,&vertexlids[0],&extrudedvalues[0]);
			pentainput2->SetInput(P1Enum,NUMVERTICES,&vertexlids[0],&extrudedvalues2[0]);
			if(start==-1 && !penta->IsOnBase()){
				pentainput3->SetInput(P1Enum,NUMVERTICES,&vertexlids[0],&extrudedvalues3[0]);
			}

			/*Stop if we have reached the surface/base*/
			if(start==-1 && penta->IsOnSurface()) break;
			if(start==+1 && penta->IsOnBase())    break;

			/*get upper/lower Penta*/
			if(start==-1) penta=penta->GetUpperPenta();
			else          penta=penta->GetLowerPenta();
			_assert_(penta->Id()!=this->id);
		}
	}
	else{
		_error_("not implemented yet");
	}
}
/*}}}*/
void       Penta::InputExtrude(int enum_type,int start){/*{{{*/

	_assert_(start==-1 || start==+1);
	_assert_(this->inputs);

	/*Are we on the the boundary we want to be?*/
	if(start==-1 && !IsOnBase())    return;
	if(start==+1 && !IsOnSurface()) return;

	/*Get original input*/
	Input* input = this->GetInput(enum_type);
	if(input->ObjectEnum()!=PentaInputEnum) _error_("not supported yet");
	PentaInput* pentainput = xDynamicCast<PentaInput*>(input);

	if(pentainput->GetInterpolation()==P1Enum || pentainput->GetInterpolation()==P1DGEnum){
		/*Extrude values first*/
		IssmDouble extrudedvalues[NUMVERTICES];

		Element::GetInputListOnVertices(&extrudedvalues[0],enum_type);
		if(start==-1){
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i+NUMVERTICES2D]=extrudedvalues[i];
		}
		else{
			for(int i=0;i<NUMVERTICES2D;i++) extrudedvalues[i]=extrudedvalues[i+NUMVERTICES2D];
		}

		/*Propagate to other Pentas*/
		Penta* penta=this;
		for(;;){

			/*Add input of the basal element to penta->inputs*/
			penta->AddInput(enum_type,&extrudedvalues[0],pentainput->GetInterpolation());

			/*Stop if we have reached the surface/base*/
			if(start==-1 && penta->IsOnSurface()) break;
			if(start==+1 && penta->IsOnBase())    break;

			/*get upper/lower Penta*/
			if(start==-1) penta=penta->GetUpperPenta();
			else          penta=penta->GetLowerPenta();
			_assert_(penta->Id()!=this->id);
		}
	}
	else{
		_error_("interpolation "<<EnumToStringx(pentainput->GetInterpolation())<<" not implemented yet (while trying to extrude "<<EnumToStringx(enum_type)<<")");
	}
}
/*}}}*/
void       Penta::InputUpdateFromIoModel(int index,IoModel* iomodel){ /*{{{*/

	/*Intermediaries*/
	int         i,j;
	int         penta_vertex_ids[NUMVERTICES];
	IssmDouble  nodeinputs[NUMVERTICES];
	IssmDouble  cmmininputs[NUMVERTICES];
	IssmDouble  cmmaxinputs[NUMVERTICES];

	IssmDouble  yts;
	bool    control_analysis;
	char**  controls = NULL;
	int     num_control_type,num_responses;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	if(control_analysis) iomodel->FindConstant(&num_control_type,"md.inversion.num_control_parameters");
	if(control_analysis) iomodel->FindConstant(&num_responses,"md.inversion.num_cost_functions");

	/*Recover vertices ids needed to initialize inputs*/
	_assert_(iomodel->elements);
	for(i=0;i<NUMVERTICES;i++){
		penta_vertex_ids[i]=iomodel->elements[NUMVERTICES*index+i]; //ids for vertices are in the elements array from Matlab
	}
}
/*}}}*/
void       Penta::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(this->element_type);

	/*Fetch dof list and allocate solution vector*/
	GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	this->AddInput(enum_type,values,this->element_type);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
void       Penta::InputUpdateFromSolutionOneDofCollapsed(IssmDouble* solution,int enum_type){/*{{{*/

	const int  numdof   = NUMVERTICES;
	const int  numdof2d = NUMVERTICES2D;

	IssmDouble  values[numdof];
	int*    doflist = NULL;
	Penta  *penta   = NULL;

	/*If not on bed, return*/
	if (!IsOnBase()) return;

	/*Get dof list: */
	GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);

	/*Use the dof list to index into the solution vector and extrude it */
	for(int i=0;i<numdof2d;i++){
		values[i]         =solution[doflist[i]];
		values[i+numdof2d]=values[i];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Start looping over all elements above current element and update all inputs*/
	penta=this;
	for(;;){
		/*Add input to the element: */
		penta->AddInput(enum_type,values,P1Enum);

		/*Stop if we have reached the surface*/
		if (penta->IsOnSurface()) break;

		/* get upper Penta*/
		penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
	}

	/*Free resources:*/
	xDelete<int>(doflist);
}
/*}}}*/
void       Penta::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

	const int   numdof         = NUMVERTICES;
	int        *doflist        = NULL;
	IssmDouble  values[numdof];
	int         lidlist[NUMVERTICES];

	GetVerticesLidList(&lidlist[0]);

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){
		case VertexLIdEnum:
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Lid()];
			}
			/*update input*/
			inputs->SetPentaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			return;

		case VertexPIdEnum:
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Pid()];
			}
			/*update input*/
			inputs->SetPentaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			return;

		case VertexSIdEnum:
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Sid()];
			}
			/*update input*/
			inputs->SetPentaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			return;

		case NodesEnum:
			/*Get dof list: */
			GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

			/*Use the dof list to index into the vector: */
			for(int i=0;i<numdof;i++){
				values[i]=vector[doflist[i]];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
			}
			/*Add input to the element: */
			inputs->SetPentaInput(name,P1Enum,NUMVERTICES,lidlist,values);

			/*Free resources:*/
			xDelete<int>(doflist);
			return;

	  case NodeSIdEnum:
			for(int i=0;i<NUMVERTICES;i++){
				values[i]=vector[nodes[i]->Sid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
			}
			/*Add input to the element: */
			inputs->SetPentaInput(name,P1Enum,NUMVERTICES,lidlist,values);

			/*Free resources:*/
			xDelete<int>(doflist);
			return;

	  default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}
}
/*}}}*/
bool       Penta::IsIcefront(void){/*{{{*/

	bool isicefront;
	int i,nrice;
   IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	isicefront=false;
	if(IsIceInElement()){
		nrice=0;
		for(i=0;i<NUMVERTICES2D;i++)
			if(ls[i]<0.) nrice++;
		if(nrice==1) isicefront= true;
	}
	return isicefront;
}/*}}}*/
bool       Penta::IsNodeOnShelfFromFlags(IssmDouble* flags){/*{{{*/

	int  i;
	bool shelf=false;

	for(i=0;i<NUMVERTICES;i++){
		if (flags[vertices[i]->Pid()]<0.){
			shelf=true;
			break;
		}
	}
	return shelf;
}
/*}}}*/
bool       Penta::IsZeroLevelset(int levelset_enum){/*{{{*/

	bool        iszerols;
	IssmDouble  ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&ls[0],levelset_enum);

	/*If the level set has always same sign, there is no ice front here*/
	iszerols = false;
	if(IsIceInElement()){
		if(ls[0]*ls[1]<0. || ls[0]*ls[2]<0. || (ls[0]*ls[1]*ls[2]==0. && ls[0]*ls[1]+ls[0]*ls[2]+ls[1]*ls[2]<=0.)){
			iszerols = true;
		}
	}
	return iszerols;
}
/*}}}*/
void       Penta::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetTriaJacobianDeterminant(pJdet,xyz_list_base,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantLine(IssmDouble* pJdet,IssmDouble* xyz_list_line,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_line,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list_quad,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetQuadJacobianDeterminant(pJdet,xyz_list_quad,(GaussPenta*)gauss);

}
/*}}}*/
void       Penta::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_top,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetTriaJacobianDeterminant(pJdet,xyz_list_top,(GaussPenta*)gauss);

}
/*}}}*/
IssmDouble Penta::MassFlux(IssmDouble* segment){/*{{{*/

	IssmDouble mass_flux=0;

	if(!IsOnBase()) return mass_flux;

	/*Depth Averaging Vx and Vy*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Spawn Tria element from the base of the Penta: */
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	mass_flux=tria->MassFlux(segment);
	delete tria->material; delete tria;

	/*clean up and return*/
	return mass_flux;
}
/*}}}*/
IssmDouble Penta::MassFlux(IssmDouble x1, IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id){/*{{{*/

	IssmDouble mass_flux=0;

	if(!IsOnBase()) return mass_flux;

	/*Depth Averaging Vx and Vy*/
	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Spawn Tria element from the base of the Penta: */
	Tria* tria=(Tria*)SpawnTria(0,1,2);
	mass_flux=tria->MassFlux(x1,y1,x2,y2,segment_id);
	delete tria->material; delete tria;

	/*clean up and return*/
	return mass_flux;
}
/*}}}*/
IssmDouble Penta::MinEdgeLength(IssmDouble* xyz_list){/*{{{*/
	/*Return the minimum lenght of the nine egdes of the penta*/

	int    i,node0,node1;
	int    edges[9][2]={{0,1},{0,2},{1,2},{3,4},{3,5},{4,5},{0,3},{1,4},{2,5}}; //list of the nine edges
	IssmDouble length;
	IssmDouble minlength=-1;

	for(i=0;i<9;i++){
		/*Find the two nodes for this edge*/
		node0=edges[i][0];
		node1=edges[i][1];

		/*Compute the length of this edge and compare it to the minimal length*/
		length=sqrt(pow(xyz_list[node0*3+0]-xyz_list[node1*3+0],2)+pow(xyz_list[node0*3+1]-xyz_list[node1*3+1],2)+pow(xyz_list[node0*3+2]-xyz_list[node1*3+2],2));
		if(length<minlength || minlength<0) minlength=length;
	}

	return minlength;
}
/*}}}*/
void	      Penta::MovingFrontalVelocity(void){/*{{{*/

	if(!this->IsOnBase()) return;
	int        domaintype, calvinglaw, i;
	IssmDouble v[3],w[3],c[3],m[3],dlsf[3];
	IssmDouble norm_dlsf, norm_calving, calvingrate, meltingrate, groundedice;
	IssmDouble migrationmax, calvinghaf, heaviside, haf_eps;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble movingfrontvx[NUMVERTICES];
	IssmDouble movingfrontvy[NUMVERTICES];
	IssmDouble vel;
	int        dim=2;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* calvingratex_input = NULL;
	Input* calvingratey_input = NULL;
	Input* calvingrate_input  = NULL;
	Input* meltingrate_input  = NULL;

	/*Load levelset function gradients*/
	Input *lsf_slopex_input = this->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
	Input *lsf_slopey_input = this->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);
	Input *vx_input         = this->GetInput(VxAverageEnum);              _assert_(vx_input);
	Input *vy_input         = this->GetInput(VyAverageEnum);              _assert_(vy_input);
	Input *gr_input         = this->GetInput(MaskOceanLevelsetEnum);      _assert_(gr_input);

	/*Get problem dimension and whether there is moving front or not*/
	this->FindParam(&calvinglaw,CalvingLawEnum);
	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingVonmisesEnum:
		case CalvingLevermannEnum:
		case CalvingCalvingMIPEnum:
			calvingratex_input=this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
			calvingratey_input=this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
			meltingrate_input = this->GetInput(CalvingMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		case CalvingMinthicknessEnum:
		case CalvingHabEnum:
		case CalvingCrevasseDepthEnum:
			meltingrate_input = this->GetInput(CalvingMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		case CalvingDev2Enum:
			this->FindParam(&calvinghaf,CalvingHeightAboveFloatationEnum);
			calvingrate_input = this->GetInput(CalvingCalvingrateEnum);     _assert_(calvingrate_input);
			meltingrate_input = this->GetInput(CalvingMeltingrateEnum);     _assert_(meltingrate_input);
			break;
		default:
			_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
	}

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Advection */
		vx_input->GetInputValue(&v[0],&gauss);
		vy_input->GetInputValue(&v[1],&gauss);
		gr_input->GetInputValue(&groundedice,&gauss);
      lsf_slopex_input->GetInputValue(&dlsf[0],&gauss);
      lsf_slopey_input->GetInputValue(&dlsf[1],&gauss);
      norm_dlsf=sqrt(dlsf[0]*dlsf[0] + dlsf[1]*dlsf[1]);

		/*Get calving speed*/
		switch(calvinglaw){
			/*"Contiuous" calving*/
			case DefaultCalvingEnum:
			case CalvingVonmisesEnum:
			case CalvingLevermannEnum:
			case CalvingCalvingMIPEnum:
				calvingratex_input->GetInputValue(&c[0],&gauss);
				calvingratey_input->GetInputValue(&c[1],&gauss);
				meltingrate_input->GetInputValue(&meltingrate,&gauss);
				if(groundedice<0) meltingrate = 0.;
				m[0]=meltingrate*dlsf[0]/norm_dlsf;
				m[1]=meltingrate*dlsf[1]/norm_dlsf;
				break;

			/*"Discrete" calving*/
			case CalvingMinthicknessEnum:
			case CalvingHabEnum:
			case CalvingCrevasseDepthEnum:
				meltingrate_input->GetInputValue(&meltingrate,&gauss);
				if(groundedice<0) meltingrate = 0.;

				if(norm_dlsf>1.e-10)
				 for(i=0;i<dim;i++){
					 c[i]=0.;
					 m[i]=meltingrate*dlsf[i]/norm_dlsf;
				 }
				else
				 for(i=0;i<dim;i++){
					 c[i]=0.;
					 m[i]=0.;
				 }
				break;

			case CalvingDev2Enum:
				  {
					calvingrate_input->GetInputValue(&calvingrate,&gauss);
					meltingrate_input->GetInputValue(&meltingrate,&gauss);
					gr_input->GetInputValue(&groundedice,&gauss);

					//idea: no retreat on ice above critical calving height "calvinghaf" . Limit using regularized Heaviside function.
					vel=sqrt(v[0]*v[0] + v[1]*v[1]);
					haf_eps=10.;
					if(groundedice-calvinghaf<=-haf_eps){
						// ice floats freely below calvinghaf: calve freely
						// undercutting has no effect:
						meltingrate=0.;
					}
					else if(groundedice-calvinghaf>=haf_eps){
						// ice is well above calvinghaf -> no calving back, i.e. limit calving rate to ice velocity
						calvingrate=min(calvingrate,vel);
						// ice is almost grounded: frontal undercutting has maximum effect (do nothing).
					}
					else{ // ice is close to calvinghaf: smooth transition between limitation and free calving.
						//heaviside: 0 for floating, 1 for grounded
						heaviside=(groundedice-calvinghaf+haf_eps)/(2.*haf_eps) + sin(M_PI*(groundedice-calvinghaf)/haf_eps)/(2.*M_PI);
						calvingrate=heaviside*(min(calvingrate,vel)-calvingrate)+calvingrate;
						meltingrate=heaviside*meltingrate+0.;
					}

					if(norm_dlsf>1.e-10)
					 for(i=0;i<dim;i++){
						 c[i]=calvingrate*dlsf[i]/norm_dlsf;
						 m[i]=meltingrate*dlsf[i]/norm_dlsf;
					 }
					else
					 for(i=0;i<dim;i++){
						 c[i]=0.;
						 m[i]=0.;
					 }
					break;
				  }

			default:
				_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
		}

		for(i=0;i<dim;i++) w[i]=v[i]-c[i]-m[i];

		movingfrontvx[iv] = w[0];
		movingfrontvy[iv] = w[1];		
	}

	/*Add input*/
	this->AddInput(MovingFrontalVxEnum,&movingfrontvx[0],P1DGEnum);
	this->AddInput(MovingFrontalVyEnum,&movingfrontvy[0],P1DGEnum);
	this->InputExtrude(MovingFrontalVxEnum,-1);
	this->InputExtrude(MovingFrontalVyEnum,-1);
}
/*}}}*/
Gauss*     Penta::NewGauss(void){/*{{{*/
	return new GaussPenta();
}
/*}}}*/
Gauss*     Penta::NewGauss(int order){/*{{{*/
	return new GaussPenta(order,order);
}
/*}}}*/
Gauss*     Penta::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/

	IssmDouble  area_coordinates[4][3];

	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,4);

	return new GaussPenta(area_coordinates,order_horiz,order_vert);
}
/*}}}*/
Gauss*     Penta::NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order){/*{{{*/
	return new GaussPenta(point1,fraction1,fraction2,mainlyfloating,order);
}
/*}}}*/
Gauss*     Penta::NewGaussBase(int order){/*{{{*/
	return new GaussPenta(0,1,2,order);
}
/*}}}*/
Gauss*     Penta::NewGaussBase(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz){/*{{{*/

	IssmDouble  area_coordinates[2][3];

	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,2);

	return new GaussPenta(area_coordinates,order_horiz);
}
/*}}}*/
Gauss*     Penta::NewGaussLine(int vertex1,int vertex2,int order){/*{{{*/
	return new GaussPenta(vertex1,vertex2,order);
}
/*}}}*/
Gauss*     Penta::NewGaussTop(int order){/*{{{*/
	return new GaussPenta(3,4,5,order);
}
/*}}}*/
void       Penta::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->element_type);

}
/*}}}*/
void       Penta::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,this->element_type);

}
/*}}}*/
void       Penta::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsMINIDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,P1bubbleEnum);

}
/*}}}*/
void       Penta::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->PressureInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsP1(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,P1Enum);

}
/*}}}*/
void       Penta::NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussPenta*)gauss,P1Enum);

}
/*}}}*/
void       Penta::NodalFunctionsP2(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,P2Enum);

}
/*}}}*/
void       Penta::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Penta::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->GetNodalFunctions(basis,(GaussPenta*)gauss,this->TensorInterpolation());

}
/*}}}*/
int        Penta::NodalValue(IssmDouble* pvalue, int index, int natureofdataenum){/*{{{*/

	int i;
	int found=0;
	IssmDouble value;
	GaussPenta gauss;

	/*First, serarch the input: */
	Input* data=this->GetInput(natureofdataenum);

	/*figure out if we have the vertex id: */
	found=0;
	for(i=0;i<NUMVERTICES;i++){
		if(index==vertices[i]->Id()){
			/*Do we have natureofdataenum in our inputs? :*/
			if(data){
				/*ok, we are good. retrieve value of input at vertex :*/
				gauss.GaussVertex(i);
				data->GetInputValue(&value,&gauss);
				found=1;
				break;
			}
		}
	}

	if(found)*pvalue=value;
	return found;
}
/*}}}*/
void       Penta::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/
	TriangleFacetNormal(bed_normal, xyz_list);
	/*Bed normal is opposite to surface normal*/
	for (int i = 0; i < 3; ++i) bed_normal[i] *= (-1);
}
/*}}}*/
void       Penta::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/

	TriangleFacetNormal(normal, xyz_list);
}
/*}}}*/
void       Penta::NormalSectionBase(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/
	LineSectionNormal(normal, xyz_list);
}
/*}}}*/
void       Penta::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/
	TriangleFacetNormal(top_normal, xyz_list);
}
/*}}}*/
int        Penta::NumberofNodesPressure(void){/*{{{*/
	return PentaRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int        Penta::NumberofNodesVelocity(void){/*{{{*/
	return PentaRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
int        Penta::ObjectEnum(void){/*{{{*/

	return PentaEnum;

}
/*}}}*/
void       Penta::PotentialUngrounding(Vector<IssmDouble>* potential_ungrounding){/*{{{*/

	IssmDouble  h[NUMVERTICES],r[NUMVERTICES],gl[NUMVERTICES];
	IssmDouble  bed_hydro;
	IssmDouble  rho_water,rho_ice,density;

	/*material parameters: */
	rho_water=FindParam(MaterialsRhoSeawaterEnum);
	rho_ice=FindParam(MaterialsRhoIceEnum);
	density=rho_ice/rho_water;
	Element::GetInputListOnVertices(&h[0],ThicknessEnum);
	Element::GetInputListOnVertices(&r[0],BedEnum);
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*go through vertices, and figure out which ones are on the ice sheet, and want to unground: */
	for(int i=0;i<NUMVERTICES;i++){
		/*Find if grounded vertices want to start floating*/
		if (gl[i]>0.){
			bed_hydro=-density*h[i];
			if(bed_hydro>r[i]){
				/*Vertex that could potentially unground, flag it*/
				potential_ungrounding->SetValue(vertices[i]->Pid(),1,INS_VAL);
			}
		}
	}
}
/*}}}*/
int        Penta::PressureInterpolation(void){/*{{{*/
	return PentaRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void       Penta::Recover3DMOLHOInput(int targetVel_enum, int numnodes, IssmDouble* vb,  IssmDouble* vsh, IssmDouble* n, IssmDouble* H, IssmDouble* s){/*{{{*/
   /* Recover the velocity acording to v=vb+(1-\zeta^{n+1})vsh, where \zeta=(s-z)/H
    * The variables vb, vsh, n, H and s are all from the 2D horizontal mesh(Tria), with "numnodes" DOFs
    * To project to penta the DOFs are doubled in size
    *
    */
   _assert_(this->inputs);
   if(!IsOnBase()) return;
   else{
      if(targetVel_enum==VxEnum || targetVel_enum==VyEnum){
         IssmDouble vel[NUMVERTICES2D*5];
         IssmDouble* xyz_list = NULL;
			IssmDouble zi;
         Penta* penta = this;
			numnodes = penta->NumberofNodes(P1xP4Enum);
         _assert_(NUMVERTICES2D*5-numnodes==0);
			GaussPenta gauss;

         for(;;){
            penta->GetVerticesCoordinates(&xyz_list);
            for(int i=0;i<NUMVERTICES2D;i++) {
					for (int j=0;j<5;j++){
						/* Get z-coordinate of the current node */
						gauss.GaussNode(P1xP4Enum, i+j*NUMVERTICES2D);
						zi = this->GetZcoord(xyz_list, &gauss);
	               vel[i+j*NUMVERTICES2D] = vb[i] + vsh[i]*(1.0-pow((s[i]-zi)/H[i], (n[i]+1)));
					}
            }
				xDelete<IssmDouble>(xyz_list);

            /*Add to the bottom side of the element*/
            penta->AddInput(targetVel_enum,&vel[0],P1xP4Enum);
            if (penta->IsOnSurface()) break;
            penta=penta->GetUpperPenta(); _assert_(penta->Id()!=this->id);
         }
      }
      else _error_("not implemented yet");
   }

}
/*}}}*/
void       Penta::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	int analysis_type;
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	if(pe){
		if(analysis_type==StressbalanceAnalysisEnum){
			if(this->element_type==MINIcondensedEnum){
				int approximation;
				this->Element::GetInputValue(&approximation,ApproximationEnum);
				if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
					//Do nothing, condensation already done in PVectorCoupling
				}
				else{
					int indices[3]={18,19,20};
					pe->StaticCondensation(Ke,3,&indices[0]);
				}
			}
			else if(this->element_type==P1bubblecondensedEnum){
				int size   = nodes[6]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
				int offset = 0;
				for(int i=0;i<6;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
				int* indices=xNew<int>(size);
				for(int i=0;i<size;i++) indices[i] = offset+i;
				pe->StaticCondensation(Ke,size,indices);
				xDelete<int>(indices);
			}
		}
	}

	if(Ke){
		if(analysis_type==StressbalanceAnalysisEnum){
			int approximation;
			this->Element::GetInputValue(&approximation,ApproximationEnum);
			if(approximation==HOFSApproximationEnum || approximation==SSAFSApproximationEnum){
				//Do nothing condensatino already done for Stokes part
			}
			else{
				if(this->element_type==MINIcondensedEnum){
					int indices[3]={18,19,20};
					Ke->StaticCondensation(3,&indices[0]);
				}
				else if(this->element_type==P1bubblecondensedEnum){
					int size   = nodes[6]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
					int offset = 0;
					for(int i=0;i<6;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
					int* indices=xNew<int>(size);
					for(int i=0;i<size;i++) indices[i] = offset+i;
					Ke->StaticCondensation(size,indices);
					xDelete<int>(indices);
				}
			}
		}
	}
}
/*}}}*/
void       Penta::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int          approximation;
	int          numindices;
	int         *indices = NULL;
	bool		isNitsche;
	IssmDouble   slopex,slopey,groundedice;
	IssmDouble   xz_plane[6];
	IssmDouble*  vertexapproximation= NULL;

	/*For FS only: we want the CS to be tangential to the bedrock*/
	this->Element::GetInputValue(&approximation,ApproximationEnum);
	this->FindParam(&isNitsche,FlowequationIsNitscheEnum);

	if(!IsOnBase() || (approximation!=FSApproximationEnum && approximation!=SSAFSApproximationEnum &&  approximation!=HOFSApproximationEnum)) return;

	/*Get number of nodes for velocity only and base*/
	BasalNodeIndices(&numindices,&indices,this->VelocityInterpolation());

	/*Get inputs*/
	Input* slopex_input=this->GetInput(BedSlopeXEnum); _assert_(slopex_input);
	Input* slopey_input=this->GetInput(BedSlopeYEnum); _assert_(slopey_input);
	Input* groundedicelevelset_input=this->GetInput(MaskOceanLevelsetEnum); _assert_(groundedicelevelset_input);

	/*Loop over basal nodes and update their CS*/
	GaussPenta gauss;
	for(int i=0;i<numindices;i++){//FIXME
		gauss.GaussNode(this->VelocityInterpolation(),indices[i]);

		slopex_input->GetInputValue(&slopex,&gauss);
		slopey_input->GetInputValue(&slopey,&gauss);
		groundedicelevelset_input->GetInputValue(&groundedice,&gauss);

		/*New X axis          New Z axis*/
		xz_plane[0]=1.;       xz_plane[3]=-slopex;
		xz_plane[1]=0.;       xz_plane[4]=-slopey;
		xz_plane[2]=slopex;   xz_plane[5]=1.;

		if (!isNitsche) {
			if(groundedice>=0){
				if(this->nodes[indices[i]]->GetApproximation()==FSvelocityEnum){
					this->nodes[indices[i]]->DofInSSet(2); //vz
				}
				else if(this->nodes[indices[i]]->GetApproximation()==SSAFSApproximationEnum || this->nodes[indices[i]]->GetApproximation()==HOFSApproximationEnum){
					this->nodes[indices[i]]->DofInSSet(4); //vz
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[indices[i]]->GetApproximation())<<" not supported yet");
			}
			else{
				if(this->nodes[indices[i]]->GetApproximation()==FSvelocityEnum){
					this->nodes[indices[i]]->DofInFSet(2); //vz
				}
				else if(this->nodes[indices[i]]->GetApproximation()==SSAFSApproximationEnum || this->nodes[indices[i]]->GetApproximation()==HOFSApproximationEnum){
					this->nodes[indices[i]]->DofInFSet(4); //vz
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[indices[i]]->GetApproximation())<<" not supported yet");
			}
		}
		XZvectorsToCoordinateSystem(&this->nodes[indices[i]]->coord_system[0][0],&xz_plane[0]);
		this->nodes[indices[i]]->isrotated = true;
	}

	/*cleanup*/
	xDelete<int>(indices);
}
/*}}}*/
void       Penta::ResetHooks(){/*{{{*/

	if(this->nodes) xDelete<Node*>(this->nodes);
	this->nodes=NULL;
	this->vertices=NULL;
	this->material=NULL;
	this->verticalneighbors=NULL;
	this->parameters=NULL;

	//deal with ElementHook mother class
	for(int i=0;i<this->numanalyses;i++) if(this->hnodes[i]) this->hnodes[i]->reset();
	this->hvertices->reset();
	this->hmaterial->reset();
	if(this->hneighbors) this->hneighbors->reset();

}
/*}}}*/
void       Penta::SetElementInput(int enum_in,IssmDouble value){/*{{{*/

	this->SetElementInput(this->inputs,enum_in,value);

}
/*}}}*/
void       Penta::SetElementInput(Inputs* inputs,int enum_in,IssmDouble value){/*{{{*/

	_assert_(inputs);
	inputs->SetPentaInput(enum_in,P0Enum,this->lid,value);

}
/*}}}*/
void       Penta::SetElementInput(Inputs* inputs,int numindices,int* indices,IssmDouble* values,int enum_in){/*{{{*/

	_assert_(inputs);
	inputs->SetPentaInput(enum_in,P1Enum,numindices,indices,values);

}
/*}}}*/
void       Penta::SetElementInput(int enum_in,IssmDouble value,int type){/*{{{*/

	if(type==P0Enum){
		this->inputs->SetPentaInput(enum_in,P0Enum,this->lid,value);
	}
	else if(type==P1Enum){
		IssmDouble values[6]; 
		for(int i=0;i<6;i++)values[i]=value;
		int lidlist[6];
		this->GetVerticesLidList(&lidlist[0]);
		this->inputs->SetPentaInput(enum_in,P1Enum,6,&lidlist[0],&values[0]);
	}
	else _error_("interpolation type not supported yet");
}
/*}}}*/
void       Penta::SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index,int offset, int M, int N){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         lidlist[NUMVERTICES];
	int         idlist[NUMVERTICES],control_init;

	/*Specific case for depth averaged quantities*/
	control_init=control_enum;
	if(control_enum==MaterialsRheologyBbarEnum){
		control_enum=MaterialsRheologyBEnum;
		if(!IsOnBase()) return;
	}
	if(control_enum==DamageDbarEnum){
		control_enum=DamageDEnum;
		if(!IsOnBase()) return;
	}

	/*Prepare index list*/
	ElementInput* input=this->inputs->GetControlInputData(control_enum,"value");   _assert_(input);
	this->GetVerticesLidList(&lidlist[0]);
	GradientIndexing(&idlist[0],control_index);

	/*Get values on vertices*/
	if(input->ObjectEnum()==PentaInputEnum && input->GetInputInterpolationType()==P1Enum){
		_assert_(N==1);
		for(int i=0;i<NUMVERTICES;i++){
			values[i] = vector[idlist[i]];
		}
		input->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
	}
	else if(input->ObjectEnum()==PentaInputEnum && input->GetInputInterpolationType()==P0Enum){
		_assert_(N==1);
		input->SetInput(P0Enum,this->lid,vector[idlist[0]]);
	}
	else if(input->ObjectEnum()==TransientInputEnum){
		for(int n=0;n<N;n++){
			_error_("not implemented");
			//Input* new_input = new PentaInput(control_enum,values,P1Enum);
			//controlinput->SetInput(new_input,n);
			//controlinput->Configure(parameters);
		}
	}
	else _error_("Type not supported");

	/*Extrude depending on the control*/
	if(control_init==MaterialsRheologyBbarEnum){
		this->ControlInputExtrude(control_enum,-1);
	}
	if(control_init==DamageDbarEnum){
		this->ControlInputExtrude(control_enum,-1);
	}
}
/*}}}*/
void       Penta::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	this->element_type=this->element_type_list[analysis_counter];

	/*Pick up nodes*/
	if(this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;

}
/*}}}*/
void       Penta::SetTemporaryElementType(int element_type_in){/*{{{*/
	this->element_type=element_type_in;
}
/*}}}*/
Element*   Penta::SpawnBasalElement(bool depthaverage_materials){/*{{{*/

	_assert_(this->IsOnBase());

	if(depthaverage_materials){
		switch(this->material->ObjectEnum()){
			case MaticeEnum:
				this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
				if(this->material->IsDamage())this->InputDepthAverageAtBase(DamageDEnum,DamageDbarEnum);
				if(this->material->IsEnhanced())this->InputDepthAverageAtBase(MaterialsRheologyEEnum,MaterialsRheologyEbarEnum);
				break;
			case MatestarEnum:
				this->InputDepthAverageAtBase(MaterialsRheologyBEnum,MaterialsRheologyBbarEnum);
				this->InputDepthAverageAtBase(MaterialsRheologyEcEnum,MaterialsRheologyEcbarEnum);
				this->InputDepthAverageAtBase(MaterialsRheologyEsEnum,MaterialsRheologyEsbarEnum);
				break;
			default:
				_error_("not supported yet");
		}
	}

	return SpawnTria(0,1,2);
}
/*}}}*/
Element*   Penta::SpawnTopElement(void){/*{{{*/

	_assert_(this->IsOnSurface());
	return SpawnTria(3,4,5);
}
/*}}}*/
Tria*      Penta::SpawnTria(int index1,int index2,int index3){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	this->parameters->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Create Tria*/
	Tria* tria=new Tria();
	tria->id=this->id;
	tria->sid=this->sid;
	tria->lid=this->lid;
	tria->parameters=this->parameters;
	tria->inputs=this->inputs;
	tria->element_type=P1Enum; //Only P1 CG for now (TO BE CHANGED)

	if(index1==0 && index2==1 && index3==2){
		tria->iscollapsed = 1;
		if(this->nodes)tria->nodes=this->nodes;
		if(this->vertices)tria->vertices=this->vertices;
	}
	else if(index1==3 && index2==4 && index3==5){
		tria->iscollapsed = 2;
		if(this->nodes)tria->nodes=&this->nodes[3];
		if(this->vertices)tria->vertices=&this->vertices[3];
	}
	else _error_("Spawn supported only for surface and basal elements");
	/*Spawn material*/
	tria->material=(Material*)this->material->copy2(tria);

	/*Return new Tria*/
	return tria;
}
/*}}}*/
bool       Penta::IsSpawnedElement(void){/*{{{*/

	/*Penta cannot be collapsed elements*/
	return false;

}/*}}}*/
IssmDouble Penta::StabilizationParameter(IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble diameter, IssmDouble kappa){/*{{{*/
	/*Compute stabilization parameter*/
	/*kappa=thermalconductivity/(rho_ice*heatcapacity) for thermal model*/
	/*kappa=enthalpydiffusionparameter for enthalpy model*/

	IssmDouble normu;
	IssmDouble tau_parameter;

	normu=pow(pow(u,2)+pow(v,2)+pow(w,2),0.5);
	if(normu*diameter/(3*2*kappa)<1){
		tau_parameter=pow(diameter,2)/(3*2*2*kappa);
	}
	else tau_parameter=diameter/(2*normu);

	return tau_parameter;
}/*}}}*/
void Penta::StabilizationParameterAnisotropic(IssmDouble* tau_parameter_anisotropic, IssmDouble u, IssmDouble v, IssmDouble w, IssmDouble hx, IssmDouble hy, IssmDouble hz, IssmDouble kappa){/*{{{*/
	/*Compute stabilization parameter*/
	/*kappa=thermalconductivity/(rho_ice*heatcapacity) for thermal model*/
	/*kappa=enthalpydiffusionparameter for enthalpy model*/

	IssmDouble normu,hk,C,area,p;

	/* compute tau for the horizontal direction */
	p=2.; C=3.;
	normu=pow(pow(u,p)+pow(v,p)+pow(w,p),1./p);
	hk=sqrt(pow(hx,2)+pow(hy,2));

	if(normu*hk/(C*2*kappa)<1){
		tau_parameter_anisotropic[0]=pow(hk,2)/(C*2*2*kappa);
	}
	else tau_parameter_anisotropic[0]=hk/(2*normu);

	/* compute tau for the vertical direction */
	hk=hz;
	if(normu*hk/(C*2*kappa)<1){
		tau_parameter_anisotropic[1]=pow(hk,2)/(C*2*2*kappa);
	}
	else{
		tau_parameter_anisotropic[1]=hk/(2*normu);
	}
}
/*}}}*/
void       Penta::StrainRateparallel(){/*{{{*/

	IssmDouble  epsilon[6];
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainparallel[NUMVERTICES];

	/* Get node coordinates and dof list: */
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will need*/
	Input* vx_input=this->GetInput(VxEnum);                                  _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);                                  _assert_(vy_input);
	Input* vz_input=this->GetInput(VzEnum);												_assert_(vz_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate and viscosity: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[3];

		/*strainparallel= Strain rate along the ice flow direction */
		strainparallel[iv]=(vx*vx*(strainxx)+vy*vy*(strainyy)+2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->AddInput(StrainRateparallelEnum,&strainparallel[0],P1DGEnum);
}
/*}}}*/
void       Penta::StrainRateperpendicular(){/*{{{*/

	IssmDouble  epsilon[6];
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainperpendicular[NUMVERTICES];

	/* Get node coordinates and dof list: */
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will need*/
	Input* vx_input=this->GetInput(VxEnum);                                  _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);                                  _assert_(vy_input);
	Input* vz_input=this->GetInput(VzEnum);												_assert_(vz_input);

	/* Start looping on the number of vertices: */
	GaussPenta gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate and viscosity: */
		this->StrainRateFS(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input,vz_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[3];

		/*strainperpendicular= Strain rate perpendicular to the ice flow direction */
		strainperpendicular[iv]=(vx*vx*(strainyy)+vy*vy*(strainxx)-2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->AddInput(StrainRateperpendicularEnum,&strainperpendicular[0],P1DGEnum);
}
/*}}}*/
void       Penta::StressIntensityFactor(){/*{{{*/

	/* Check if we are on the base */
	if(!IsOnBase()) return;

	IssmDouble  ki[6]={0.};
	IssmDouble  const_grav=9.81;
	IssmDouble  rho_ice=900;
	IssmDouble  rho_water=1000;
	IssmDouble  Jdet[3];
	IssmDouble  pressure,vx,vy,vel,deviaxx,deviaxy,deviayy,water_depth,prof,stress_xx,thickness;

	Penta* penta=this;
	for(;;){

		IssmDouble  xyz_list[NUMVERTICES][3];
		/* Get node coordinates and dof list: */
		::GetVerticesCoordinates(&xyz_list[0][0],penta->vertices,NUMVERTICES);

		///*Compute the Jacobian for the vertical integration*/
		Jdet[0]=(xyz_list[3][2]-xyz_list[0][2])*0.5;
		Jdet[1]=(xyz_list[4][2]-xyz_list[1][2])*0.5;
		Jdet[2]=(xyz_list[5][2]-xyz_list[2][2])*0.5;

		/*Retrieve all inputs we will need*/
		Input* vx_input=this->GetInput(VxEnum);                                  _assert_(vx_input);
		Input* vy_input=this->GetInput(VyEnum);                                  _assert_(vy_input);
		Input* vel_input=this->GetInput(VelEnum);                                _assert_(vel_input);
		Input* pressure_input=this->GetInput(PressureEnum);                      _assert_(pressure_input);
		Input* deviaxx_input=this->GetInput(DeviatoricStressxxEnum);             _assert_(deviaxx_input);
		Input* deviaxy_input=this->GetInput(DeviatoricStressxyEnum);             _assert_(deviaxy_input);
		Input* deviayy_input=this->GetInput(DeviatoricStressyyEnum);             _assert_(deviayy_input);
		Input* surface_input=this->GetInput(SurfaceEnum);								_assert_(surface_input);
		Input* thickness_input=this->GetInput(ThicknessEnum);							_assert_(thickness_input);

		/* Start looping on the number of 2D vertices: */
		for(int ig=0;ig<3;ig++){
			GaussPenta* gauss=new GaussPenta(ig,3+ig,11);
			while(gauss->next()){

				/* Get the value we need*/
				pressure_input->GetInputValue(&pressure,gauss);
				vx_input->GetInputValue(&vx,gauss);
				vy_input->GetInputValue(&vy,gauss);
				vel_input->GetInputValue(&vel,gauss);
				deviaxx_input->GetInputValue(&deviaxx,gauss);
				deviaxy_input->GetInputValue(&deviaxy,gauss);
				deviayy_input->GetInputValue(&deviayy,gauss);
				surface_input->GetInputValue(&water_depth,gauss);
				thickness_input->GetInputValue(&thickness,gauss);
				prof=water_depth-penta->GetZcoord(&xyz_list[0][0],gauss);

				/*stress_xx= Deviatoric stress along the ice flow direction plus cryostatic pressure */
				stress_xx=(vx*vx*(deviaxx)+vy*vy*(deviayy)+2*vy*vx*deviaxy)/(vel*vel+1.e-6);

				if(prof<water_depth&prof<thickness){
					/* Compute the local stress intensity factor*/
					ki[ig]+=Jdet[ig]*gauss->weight*stress_xx*StressIntensityIntegralWeight(prof,min(water_depth,thickness),thickness);
				}
			}
			delete gauss;
		}

		/*Stop if we have reached the surface/base*/
		if(penta->IsOnSurface()) break;

		/*get upper Penta*/
		penta=penta->GetUpperPenta();
		_assert_(penta->Id()!=this->id);
	}

	/*Add input*/
	this->AddInput(StressIntensityFactorEnum,&ki[0],P1Enum);
	this->InputExtrude(StressIntensityFactorEnum,-1);
}
/*}}}*/
IssmDouble Penta::SurfaceArea(void){/*{{{*/

	int    approximation;
	IssmDouble S;
	Tria*  tria=NULL;

	/*retrieve inputs :*/
	this->Element::GetInputValue(&approximation,ApproximationEnum);

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Bail out if this element if:
	 * -> Non SSA not on the surface
	 * -> SSA (2d model) and not on bed) */
	if ((approximation!=SSAApproximationEnum && !IsOnSurface()) || (approximation==SSAApproximationEnum && !IsOnBase())){
		return 0;
	}
	else if (approximation==SSAApproximationEnum){

		/*This element should be collapsed into a tria element at its base. Create this tria element,
		 * and compute SurfaceArea*/
		tria=(Tria*)SpawnTria(0,1,2);
		S=tria->SurfaceArea();
		delete tria->material; delete tria;
		return S;
	}
	else{

		tria=(Tria*)SpawnTria(3,4,5);
		S=tria->SurfaceArea();
		delete tria->material; delete tria;
		return S;
	}
}
/*}}}*/
void       Penta::TangentBase(IssmDouble* bed_tangent,IssmDouble* bed_normal){/*{{{*/
	/*
	 To compute the two tangent bed_tangent[0:2] and bed_tangent[3:5] from the given normal vecotr.
	*/
	IssmDouble n1, n2, n3;
	IssmDouble tangent_norm;

	n1 = fabs(bed_normal[0]);
	n2 = fabs(bed_normal[1]);
	n3 = fabs(bed_normal[2]);

	/* For the first tangent, on x-y plane in most cases*/
	if ((n1<=n3) && (n2<=n3)) {
		bed_tangent[0] = 0.0;
		bed_tangent[1] = -bed_normal[2];
		bed_tangent[2] = bed_normal[1];
	}
	else {
		bed_tangent[0] = -bed_normal[1];
		bed_tangent[1] = bed_normal[0];
		bed_tangent[2] = 0.0;
	}
	tangent_norm = sqrt(bed_tangent[0]*bed_tangent[0]+bed_tangent[1]*bed_tangent[1]+bed_tangent[2]*bed_tangent[2]);
	for(int i=0;i<3;i++) bed_tangent[i] = bed_tangent[i]/tangent_norm;
	/* The second tangent*/
	bed_tangent[3] = bed_normal[1]*bed_tangent[2] - bed_normal[2]*bed_tangent[1];
	bed_tangent[4] = bed_normal[2]*bed_tangent[0] - bed_normal[0]*bed_tangent[2];
	bed_tangent[5] = bed_normal[0]*bed_tangent[1] - bed_normal[1]*bed_tangent[0];
	tangent_norm = sqrt(bed_tangent[3]*bed_tangent[3]+bed_tangent[4]*bed_tangent[4]+bed_tangent[5]*bed_tangent[5]);
	for(int i=3;i<6;i++) bed_tangent[i] = bed_tangent[i]/tangent_norm;

	IssmDouble checksum = 0.0;
	for (int i=0;i<3;i++) {
		checksum += bed_normal[i]*bed_tangent[i];
		checksum += bed_normal[i]*bed_tangent[i+3];
		checksum += bed_tangent[i]*bed_tangent[i+3];
	}
	_assert_(fabs(checksum)<1e-10);
}
/*}}}*/
IssmDouble Penta::TimeAdapt(void){/*{{{*/

	/*intermediary: */
	bool       ishydro;
	IssmDouble C;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*get CFL coefficient:*/
	this->parameters->FindParam(&C,TimesteppingCflCoefficientEnum);
	this->parameters->FindParam(&ishydro,TransientIshydrologyEnum);

	/*Get for Vx and Vy, the max of abs value: */
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* vz_input = this->GetInput(VzEnum); _assert_(vz_input);
	IssmDouble maxabsvx = vx_input->GetInputMaxAbs();
	IssmDouble maxabsvy = vy_input->GetInputMaxAbs();
	IssmDouble maxabsvz = vz_input->GetInputMaxAbs();

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	IssmDouble minx=xyz_list[0][0];
	IssmDouble maxx=xyz_list[0][0];
	IssmDouble miny=xyz_list[0][1];
	IssmDouble maxy=xyz_list[0][1];
	IssmDouble minz=xyz_list[0][2];
	IssmDouble maxz=xyz_list[0][2];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<minx) minx=xyz_list[i][0];
		if(xyz_list[i][0]>maxx) maxx=xyz_list[i][0];
		if(xyz_list[i][1]<miny) miny=xyz_list[i][1];
		if(xyz_list[i][1]>maxy) maxy=xyz_list[i][1];
		if(xyz_list[i][2]<minz) minz=xyz_list[i][2];
		if(xyz_list[i][2]>maxz) maxz=xyz_list[i][2];
	}
	IssmDouble dx=maxx-minx;
	IssmDouble dy=maxy-miny;
	IssmDouble dz=maxz-minz;

	/*CFL criterion: */
	IssmDouble dt = C/(maxabsvx/dx+maxabsvy/dy+maxabsvz/dz + 1.e-18);

	/*Check hydro timestep also and take the minimum*/
	if(ishydro){
		if(IsOnBase()){
			Input* vx_input = this->GetInput(HydrologyWaterVxEnum); _assert_(vx_input);
			Input* vy_input = this->GetInput(HydrologyWaterVyEnum); _assert_(vy_input);
			IssmDouble maxabsvx = vx_input->GetInputMaxAbs();
			IssmDouble maxabsvy = vy_input->GetInputMaxAbs();

			/*CFL criterion: */
			IssmDouble dt2 = C/(maxabsvx/dx+maxabsvy/dy + 1.e-18);
			if(dt2<dt) dt = dt2;
		}
	}

	return dt;
}/*}}}*/
IssmDouble Penta::TotalCalvingFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum) || !IsOnBase()) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	int pt1 = 0;
	int pt2 = 1;
	if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

		/*Portion of the segments*/
		s1=gl[2]/(gl[2]-gl[1]);
		s2=gl[2]/(gl[2]-gl[0]);
		if(gl[2]<0.){
			pt1 = 1; pt2 = 0;
		}
		xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
		xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
		xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
		xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
		xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
		xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
	}
	else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

		/*Portion of the segments*/
		s1=gl[0]/(gl[0]-gl[1]);
		s2=gl[0]/(gl[0]-gl[2]);
		if(gl[0]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
		xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
		xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
		xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
		xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
		xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
	}
	else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

		/*Portion of the segments*/
		s1=gl[1]/(gl[1]-gl[0]);
		s2=gl[1]/(gl[1]-gl[2]);
		if(gl[1]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
		xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
		xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
		xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
		xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
		xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
	}
	else{
		_error_("case not possible");
	}

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble calvingratex,calvingratey,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* calvingratex_input=NULL;
	Input* calvingratey_input=NULL;
	calvingratex_input=this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
	calvingratey_input=this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		calvingratex_input->GetInputValue(&calvingratex,gauss);
		calvingratey_input->GetInputValue(&calvingratey,gauss);
		this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(calvingratex*normal[0] + calvingratey*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Penta::TotalCalvingMeltingFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum) || !IsOnBase()) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	int pt1 = 0;
	int pt2 = 1;
	if(gl[0]*gl[1]>0){ //Nodes 0 and 1 are similar, so points must be found on segment 0-2 and 1-2

		/*Portion of the segments*/
		s1=gl[2]/(gl[2]-gl[1]);
		s2=gl[2]/(gl[2]-gl[0]);
		if(gl[2]<0.){
			pt1 = 1; pt2 = 0;
		}
		xyz_front[pt2][0]=xyz_list[2][0]+s1*(xyz_list[1][0]-xyz_list[2][0]);
		xyz_front[pt2][1]=xyz_list[2][1]+s1*(xyz_list[1][1]-xyz_list[2][1]);
		xyz_front[pt2][2]=xyz_list[2][2]+s1*(xyz_list[1][2]-xyz_list[2][2]);
		xyz_front[pt1][0]=xyz_list[2][0]+s2*(xyz_list[0][0]-xyz_list[2][0]);
		xyz_front[pt1][1]=xyz_list[2][1]+s2*(xyz_list[0][1]-xyz_list[2][1]);
		xyz_front[pt1][2]=xyz_list[2][2]+s2*(xyz_list[0][2]-xyz_list[2][2]);
	}
	else if(gl[1]*gl[2]>0){ //Nodes 1 and 2 are similar, so points must be found on segment 0-1 and 0-2

		/*Portion of the segments*/
		s1=gl[0]/(gl[0]-gl[1]);
		s2=gl[0]/(gl[0]-gl[2]);
		if(gl[0]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt1][0]=xyz_list[0][0]+s1*(xyz_list[1][0]-xyz_list[0][0]);
		xyz_front[pt1][1]=xyz_list[0][1]+s1*(xyz_list[1][1]-xyz_list[0][1]);
		xyz_front[pt1][2]=xyz_list[0][2]+s1*(xyz_list[1][2]-xyz_list[0][2]);
		xyz_front[pt2][0]=xyz_list[0][0]+s2*(xyz_list[2][0]-xyz_list[0][0]);
		xyz_front[pt2][1]=xyz_list[0][1]+s2*(xyz_list[2][1]-xyz_list[0][1]);
		xyz_front[pt2][2]=xyz_list[0][2]+s2*(xyz_list[2][2]-xyz_list[0][2]);
	}
	else if(gl[0]*gl[2]>0){ //Nodes 0 and 2 are similar, so points must be found on segment 1-0 and 1-2

		/*Portion of the segments*/
		s1=gl[1]/(gl[1]-gl[0]);
		s2=gl[1]/(gl[1]-gl[2]);
		if(gl[1]<0.){
			pt1 = 1; pt2 = 0;
		}

		xyz_front[pt2][0]=xyz_list[1][0]+s1*(xyz_list[0][0]-xyz_list[1][0]);
		xyz_front[pt2][1]=xyz_list[1][1]+s1*(xyz_list[0][1]-xyz_list[1][1]);
		xyz_front[pt2][2]=xyz_list[1][2]+s1*(xyz_list[0][2]-xyz_list[1][2]);
		xyz_front[pt1][0]=xyz_list[1][0]+s2*(xyz_list[2][0]-xyz_list[1][0]);
		xyz_front[pt1][1]=xyz_list[1][1]+s2*(xyz_list[2][1]-xyz_list[1][1]);
		xyz_front[pt1][2]=xyz_list[1][2]+s2*(xyz_list[2][2]-xyz_list[1][2]);
	}
	else{
		_error_("case not possible");
	}

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSectionBase(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	this->InputDepthAverageAtBase(VxEnum,VxAverageEnum);
	this->InputDepthAverageAtBase(VyEnum,VyAverageEnum);

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble calvingratex,calvingratey,vx,vy,vel,meltingrate,meltingratex,meltingratey,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* calvingratex_input=NULL;
	Input* calvingratey_input=NULL;
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	Input* meltingrate_input=NULL;
	calvingratex_input=this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
	calvingratey_input=this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
	vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
	vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	meltingrate_input=this->GetInput(CalvingMeltingrateEnum); _assert_(meltingrate_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGaussBase(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		calvingratex_input->GetInputValue(&calvingratex,gauss);
		calvingratey_input->GetInputValue(&calvingratey,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;
		meltingrate_input->GetInputValue(&meltingrate,gauss);
		meltingratex=meltingrate*vx/(sqrt(vel)+1.e-14);
		meltingratey=meltingrate*vy/(sqrt(vel)+1.e-14);
		this->JacobianDeterminantLine(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*((calvingratex+meltingratex)*normal[0] + (calvingratey+meltingratey)*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Penta::TotalFloatingBmb(bool scaled){/*{{{*/

	/*The fbmb[kg yr-1] of one element is area[m2] * melting_rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble fbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,floatingmelt,Jdet,scalefactor;
	IssmDouble Total_Fbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement() || !IsOnBase())return 0;

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* floatingmelt_input = this->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(floatingmelt_input);
	Input* gllevelset_input = this->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
	Input* scalefactor_input = NULL;
	if(scaled==true){
		scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
	}
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,1-mainlyfloating,3);
	while(gauss->next()){
		this->JacobianDeterminantBase(&Jdet,&xyz_list[0][0],gauss);
		floatingmelt_input->GetInputValue(&floatingmelt,gauss);
		if(scaled==true){
			scalefactor_input->GetInputValue(&scalefactor,gauss);
		}
		else scalefactor=1;
		fbmb+=floatingmelt*Jdet*gauss->weight*scalefactor;
	}

   Total_Fbmb=rho_ice*fbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Fbmb;
}
/*}}}*/
IssmDouble Penta::TotalGroundedBmb(bool scaled){/*{{{*/

	/*The gbmb[kg yr-1] of one element is area[m2] * gounded melting rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble gbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,groundedmelt,Jdet,scalefactor;
	IssmDouble Total_Gbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement() || !IsOnBase())return 0;

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* groundedmelt_input = this->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedmelt_input);
	Input* gllevelset_input   = this->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
	Input* scalefactor_input  = NULL;
	if(scaled==true){
		scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
	}
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,mainlyfloating,3);
	while(gauss->next()){
		this->JacobianDeterminantBase(&Jdet,&xyz_list[0][0],gauss);
		groundedmelt_input->GetInputValue(&groundedmelt,gauss);
		if(scaled==true){
			scalefactor_input->GetInputValue(&scalefactor,gauss);
		}
		else scalefactor=1;
		gbmb+=groundedmelt*Jdet*gauss->weight*scalefactor;
	}

   Total_Gbmb=rho_ice*gbmb;	        // from volume to mass

	/*Return: */
	delete gauss;
	return Total_Gbmb;
}
/*}}}*/
IssmDouble Penta::TotalSmb(bool scaled){/*{{{*/

	/*The smb[Gt yr-1] of one element is area[m2] * smb [ m ice yr^-1] * rho_ice [kg m-3] / 1e+10^12 */
	IssmDouble base,smb,rho_ice,scalefactor;
	IssmDouble Total_Smb=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

	if(!IsIceInElement() || !IsOnSurface()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
   if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* smb_vertices = xNew<IssmDouble>(NUMVERTICES);
		IssmDouble  weights[NUMVERTICES2D];
		IssmDouble  lsf2d[NUMVERTICES2D];
      IssmDouble f1,f2,phi;
      Element::GetInputListOnVertices(&smb_vertices[0],SmbMassBalanceEnum);
		for(int i=0;i<NUMVERTICES2D;i++) lsf2d[i] = lsf[i];
		GetFractionGeometry2D(weights,&phi,&point,&f1,&f2,&mainlyice,lsf2d);
		smb = 0.0;
		for(int i=0;i<NUMVERTICES2D;i++) smb += weights[i]*smb_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices   = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         /*Compute loop only over lower vertices: i<NUMVERTICES2D*/
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES2D;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
		}
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(smb_vertices);
	}

	else{
		/*Fully ice-covered element*/
		Input* smb_input = this->GetInput(SmbMassBalanceEnum); _assert_(smb_input);
		smb_input->GetInputAverage(&smb);

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

	Total_Smb=rho_ice*base*smb*scalefactor;// smb on element in kg s-1

	/*Return: */
	return Total_Smb;
}
/*}}}*/
IssmDouble Penta::TotalSmbMelt(bool scaled){/*{{{*/

	/*The smbmelt[Gt yr-1] of one element is area[m2] * smb [ m ice yr^-1] * rho_ice [kg m-3] / 1e+10^12 */
	IssmDouble base,smbmelt,rho_ice,scalefactor;
	IssmDouble Total_SmbMelt=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

	if(!IsIceInElement() || !IsOnSurface()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
   if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* smbmelt_vertices = xNew<IssmDouble>(NUMVERTICES);
		IssmDouble  weights[NUMVERTICES2D];
		IssmDouble  lsf2d[NUMVERTICES2D];
      IssmDouble f1,f2,phi;
      Element::GetInputListOnVertices(&smbmelt_vertices[0],SmbMeltEnum);
		for(int i=0;i<NUMVERTICES2D;i++) lsf2d[i] = lsf[i];
		GetFractionGeometry2D(weights,&phi,&point,&f1,&f2,&mainlyice,lsf2d);
		smbmelt = 0.0;
		for(int i=0;i<NUMVERTICES2D;i++) smbmelt += weights[i]*smbmelt_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices   = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         /*Compute loop only over lower vertices: i<NUMVERTICES2D*/
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES2D;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
		}
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(smbmelt_vertices);
	}

	else{
		/*Fully ice-covered element*/
		Input* smbmelt_input = this->GetInput(SmbMeltEnum); _assert_(smbmelt_input);
		smbmelt_input->GetInputAverage(&smbmelt);

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

	Total_SmbMelt=rho_ice*base*smbmelt*scalefactor;// smbmelt on element in kg s-1

	/*Return: */
	return Total_SmbMelt;
}
/*}}}*/
IssmDouble Penta::TotalSmbRefreeze(bool scaled){/*{{{*/

	/*The smbrefreeze[Gt yr-1] of one element is area[m2] * smb [ m ice yr^-1] * rho_ice [kg m-3] / 1e+10^12 */
	IssmDouble base,smbrefreeze,rho_ice,scalefactor;
	IssmDouble Total_SmbRefreeze=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

	if(!IsIceInElement() || !IsOnSurface()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
   if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* smbrefreeze_vertices = xNew<IssmDouble>(NUMVERTICES);
		IssmDouble  weights[NUMVERTICES2D];
		IssmDouble  lsf2d[NUMVERTICES2D];
      IssmDouble f1,f2,phi;
      Element::GetInputListOnVertices(&smbrefreeze_vertices[0],SmbRefreezeEnum);
		for(int i=0;i<NUMVERTICES2D;i++) lsf2d[i] = lsf[i];
		GetFractionGeometry2D(weights,&phi,&point,&f1,&f2,&mainlyice,lsf2d);
		smbrefreeze = 0.0;
		for(int i=0;i<NUMVERTICES2D;i++) smbrefreeze += weights[i]*smbrefreeze_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices   = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         /*Compute loop only over lower vertices: i<NUMVERTICES2D*/
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES2D;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
		}
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(smbrefreeze_vertices);
	}

	else{
		/*Fully ice-covered element*/
		Input* smbrefreeze_input = this->GetInput(SmbRefreezeEnum); _assert_(smbrefreeze_input);
		smbrefreeze_input->GetInputAverage(&smbrefreeze);

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

	Total_SmbRefreeze=rho_ice*base*smbrefreeze*scalefactor;// smbrefreeze on element in kg s-1

	/*Return: */
	return Total_SmbRefreeze;
}
/*}}}*/
void       Penta::Update(Inputs* inputs,int index,IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){ /*{{{*/

	/*Intermediaries*/
	int        i;
	int        penta_vertex_ids[6];
	IssmDouble nodeinputs[6];
	IssmDouble yts;
	bool       dakota_analysis;
	int        numnodes;
	int*       penta_node_ids = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");

	/*Checks if debuging*/
	_assert_(iomodel->elements);
	_assert_(index==this->sid);

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<6;i++) penta_vertex_ids[i]=iomodel->elements[6*index+i]; //ids for vertices are in the elements array from Matlab

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P1Enum:
			numnodes         = 6;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[0]=iomodel->elements[6*index+0];
			penta_node_ids[1]=iomodel->elements[6*index+1];
			penta_node_ids[2]=iomodel->elements[6*index+2];
			penta_node_ids[3]=iomodel->elements[6*index+3];
			penta_node_ids[4]=iomodel->elements[6*index+4];
			penta_node_ids[5]=iomodel->elements[6*index+5];
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes         = 7;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[0]=iomodel->elements[6*index+0];
			penta_node_ids[1]=iomodel->elements[6*index+1];
			penta_node_ids[2]=iomodel->elements[6*index+2];
			penta_node_ids[3]=iomodel->elements[6*index+3];
			penta_node_ids[4]=iomodel->elements[6*index+4];
			penta_node_ids[5]=iomodel->elements[6*index+5];
			penta_node_ids[6]=iomodel->numberofvertices+index+1;
			break;
		case P1xP2Enum:
			numnodes         = 9;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoverticaledgeconnectivity[3*index+2]+1;
			break;
		case P1xP3Enum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+2;
			penta_node_ids[10]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+2;
			penta_node_ids[11]=iomodel->numberofvertices+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+2;
			break;
		case P2xP1Enum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtohorizontaledgeconnectivity[6*index+5]+1;
			break;
		case P1xP4Enum:
			numnodes         = 15;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0]; /*Vertex 1*/
			penta_node_ids[ 1]=iomodel->elements[6*index+1]; /*Vertex 2*/
			penta_node_ids[ 2]=iomodel->elements[6*index+2]; /*Vertex 3*/
			penta_node_ids[ 3]=iomodel->elements[6*index+3]; /*Vertex 4*/
			penta_node_ids[ 4]=iomodel->elements[6*index+4]; /*Vertex 5*/
			penta_node_ids[ 5]=iomodel->elements[6*index+5]; /*Vertex 6*/
			penta_node_ids[ 6]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+0]+1; /*mid vertical edge 1*/
			penta_node_ids[ 7]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+1]+1; /*mid vertical edge 2*/
			penta_node_ids[ 8]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+2]+1; /*mid vertical edge 3*/
			penta_node_ids[ 9]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+0]+2; /* 1/4 vertical edge 1*/
			penta_node_ids[10]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+1]+2; /* 1/4 vertical edge 2*/
			penta_node_ids[11]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+2]+2; /* 1/4 vertical edge 3*/
			penta_node_ids[12]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+0]+3; /* 3/4 vertical edge 1*/
			penta_node_ids[13]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+1]+3; /* 3/4 vertical edge 2*/
			penta_node_ids[14]=iomodel->numberofvertices+3*iomodel->elementtoverticaledgeconnectivity[3*index+2]+3; /* 3/4 vertical edge 3*/
			break;
		case P2xP4Enum:
			numnodes         = 30;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0]; /*Vertex 1*/
			penta_node_ids[ 1]=iomodel->elements[6*index+1]; /*Vertex 2*/
			penta_node_ids[ 2]=iomodel->elements[6*index+2]; /*Vertex 3*/
			penta_node_ids[ 3]=iomodel->elements[6*index+3]; /*Vertex 4*/
			penta_node_ids[ 4]=iomodel->elements[6*index+4]; /*Vertex 5*/
			penta_node_ids[ 5]=iomodel->elements[6*index+5]; /*Vertex 6*/
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1; /*mid vertical edge 1*/
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1; /*mid vertical edge 2*/
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1; /*mid vertical edge 3*/
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1; /*mid basal edge 1*/
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1; /*mid basal edge 2*/
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1; /*mid basal edge 3*/
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1; /*mid top edge 1*/
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1; /*mid top edge 2*/
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1; /*mid top edge 3*/
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+1; /* 1/4 vertical edge 1*/
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+1; /* 1/4 vertical edge 2*/
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+1; /* 1/4 vertical edge 3*/
			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+2; /* 3/4 vertical edge 1*/
			penta_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+2; /* 3/4 vertical edge 2*/
			penta_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+2; /* 3/4 vertical edge 3*/
			penta_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+1; /* 1/4 vertical face 1*/
			penta_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+1; /* 1/4 vertical face 2*/
			penta_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+1; /* 1/4 vertical face 3*/
			penta_node_ids[24]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+2; /* 2/4 vertical face 1*/
			penta_node_ids[25]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+2; /* 2/4 vertical face 2*/
			penta_node_ids[26]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+2; /* 2/4 vertical face 3*/
			penta_node_ids[27]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+3; /* 3/4 vertical face 1*/
			penta_node_ids[28]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+3; /* 3/4 vertical face 2*/
			penta_node_ids[29]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+3; /* 3/4 vertical face 3*/
			break;
		case P2Enum:
			numnodes         = 18;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+0]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+1]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+2]+1;
			break;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			numnodes         = 19;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;
			break;
		case P1P1Enum: case P1P1GLSEnum:
			numnodes         = 12;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];

			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elements[6*index+0];
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elements[6*index+1];
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elements[6*index+2];
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elements[6*index+3];
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elements[6*index+4];
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elements[6*index+5];
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes         = 13;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+index+1;

			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+0];
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+1];
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+2];
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+3];
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+4];
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[6*index+5];
			break;
		case TaylorHoodEnum:
			numnodes         = 24;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+0]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+1]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtoverticalfaceconnectivity[3*index+2]+1;

			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+0];
			penta_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+1];
			penta_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+2];
			penta_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+3];
			penta_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+4];
			penta_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofverticalfaces+iomodel->elements[6*index+5];
			break;
		case LATaylorHoodEnum:
			numnodes         = 18;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			break;
		case OneLayerP4zEnum:
			numnodes         = 30+6;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0]; /*Vertex 1*/
			penta_node_ids[ 1]=iomodel->elements[6*index+1]; /*Vertex 2*/
			penta_node_ids[ 2]=iomodel->elements[6*index+2]; /*Vertex 3*/
			penta_node_ids[ 3]=iomodel->elements[6*index+3]; /*Vertex 4*/
			penta_node_ids[ 4]=iomodel->elements[6*index+4]; /*Vertex 5*/
			penta_node_ids[ 5]=iomodel->elements[6*index+5]; /*Vertex 6*/
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1; /*mid vertical edge 1*/
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1; /*mid vertical edge 2*/
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1; /*mid vertical edge 3*/
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1; /*mid basal edge 1*/
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1; /*mid basal edge 2*/
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1; /*mid basal edge 3*/
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1; /*mid top edge 1*/
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1; /*mid top edge 2*/
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1; /*mid top edge 3*/
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+1; /* 1/4 vertical edge 1*/
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+1; /* 1/4 vertical edge 2*/
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+1; /* 1/4 vertical edge 3*/
			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+0]+2; /* 3/4 vertical edge 1*/
			penta_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+1]+2; /* 3/4 vertical edge 2*/
			penta_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->elementtoverticaledgeconnectivity[3*index+2]+2; /* 3/4 vertical edge 3*/
			penta_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+1; /* 1/4 vertical face 1*/
			penta_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+1; /* 1/4 vertical face 2*/
			penta_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+1; /* 1/4 vertical face 3*/
			penta_node_ids[24]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+2; /* 2/4 vertical face 1*/
			penta_node_ids[25]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+2; /* 2/4 vertical face 2*/
			penta_node_ids[26]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+2; /* 2/4 vertical face 3*/
			penta_node_ids[27]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+0]+3; /* 3/4 vertical face 1*/
			penta_node_ids[28]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+1]+3; /* 3/4 vertical face 2*/
			penta_node_ids[29]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->elementtoverticalfaceconnectivity[3*index+2]+3; /* 3/4 vertical face 3*/

			penta_node_ids[30]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+0];
			penta_node_ids[31]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+1];
			penta_node_ids[32]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+2];
			penta_node_ids[33]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+3];
			penta_node_ids[34]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+4];
			penta_node_ids[35]=iomodel->numberofvertices+iomodel->numberofedges+2*iomodel->numberofverticaledges+3*iomodel->numberofverticalfaces+iomodel->elements[6*index+5];
			break;
		case CrouzeixRaviartEnum:
			numnodes         = 25;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;

			penta_node_ids[19]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+1;
			penta_node_ids[20]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+2;
			penta_node_ids[21]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+3;
			penta_node_ids[22]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+4;
			penta_node_ids[23]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+5;
			penta_node_ids[24]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+iomodel->numberofelements+6*index+6;
			break;
		case LACrouzeixRaviartEnum:
			numnodes         = 19;
			penta_node_ids   = xNew<int>(numnodes);
			penta_node_ids[ 0]=iomodel->elements[6*index+0];
			penta_node_ids[ 1]=iomodel->elements[6*index+1];
			penta_node_ids[ 2]=iomodel->elements[6*index+2];
			penta_node_ids[ 3]=iomodel->elements[6*index+3];
			penta_node_ids[ 4]=iomodel->elements[6*index+4];
			penta_node_ids[ 5]=iomodel->elements[6*index+5];
			penta_node_ids[ 6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+0]+1;
			penta_node_ids[ 7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+1]+1;
			penta_node_ids[ 8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+2]+1;
			penta_node_ids[ 9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+3]+1;
			penta_node_ids[10]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+4]+1;
			penta_node_ids[11]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+5]+1;
			penta_node_ids[12]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+6]+1;
			penta_node_ids[13]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+7]+1;
			penta_node_ids[14]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[9*index+8]+1;
			penta_node_ids[15]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+2]+1;
			penta_node_ids[16]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+3]+1;
			penta_node_ids[17]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elementtofaceconnectivity[5*index+4]+1;
			penta_node_ids[18]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberoffaces+index+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(penta_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(penta_node_ids);

	/*Defaults if not provided in iomodel*/
	switch(analysis_type){

		case StressbalanceAnalysisEnum:
			_assert_(iomodel->Data("md.flowequation.element_equation"));

			if((IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index])))==HOFSApproximationEnum){
				int vertexlids[NUMVERTICES];
				for(i=0;i<NUMVERTICES;i++) vertexlids[i]=iomodel->my_vertices_lids[penta_vertex_ids[i]-1];
				/*Create VzHO and VzFS Enums*/
				if(iomodel->Data("md.initialization.vz") && iomodel->Data("md.flowequation.borderFS")){
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1];
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzFSEnum);
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*(1-iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1]);
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzHOEnum);
				}
				else{
					for(i=0;i<6;i++)nodeinputs[i]=0;
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzFSEnum);
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzHOEnum);
				}
			}
			if((IoCodeToEnumElementEquation(reCast<int>(iomodel->Data("md.flowequation.element_equation")[index])))==SSAFSApproximationEnum){
				int vertexlids[NUMVERTICES];
				for(i=0;i<NUMVERTICES;i++) vertexlids[i]=iomodel->my_vertices_lids[penta_vertex_ids[i]-1];
				/*Create VzSSA and VzFS Enums*/
				if(iomodel->Data("md.initialization.vz") && iomodel->Data("md.flowequation.borderFS")){
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1];
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzFSEnum);
					for(i=0;i<6;i++) nodeinputs[i]=iomodel->Data("md.initialization.vz")[penta_vertex_ids[i]-1]*(1-iomodel->Data("md.flowequation.borderFS")[penta_vertex_ids[i]-1]);
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzSSAEnum);
				}
				else{
					for(i=0;i<6;i++)nodeinputs[i]=0;
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzFSEnum);
					this->SetElementInput(inputs,NUMVERTICES,vertexlids,nodeinputs,VzSSAEnum);
				}
			}
			break;
		default:
			/*No update for other solution types*/
			break;
	}
}
/*}}}*/
void       Penta::UpdateConstraintsExtrudeFromBase(void){/*{{{*/

	if(!IsOnBase()) return;

	int        extrusioninput;
	IssmDouble value,isonbase;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = this->GetInput(extrusioninput);      _assert_(extrusioninput);
	Input* onbase = this->GetInput(MeshVertexonbaseEnum); _assert_(onbase);

	GaussPenta gauss;
	for(int iv=0;iv<this->NumberofNodes(this->element_type);iv++){
		gauss.GaussNode(this->element_type,iv);
		onbase->GetInputValue(&isonbase,&gauss);
		if(isonbase==1.){
			input->GetInputValue(&value,&gauss);
			this->nodes[iv]->ApplyConstraint(0,value);
		}
	}
}/*}}}*/
void       Penta::UpdateConstraintsExtrudeFromTop(void){/*{{{*/

	if(!IsOnSurface()) return;

	int extrusioninput;
	int indices[3]={3,4,5};
	IssmDouble value;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = this->GetInput(extrusioninput); _assert_(extrusioninput);

	GaussPenta gauss;
	for(int i=0;i<3;i++){
		gauss.GaussNode(P1Enum,indices[i]);
		input->GetInputValue(&value,&gauss);
		this->nodes[indices[i]]->ApplyConstraint(0,value);
	}

}
/*}}}*/
int        Penta::UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf){/*{{{*/

	int i;
	int nflipped=0;

	/*Go through nodes, and whoever is on the potential_ungrounding, ends up in nodes_on_iceshelf: */
	for(i=0;i<NUMVERTICES;i++){
		if (reCast<bool,IssmDouble>(vertices_potentially_ungrounding[vertices[i]->Pid()])){
			vec_nodes_on_iceshelf->SetValue(vertices[i]->Pid(),-1.,INS_VAL);

			/*If node was not on ice shelf, we flipped*/
			if(nodes_on_iceshelf[vertices[i]->Pid()]>=0.){
				nflipped++;
			}
		}
	}
	return nflipped;
}
/*}}}*/
void       Penta::ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	PentaRef::GetInputDerivativeValue(dvalue,values,xyz_list,gauss,P1Enum);
}
/*}}}*/
void       Penta::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	PentaRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int        Penta::VelocityInterpolation(void){/*{{{*/
	return PentaRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
int        Penta::VertexConnectivity(int vertexindex){/*{{{*/
	_assert_(this->vertices);
	return this->vertices[vertexindex]->Connectivity();
}
/*}}}*/
void       Penta::VerticalSegmentIndices(int** pindices,int* pnumseg){/*{{{*/

	/*output*/
	int *indices = xNew<int>(3*2);
	indices[0*2 + 0] = 0; indices[0*2 + 1] = 3;
	indices[1*2 + 0] = 1; indices[1*2 + 1] = 4;
	indices[2*2 + 0] = 2; indices[2*2 + 1] = 5;

	/*Assign output pointers*/
	*pindices = indices;
	*pnumseg  = 3;
}
/*}}}*/
void       Penta::VerticalSegmentIndicesBase(int** pindices,int* pnumseg){/*{{{*/

	PentaRef::VerticalSegmentIndicesBase(pindices,pnumseg,this->GetElementType());

	return;
}
/*}}}*/
void       Penta::ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble phi;
	IssmDouble viscosity;
	IssmDouble epsilon[6];

	_assert_(gauss->Enum()==GaussPentaEnum);
	this->StrainRateFS(&epsilon[0],xyz_list,(GaussPenta*)gauss,vx_input,vy_input,vz_input);
	this->material->ViscosityFS(&viscosity,3,xyz_list,(GaussPenta*)gauss,vx_input,vy_input,vz_input);
	GetPhi(&phi,&epsilon[0],viscosity);

	/*Assign output pointer*/
	*pphi = phi;
}
/*}}}*/

#ifdef _HAVE_DAKOTA_
void       Penta::InputScaleFromDakota(IssmDouble* distributed_values, IssmDouble* partition, int npart, int nt, int name){/*{{{*/

	int interp;
	int type;

	/*Branch according to whether we have a transient or not input: */
	type=this->inputs->GetInputObjectEnum(name);
	if(type==PentaInputEnum){
		/*Figure out if we are P0 or P1 interpolation: */
		PentaInput* pentainput = this->inputs->GetPentaInput(name);
		PentaInput* pentainput2 = this->inputs->GetPentaInput(DummyEnum);
		interp=pentainput->GetInterpolation();

		if (interp==P0Enum){
			/*Update the value if this element belongs to the partition: */
			if(partition[this->Sid()]!=-1){
				int lid=this->lid; pentainput->Serve(1,&lid);
				/*scale P0 value  for this element, corresponding to the partition:*/
				IssmDouble value = pentainput->element_values[0];
				value*=distributed_values[(int)partition[this->Sid()]];
				pentainput2->SetInput(P0Enum,this->lid,value);
			}
		}
		else if (interp==P1Enum){
			IssmDouble values[NUMVERTICES];
			int lidlist[NUMVERTICES];
			this->GetVerticesLidList(&lidlist[0]);
			pentainput->Serve(NUMVERTICES,&lidlist[0]);
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=pentainput->element_values[i];
				if(partition[this->vertices[i]->Sid()]!=-1) values[i]*=distributed_values[(int)partition[this->vertices[i]->Sid()]];
			}
			pentainput2->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
		}
		else _error_("Penta::InputScaleFromDakota error message: input interpolation " << EnumToStringx(interp) << " not supported yet!");
	}
	else if(type==TransientInputEnum){

		IssmDouble* steps=NULL;
		int nsteps;
		TransientInput* transientinput = NULL;
		TransientInput* transientinput2 = NULL;

		/*retrieve transient input:*/
		transientinput= this->inputs->GetTransientInput(name);
		transientinput2= this->inputs->GetTransientInput(DummyEnum);

		/*retrieve time steps: */
		transientinput->GetAllTimes(&steps,&nsteps);

		/*double check:*/
		if (nsteps!=nt && nt!=1) _error_("Penta:InputScaleFromDakota error message: transient input " << EnumToStringx(name) <<
				" should have the same number of time steps as the number of time values distributed by Dakota: " << nt << "\n");

		/*needed to update inputs:*/
		int lidlist[NUMVERTICES];
		this->GetVerticesLidList(&lidlist[0]);

		/*go through the transient inputs, and update:*/
		for (int i=0;i<nsteps;i++){
			PentaInput* pentainput=transientinput->GetPentaInput(i);
			PentaInput* pentainput2=transientinput2->GetPentaInput(i);
			interp=pentainput->GetInterpolation();

			if (interp==P0Enum){
				/*Update the value if this element belongs to the partition: */
				if(partition[this->Sid()]!=-1){
					int lid=this->lid; pentainput->Serve(1,&lid);
					/*scale P0 value  for this element, corresponding to the partition:*/
					IssmDouble value = pentainput->element_values[0];
					if(nt==1) value*=distributed_values[(int)partition[this->Sid()]]; //we scale all the time steps  with the same distributed_value
					else value*=distributed_values[(int)partition[this->Sid()]*nsteps+i]; //we scale all the time steps with distributed value for each step

					pentainput2->SetInput(P0Enum,this->lid,value);
				}
			}
			else if (interp==P1Enum){
				IssmDouble values[NUMVERTICES];
				pentainput->Serve(NUMVERTICES,&lidlist[0]);
				for (int j=0;j<NUMVERTICES;j++){
					values[j]=pentainput->element_values[j];
					if(partition[this->vertices[i]->Sid()]!=-1){
						if(nt==1) values[j]*=distributed_values[(int)partition[this->vertices[j]->Sid()]];//we scale all the time steps  with the same distributed_value
						else values[j]*=distributed_values[(int)partition[this->vertices[j]->Sid()]*nsteps+i];//we scale all the time steps with distributed value for each step
					}
				}
				pentainput2->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
			}
			else _error_("Penta::InputScaleFromDakota error message: input interpolation " << EnumToStringx(interp) << " not supported yet!");
		}
	}
	else _error_("Penta::InputScaleFromDakota error message: input type " << EnumToStringx(name) << " not supported yet!");
}
/*}}}*/
void       Penta::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");
	TransientInput* transientinput = inputs->GetTransientInput(name);

	switch(type){

		case VertexEnum:

			/*Get LID lists once for all*/
			IssmDouble  values[NUMVERTICES];
			int         lidlist[NUMVERTICES];
			this->GetVerticesLidList(&lidlist[0]);

			/*Create transient input: */
			for(int t=0;t<ncols;t++){ //ncols is the number of times

				/*create input values: */
				for(int i=0;i<6;i++){
					int row=this->vertices[i]->Sid();
					values[i]=matrix[ncols*row+t];
				}

				/*time:*/
				IssmDouble time=matrix[(nrows-1)*ncols+t];

				transientinput->AddPentaTimeInput(t,NUMVERTICES,&lidlist[0],&values[0],P1Enum);
			}
			break;

		case ElementEnum:
			/*Get value for the element: */
			for(int t=0;t<ncols;t++){ //ncols is the number of times
				IssmDouble value=matrix[ncols*(this->Sid())+t];
				IssmDouble time=matrix[(nrows-1)*ncols+t];
				transientinput->AddPentaTimeInput(t,1,&(this->lid),&value,P0Enum);
			}
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
void       Penta::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/

	int i,j;

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:

			/*New PentaInput*/
			IssmDouble values[6];

			/*Get values on the 6 vertices*/
			for (i=0;i<6;i++){
				values[i]=vector[this->vertices[i]->Sid()]; //careful, vector of values here is not parallel distributed, but serial distributed (from a serial Dakota core!)
			}

			/*Branch on the specified type of update: */
			switch(name){
				case ThicknessEnum:
					/*Update thickness + surface: assume bed is constant. On ice shelves, takes hydrostatic equilibrium*/
					IssmDouble  thickness[6];
					IssmDouble  thickness_init[6];
					IssmDouble  hydrostatic_ratio[6];
					IssmDouble  surface[6];
					IssmDouble  bed[6];

					/*retrieve inputs: */
					Element::GetInputListOnVertices(&thickness_init[0],ThicknessEnum);
					Element::GetInputListOnVertices(&hydrostatic_ratio[0],GeometryHydrostaticRatioEnum);
					Element::GetInputListOnVertices(&bed[0],BaseEnum);
					Element::GetInputListOnVertices(&surface[0],SurfaceEnum);

					/*build new thickness: */
//					for(j=0;j<6;j++)thickness[j]=values[j];

					/*build new bed and surface: */
					if (this->IsAllFloating()){
						/*hydrostatic equilibrium: */
						IssmDouble rho_ice,rho_water,di;
						rho_ice=this->FindParam(MaterialsRhoIceEnum);
						rho_water=this->FindParam(MaterialsRhoSeawaterEnum);

						di=rho_ice/rho_water;

						/*build new thickness: */
						for (j=0; j<6; j++) {
						/*  for observed/interpolated/hydrostatic thickness, remove scaling from any hydrostatic thickness  */
							if     (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j]-(values[j]/thickness_init[j]-1.)*hydrostatic_ratio[j]*surface[j]/(1.-di);
						/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];

						/*  check the computed thickness and update bed  */
							if (thickness[j] < 0.)
								thickness[j]=1./(1.-di);
							bed[j]=surface[j]-thickness[j];
						}

//						for(j=0;j<6;j++){
//							surface[j]=(1-di)*thickness[j];
//							bed[j]=-di*thickness[j];
//						}
					}
					else{
						/*build new thickness: */
						for (j=0; j<6; j++) {
						/*  for observed thickness, use scaled value  */
							if(hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j];
						/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];
						}

						/*update bed on grounded ice: */
//						for(j=0;j<6;j++)surface[j]=bed[j]+thickness[j];
						for(j=0;j<6;j++)bed[j]=surface[j]-thickness[j];
					}

					/*Add new inputs: */
					this->AddInput(ThicknessEnum,thickness,P1Enum);
					this->AddInput(BaseEnum,bed,P1Enum);
					this->AddInput(SurfaceEnum,surface,P1Enum);
					break;

				default:
					this->AddInput(name,values,P1Enum);
			}
			break;

		case ElementEnum:
			IssmDouble value;
			/*Get value for the element: */
			value=vector[this->Sid()]; //careful, vector of values here is not parallel distributed, but serial distributed (from a serial Dakota core!)
			this->AddInput(name,&value,P0Enum);
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

}
/*}}}*/
#endif
