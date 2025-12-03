/*!\file Tria.cpp
 * \brief: implementation of the Tria object
 */
/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <gsl_cblas.h>
#include "../classes.h"
#include "../Inputs/TriaInput.h"
#include "../Inputs/PentaInput.h"
#include "../Inputs/ControlInput.h"
#include "../Inputs/DatasetInput.h"
#include "../Inputs/TransientInput.h"
#include "../../shared/shared.h"
#ifdef _HAVE_SEALEVELCHANGE_
#include "../../modules/GiaDeflectionCorex/GiaDeflectionCorex.h"
#endif
/*}}}*/

/*Element macros*/
#define NUMVERTICES   3
#define NUMVERTICES1D 2
//#define MICI          0 //1 = DeConto & Pollard, 2 = Anna Crawford DOMINOS

/*Constructors/destructor/copy*/
Tria::Tria(int tria_id,int tria_sid,int tria_lid,IoModel* iomodel,int nummodels)/*{{{*/
	:ElementHook(nummodels,tria_id,NUMVERTICES,iomodel){

		this->iscollapsed = 0;

		/*id: */
		this->id  = tria_id;
		this->sid = tria_sid;
		this->lid = tria_lid;

		/*this->parameters: we still can't point to it, it may not even exist. Configure will handle this.*/
		this->parameters = NULL;

		/*initialize pointers:*/
		this->nodes    = NULL;
		this->vertices = NULL;
		this->material = NULL;
		if(nummodels>0){
			this->element_type_list=xNew<int>(nummodels);
			for(int i=0;i<nummodels;i++) this->element_type_list[i] = 0;
		}
		else this->element_type_list = NULL;

		/*surface and base*/
		IssmDouble sum;
		this->isonsurface = false;
		this->isonbase    = false;
		switch(iomodel->domaintype){
			case Domain2DverticalEnum:
				_assert_(iomodel->Data("md.mesh.vertexonsurface"));
				_assert_(iomodel->Data("md.mesh.vertexonbase"));
				sum = 0.;
				for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonsurface")[reCast<int>(iomodel->elements[(tria_id-1)*NUMVERTICES+i])-1];
				_assert_(sum>=0 && sum<3);
				if(sum>1.) this->isonsurface = true;
				sum = 0.;
				for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonbase")[reCast<int>(iomodel->elements[(tria_id-1)*NUMVERTICES+i])-1];
				_assert_(sum>=0 && sum<3);
				if(sum>1.) this->isonbase = true;
				break;
			case Domain2DhorizontalEnum:
			case Domain3DsurfaceEnum:
				this->isonsurface = true;
				this->isonbase    = true;
				break;
			default: _error_("mesh "<<EnumToStringx(iomodel->domaintype)<<" not supported yet");
		}

}/*}}}*/
Tria::~Tria(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Object* Tria::copy() {/*{{{*/

	int i;
	Tria* tria=NULL;

	tria=new Tria();

	tria->iscollapsed=this->iscollapsed;

	//deal with TriaRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		tria->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++){
			if (this->element_type_list[i]) tria->element_type_list[i]=this->element_type_list[i];
			else tria->element_type_list[i] = 0;
		}
	}
	else tria->element_type_list = NULL;
	tria->element_type=this->element_type;
	tria->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		tria->hnodes=xNew<Hook*>(tria->numanalyses);
		for(i=0;i<tria->numanalyses;i++){
			if (this->hnodes[i]) tria->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else tria->hnodes[i] = NULL;
		}
	}
	else tria->hnodes = NULL;

	tria->hvertices = (Hook*)this->hvertices->copy();
	if(this->hmaterial)tria->hmaterial = (Hook*)this->hmaterial->copy();
	tria->hneighbors = NULL;

	/*deal with Tria fields: */
	tria->id  = this->id;
	tria->sid = this->sid;
	tria->lid = this->lid;
	tria->isonbase  = this->isonbase;
	tria->isonsurface  = this->isonsurface;

	/*point parameters: */
	tria->parameters=this->parameters;

	/*recover objects: */
	if (this->nodes){
		unsigned int num_nodes = 3;
		tria->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
		for(i=0;i<num_nodes;i++) if(this->nodes[i]) tria->nodes[i]=this->nodes[i]; else tria->nodes[i] = NULL;
	}
	else tria->nodes = NULL;

	tria->vertices = (Vertex**)this->hvertices->deliverp();
	if(this->hmaterial)tria->material = (Material*)this->hmaterial->delivers();

	return tria;
}
/*}}}*/
void Tria::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = TriaEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->iscollapsed);
	marshallhandle->call(this->isonsurface);
	marshallhandle->call(this->isonbase);

	/*Call parent classes: */
	ElementHook::Marshall(marshallhandle);
	Element::MarshallElement2(marshallhandle,this->numanalyses);
	TriaRef::Marshall(marshallhandle);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();

}
/*}}}*/

/*Other*/
void       Tria::AddBasalInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/*Call inputs method*/
	_assert_(this->inputs);

	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
		case Domain3DsurfaceEnum:
			this->AddInput(input_enum,values,interpolation_enum);
			break;
		case Domain2DverticalEnum:{
			_error_("not implemented yet");
										  }
			break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

}
/*}}}*/
void       Tria::AddInput(int input_enum,IssmDouble* values, int interpolation_enum){/*{{{*/

	/*Intermediaries*/
	int vertexlids[NUMVERTICES];

	/*Call inputs method*/
	if(!this->inputs){
		int* temp = xNew<int>(3);
		_error_("inputs not set");
	}
	_assert_(this->inputs);
	switch(interpolation_enum){
		case P1Enum:
			for(int i=0;i<NUMVERTICES;i++) vertexlids[i]=this->vertices[i]->lid;
			inputs->SetTriaInput(input_enum,interpolation_enum,NUMVERTICES,vertexlids,values);
			break;
		case P1DGEnum:
			for(int i=0;i<NUMVERTICES;i++) vertexlids[i]=this->vertices[i]->lid;
			inputs->SetTriaInput(input_enum,interpolation_enum,this->lid,NUMVERTICES,values);
			break;
		default:
			inputs->SetTriaInput(input_enum,interpolation_enum,this->lid,this->GetNumberOfNodes(interpolation_enum),values);
	}

}
/*}}}*/
void       Tria::AddControlInput(int input_enum,Inputs* inputs,IoModel* iomodel,IssmDouble* values,IssmDouble* values_min,IssmDouble* values_max, int interpolation_enum,int id){/*{{{*/

	/*Intermediaries*/
	int vertexlids[NUMVERTICES];

	_assert_(iomodel->elements);
	for(int i=0;i<NUMVERTICES;i++){
		int vertexid =reCast<int>(iomodel->elements[NUMVERTICES*this->Sid()+i]); //ids for vertices are in the elements array from Matlab
		vertexlids[i]=iomodel->my_vertices_lids[vertexid-1];
	}

	/*Create Control Input*/
	inputs->SetControlInput(input_enum,TriaInputEnum,interpolation_enum,id);
	ControlInput* control_input = inputs->GetControlInput(input_enum); _assert_(input_enum);

	/*Call inputs method*/
	switch(interpolation_enum){
		case P0Enum:
			control_input->SetControl(interpolation_enum,1,&this->lid,values,values_min,values_max);
			break;
		case P1Enum:
			for(int i=0;i<NUMVERTICES;i++){
				if(xIsNan<IssmDouble>(values[i])) values[i] = 0.;
				if(xIsNan<IssmDouble>(values_min[i])) values_min[i] = 0.;
				if(xIsNan<IssmDouble>(values_max[i])) values_max[i] = 0.;
			}
			control_input->SetControl(interpolation_enum,NUMVERTICES,&vertexlids[0],values,values_min,values_max);
			break;
		default:
			_error_("Cannot add \""<<EnumToStringx(input_enum)<<"\" interpolation "<<EnumToStringx(interpolation_enum)<<" not supported");
	}

}
/*}}}*/
void       Tria::DatasetInputCreate(IssmDouble* array,int M,int N,int* individual_enums,int num_inputs,Inputs* inputs,IoModel* iomodel,int input_enum){/*{{{*/

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
		inputs->SetTriaDatasetInput(input_enum,individual_enums[i],P1Enum,NUMVERTICES,vertexlids,nodeinputs);
	}
}
/*}}}*/
void       Tria::AverageOntoPartition(Vector<IssmDouble>* partition_contributions,Vector<IssmDouble>* partition_areas,IssmDouble* vertex_response,IssmDouble* qmu_part){/*{{{*/

	bool       already = false;
	int        i,j;
	int        partition[NUMVERTICES];
	int        offsetsid[NUMVERTICES];
	int        offsetdof[NUMVERTICES];
	IssmDouble area;
	IssmDouble mean;

	/*First, get the area: */
	area=this->GetArea();

	/*Figure out the average for this element: */
	this->GetVerticesSidList(&offsetsid[0]);
	this->GetVerticesPidList(&offsetdof[0]);
	mean=0;
	for(i=0;i<NUMVERTICES;i++){
		partition[i]=reCast<int>(qmu_part[offsetsid[i]]);
		mean=mean+1.0/NUMVERTICES*vertex_response[offsetdof[i]];
	}

	/*Add contribution: */
	for(i=0;i<NUMVERTICES;i++){
		already=false;
		for(j=0;j<i;j++){
			if (partition[i]==partition[j]){
				already=true;
				break;
			}
		}
		if(!already){
			partition_contributions->SetValue(partition[i],mean*area,ADD_VAL);
			partition_areas->SetValue(partition[i],area,ADD_VAL);
		};
	}
}
/*}}}*/
bool       Tria::Buttressing(IssmDouble* ptheta, IssmDouble* plength){/*{{{*/

	/*Make sure there is a grounding line here*/
	if(!IsIceInElement()) return false;
	if(!IsZeroLevelset(MaskOceanLevelsetEnum)) return true;

	int               domaintype,index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2,s_xx,s_yy,s_xy;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		_error_("not implemented");
	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],&xyz_front[0][0]);

	/*Get deviatoric stress tensor*/
	this->ComputeDeviatoricStressTensor();

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice        = this->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_seawater   = this->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble constant_g     = this->FindParam(ConstantsGEnum);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	}
	Input *s_xx_input      = this->GetInput(DeviatoricStressxxEnum); _assert_(s_xx_input);
	Input *s_xy_input      = this->GetInput(DeviatoricStressxyEnum); _assert_(s_xy_input);
	Input *s_yy_input      = this->GetInput(DeviatoricStressyyEnum); _assert_(s_yy_input);
	Input *thickness_input = this->GetInput(ThicknessEnum);          _assert_(thickness_input);

   IssmDouble theta  = 0.;
   IssmDouble length = 0.;
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		s_xx_input->GetInputValue(&s_xx,gauss);
		s_xy_input->GetInputValue(&s_xy,gauss);
		s_yy_input->GetInputValue(&s_yy,gauss);
      this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

      IssmDouble R[2][2];
      R[0][0] = 2*s_xx+s_yy;
      R[1][0] = s_xy;
      R[0][1] = s_xy;
      R[1][1] = 2*s_yy+s_xx;

      IssmDouble N  = normal[0]*(R[0][0]*normal[0]+R[0][1]*normal[1]) + normal[1]*(R[1][0]*normal[0]+R[1][1]*normal[1]);
      IssmDouble N0 = 0.5*constant_g*rho_ice*(1-rho_ice/rho_seawater)*thickness;

		theta  += Jdet*gauss->weight*N/N0;
      length += Jdet*gauss->weight;
	}

	/*Cleanup and return*/
   delete gauss;
	//_assert_(theta>0.);
	_assert_(length>0.);
	*ptheta  = theta;
	*plength = length;
	return true;
}
/*}}}*/
void       Tria::CalvingRateVonmises(){/*{{{*/

	/*First, compute Von Mises Stress*/
	this->ComputeSigmaVM();

	/*Now compute calving rate*/
	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  sigma_vm,vx,vy;
	IssmDouble  sigma_max,sigma_max_floating,sigma_max_grounded,n;
	IssmDouble  groundedice,bed,sealevel;

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input       = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input       = this->GetInput(VyEnum); _assert_(vy_input);
	Input* gr_input       = this->GetInput(MaskOceanLevelsetEnum); _assert_(gr_input);
	Input* bs_input       = this->GetInput(BaseEnum);                    _assert_(bs_input);
	Input* smax_fl_input  = this->GetInput(CalvingStressThresholdFloatingiceEnum); _assert_(smax_fl_input);
	Input* smax_gr_input  = this->GetInput(CalvingStressThresholdGroundediceEnum); _assert_(smax_gr_input);
	Input* sl_input       = this->GetInput(SealevelEnum); _assert_(sl_input);
	Input* sigma_vm_input = this->GetInput(SigmaVMEnum); _assert_(sigma_vm_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components and thickness*/
		sigma_vm_input->GetInputValue(&sigma_vm,&gauss);
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
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
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
   this->CalvingRateToVector();
}
/*}}}*/
void       Tria::CalvingRateVonmisesAD(){/*{{{*/

	/*First, compute Von Mises Stress*/
	this->ComputeSigmaVM();

	/*Now compute calving rate*/
	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  sigma_vm,vx,vy;
	IssmDouble  sigma_max,sigma_max_floating,sigma_max_grounded,n;
	IssmDouble  groundedice,bed,sealevel;
	int M;
	int basinid;
	IssmDouble* sigma_max_floating_basin=NULL;
	IssmDouble* sigma_max_grounded_basin=NULL;

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input       = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input       = this->GetInput(VyEnum); _assert_(vy_input);
	Input* gr_input       = this->GetInput(MaskOceanLevelsetEnum); _assert_(gr_input);
	Input* bs_input       = this->GetInput(BaseEnum);                    _assert_(bs_input);
	Input* sl_input       = this->GetInput(SealevelEnum); _assert_(sl_input);
	Input* sigma_vm_input = this->GetInput(SigmaVMEnum); _assert_(sigma_vm_input);

	this->Element::GetInputValue(&basinid,CalvingBasinIdEnum);

	parameters->FindParam(&sigma_max_floating_basin,&M,CalvingADStressThresholdFloatingiceEnum);
	parameters->FindParam(&sigma_max_grounded_basin,&M,CalvingADStressThresholdGroundediceEnum);

	sigma_max_floating = sigma_max_floating_basin[basinid];
	sigma_max_grounded = sigma_max_grounded_basin[basinid];

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components and thickness*/
		sigma_vm_input->GetInputValue(&sigma_vm,&gauss);
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		gr_input->GetInputValue(&groundedice,&gauss);
		bs_input->GetInputValue(&bed,&gauss);
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
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
   this->CalvingRateToVector();
}
/*}}}*/
void       Tria::CalvingRateTest(){/*{{{*/

	IssmDouble  calvingratex[NUMVERTICES];
	IssmDouble  calvingratey[NUMVERTICES];
	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  vx,vy,vel;
	IssmDouble  dphidx, dphidy, dphi;
	IssmDouble  time;
	IssmDouble  coeff, indrate;
	IssmDouble  bed, bedrate = 1.0;

	/*Retrieve all inputs and parameters we will need*/
	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&coeff,CalvingTestSpeedfactorEnum,time);
	parameters->FindParam(&indrate,CalvingTestIndependentRateEnum,time);

	Input *bs_input = this->GetInput(BedEnum);_assert_(bs_input);
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input *lsf_slopex_input  = this->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
   Input *lsf_slopey_input  = this->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Get velocity components and thickness*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);

		lsf_slopex_input->GetInputValue(&dphidx,&gauss);
      lsf_slopey_input->GetInputValue(&dphidy,&gauss);

		bs_input->GetInputValue(&bed,&gauss);
		bedrate = (bed>0)?0.0:1.0;

      vel=sqrt(vx*vx + vy*vy) + 1e-14;
      dphi=sqrt(dphidx*dphidx+dphidy*dphidy)+ 1e-14;

		calvingratex[iv]= coeff*vx + bedrate*indrate*dphidx/dphi;
		calvingratey[iv]= coeff*vy + bedrate*indrate*dphidy/dphi;
		calvingrate[iv] = sqrt(calvingratex[iv]*calvingratex[iv] + calvingratey[iv]*calvingratey[iv]);
	}

	/*Add input*/
	this->AddInput(CalvingratexEnum,&calvingratex[0],P1DGEnum);
	this->AddInput(CalvingrateyEnum,&calvingratey[0],P1DGEnum);
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
}
/*}}}*/
void       Tria::CalvingCrevasseDepth(){/*{{{*/

	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  vx,vy;
	IssmDouble  water_height, bed,Hab,thickness,surface;
	IssmDouble  surface_crevasse[NUMVERTICES], basal_crevasse[NUMVERTICES], crevasse_depth[NUMVERTICES], H_surf, H_surfbasal;
	IssmDouble  strainparallel, straineffective,B,n;
	IssmDouble  s_xx,s_xy,s_yy,s1,s2,stmp,vH,Kmax;
	int         crevasse_opening_stress;

	/*reset if no ice in element*/
	if(!this->IsIceInElement()){
		for(int i=0;i<NUMVERTICES;i++){
			surface_crevasse[i] = 0.;
			basal_crevasse[i] = 0.;
			crevasse_depth[i] = 0.;
		}
		this->AddInput(SurfaceCrevasseEnum,&surface_crevasse[0],P1DGEnum);
		this->AddInput(BasalCrevasseEnum,&basal_crevasse[0],P1DGEnum);
		this->AddInput(CrevasseDepthEnum,&crevasse_depth[0],P1DGEnum);
		return;
	}

	/*retrieve the type of crevasse_opening_stress*/
	this->parameters->FindParam(&crevasse_opening_stress,CalvingCrevasseDepthEnum);

	IssmDouble rho_ice        = this->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_seawater   = this->FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble rho_freshwater = this->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble constant_g     = this->FindParam(ConstantsGEnum);

	Input*   H_input                 = this->GetInput(ThicknessEnum); _assert_(H_input);
	Input*   bed_input               = this->GetInput(BedEnum); _assert_(bed_input);
	Input*   surface_input           = this->GetInput(SurfaceEnum); _assert_(surface_input);
	Input*	strainrateparallel_input  = this->GetInput(StrainRateparallelEnum);  _assert_(strainrateparallel_input);
	Input*	strainrateeffective_input = this->GetInput(StrainRateeffectiveEnum); _assert_(strainrateeffective_input);
	Input*	vx_input                  = this->GetInput(VxEnum); _assert_(vx_input);
	Input*	vy_input                  = this->GetInput(VxEnum); _assert_(vy_input);
	Input*   waterheight_input       = this->GetInput(WaterheightEnum); _assert_(waterheight_input);
	Input*   s_xx_input              = this->GetInput(DeviatoricStressxxEnum);     _assert_(s_xx_input);
	Input*   s_xy_input              = this->GetInput(DeviatoricStressxyEnum);     _assert_(s_xy_input);
	Input*   s_yy_input              = this->GetInput(DeviatoricStressyyEnum);     _assert_(s_yy_input);
	Input*	B_input  = this->GetInput(MaterialsRheologyBbarEnum);   _assert_(B_input);
	Input*	n_input  = this->GetInput(MaterialsRheologyNEnum);   _assert_(n_input);

	/*Loop over all elements of this partition*/
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		H_input->GetInputValue(&thickness,&gauss);
		bed_input->GetInputValue(&bed,&gauss);
		surface_input->GetInputValue(&surface,&gauss);

		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		waterheight_input->GetInputValue(&water_height,&gauss);
		s_xx_input->GetInputValue(&s_xx,&gauss);
		s_xy_input->GetInputValue(&s_xy,&gauss);
		s_yy_input->GetInputValue(&s_yy,&gauss);

		/*Get longitudinal or maximum Eigen stress*/
		if(crevasse_opening_stress==0){
			/*Otero2010: balance between the tensile deviatoric stress and ice overburden pressure*/
			strainrateparallel_input->GetInputValue(&strainparallel,&gauss);
			strainrateeffective_input->GetInputValue(&straineffective,&gauss);
			B_input->GetInputValue(&B,&gauss);
			n_input->GetInputValue(&n,&gauss);
			s1 =  B * strainparallel * pow(straineffective, (1./n)-1);
		}
		else if(crevasse_opening_stress==1){
			/*Benn2017,Todd2018: maximum principal stress */
			Matrix2x2Eigen(&s1,&s2,NULL,NULL,s_xx,s_xy,s_yy);
			s1 = max(0.,max(s1,s2));
		}
		else if(crevasse_opening_stress==2){
			/*Coffey 2024, Buttressing based */
			Matrix2x2Eigen(&s1,&s2,NULL,NULL,s_xx,s_xy,s_yy);
			strainrateeffective_input->GetInputValue(&straineffective,&gauss);
			n_input->GetInputValue(&n,&gauss);
			B_input->GetInputValue(&B,&gauss);
			if((straineffective <= 0.) || (thickness <= 0.) ){
				vH = 1e14;
			}
			else{
				vH = 0.5*B/thickness*pow(straineffective, (1./n)-1);
			}
			Kmax = 1.0 - 4.0*vH*(s1+s2+min(s1,s2))/(rho_ice*constant_g*(rho_seawater-rho_ice)/rho_seawater);
			if(Kmax<0.) Kmax = 0.0;
		}
		else{
			_error_("not supported");
		}

		if(crevasse_opening_stress==2) {
			/*Coffey 2024, Buttressing based */
			surface_crevasse[iv] = thickness*(1.0-rho_ice/rho_seawater)*(1.0 - sqrt(Kmax));
			basal_crevasse[iv]   = thickness*(rho_ice/rho_seawater)*(1.0 - sqrt(Kmax));
			//_printf0_(Kmax<<", "<<basal_crevasse[iv]<<", "<<surface_crevasse[iv]<<endl);
		}
		else {
			/*Surface crevasse: sigma'_xx - rho_i g d + rho_fw g d_w = 0*/
			surface_crevasse[iv] = 2*s1 / (rho_ice*constant_g) + (rho_freshwater/rho_ice)*water_height;

			/*Basal crevasse: sigma'_xx - rho_i g (H-d) - rho_w g (b+d) = 0*/
			if(bed>0.){
				basal_crevasse[iv] = 0.;
			}
			else{
				Hab = thickness - (rho_seawater/rho_ice) * (-bed);
				if(Hab<0.)  Hab=0.;
				basal_crevasse[iv] = (rho_ice/(rho_seawater-rho_ice))* (2*s1/ (rho_ice*constant_g)-Hab);
			}
		}

		/*crevasse depths (total = surface + basal)*/
		if(surface_crevasse[iv]<0.) surface_crevasse[iv]=0.;
		if(basal_crevasse[iv]<0.)   basal_crevasse[iv]=0.;
		crevasse_depth[iv]  = surface_crevasse[iv] + basal_crevasse[iv];
	}

	this->AddInput(SurfaceCrevasseEnum,&surface_crevasse[0],P1DGEnum);
	this->AddInput(BasalCrevasseEnum,&basal_crevasse[0],P1DGEnum);
	this->AddInput(CrevasseDepthEnum,&crevasse_depth[0],P1DGEnum);
}
/*}}}*/
void       Tria::CalvingRateLevermann(){/*{{{*/

	IssmDouble  strainparallel;
	IssmDouble  propcoeff,bed;
	IssmDouble  strainperpendicular;
	IssmDouble  calvingrate[NUMVERTICES];

	/*Retrieve all inputs and parameters we will need*/
	Input *bs_input                  = this->GetInput(BaseEnum);                    _assert_(bs_input);
	Input *strainparallel_input      = this->GetInput(StrainRateparallelEnum);      _assert_(strainparallel_input);
	Input *strainperpendicular_input = this->GetInput(StrainRateperpendicularEnum); _assert_(strainperpendicular_input);
	Input *levermanncoeff_input      = this->GetInput(CalvinglevermannCoeffEnum);   _assert_(levermanncoeff_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		strainparallel_input->GetInputValue(&strainparallel,&gauss);
		strainperpendicular_input->GetInputValue(&strainperpendicular,&gauss);
		levermanncoeff_input->GetInputValue(&propcoeff,&gauss);
		bs_input->GetInputValue(&bed,&gauss);

		/*Calving rate proportionnal to the positive product of the strain rate along the ice flow direction and the strain rate perpendicular to the ice flow */
		if(strainparallel>0. && strainperpendicular>0. && bed<=0.){
			calvingrate[iv]=propcoeff*strainparallel*strainperpendicular;
		}
		else
			calvingrate[iv]=0.;
	}

	/*Add input*/
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->CalvingRateToVector();
}/*}}}*/
void       Tria::CalvingPollard(){/*{{{*/

	/*Intermediaries*/
	IssmDouble calvingrate[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble dvx[2], dvy[2];
	IssmDouble B, n, H, bed, vx, vy, vel, smb;
	IssmDouble ds, db, da, dt, dw, r, R;

	if(!IsIceInElement()){
		for (int iv=0;iv<NUMVERTICES;iv++) calvingrate[iv]=0.;
		this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
		this->CalvingRateToVector();
		return;
	}

	/*Retrieve all inputs and parameters we will need*/
	IssmDouble rc        = FindParam(CalvingRcEnum);
	IssmDouble rho_ice   = FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble gravity   = FindParam(ConstantsGEnum);
   IssmDouble mig_max   = FindParam(MigrationMaxEnum);

	/*Retrieve all inputs and parameters we will need */
	Input *bs_input  = this->GetInput(BaseEnum);                  _assert_(bs_input);
	Input *vx_input  = this->GetInput(VxEnum);                    _assert_(vx_input);
	Input *vy_input  = this->GetInput(VyEnum);                    _assert_(vy_input);
	Input *B_input   = this->GetInput(MaterialsRheologyBbarEnum); _assert_(B_input);
	Input *n_input   = this->GetInput(MaterialsRheologyNEnum);    _assert_(n_input);
	Input *H_input   = this->GetInput(ThicknessEnum);             _assert_(H_input);
	Input *smb_input = this->GetInput(SmbMassBalanceEnum);        _assert_(smb_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		bs_input->GetInputValue(&bed,&gauss);

		/*Only calve if bed is below sea level, as always*/
		if(bed<=0.){

			/*Get Triangle node coordinates*/
			::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

			/*Get strain rates*/
			vx_input->GetInputDerivativeValue(&dvx[0],&xyz_list[0][0],&gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],&xyz_list[0][0],&gauss);

			/*Get other inputs*/
			B_input->GetInputValue(&B,&gauss);
         n_input->GetInputValue(&n,&gauss);
         H_input->GetInputValue(&H,&gauss);
         vx_input->GetInputValue(&vx,&gauss);
         vy_input->GetInputValue(&vy,&gauss);
         smb_input->GetInputValue(&smb,&gauss);

			/*1. with surface crevasses, ds*/
			ds = 2./(rho_ice*gravity) * B * pow( max(0.,dvx[0]) + max(0.,dvy[1]) , 1./n);

         /*2. basal crevasses*/
         db = (rho_ice)/(rho_water - rho_ice) * ds;

         /*3. "Additional" crevasse opening*/
			vel = sqrt(vx*vx + vy*vy);
			da = H* max(0., log(vel*365*24*3600/1600.))/log(1.2);

         /*4. deal with shallow ice*/
         dt = H* max(0., min(1., (150. - H)/50.));

         /*5. water induced opening*/
         dw = 0.;
         R = smb*365.25*24*3600; //convert from m/s to m/yr
			if(R>1.5 && R<=3.){
            dw = 4*1.5*(R - 1.5);
         }
         else if(R>3.){
            dw = R*R;
         }

         /*Total calving rate*/
         r = (ds+db+da+dt+dw)/H;
			//if(this->Id()==1){
			//	printf("rc = %g\n",rc);
			//	printf("ds = %g\n",ds);
			//}
			calvingrate[iv]= mig_max * max(0., min(1., (r - rc)/(1 - rc))); //P&DC: mig_max = 3000 m/yr
			_assert_(!xIsNan<IssmDouble>(calvingrate[iv]));
			_assert_(!xIsInf<IssmDouble>(calvingrate[iv]));
		}
		else
		 calvingrate[iv]=0.;
	}

	/*Add input*/
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->CalvingRateToVector();
}/*}}}*/
void       Tria::CalvingFluxLevelset(){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)){
		IssmDouble flux_per_area=0;
		this->AddInput(CalvingFluxLevelsetEnum,&flux_per_area,P0Enum);
	}
	else{
		int               domaintype,index1,index2;
		const IssmPDouble epsilon = 1.e-15;
		IssmDouble        s1,s2;
		IssmDouble        gl[NUMVERTICES];
		IssmDouble        xyz_front[2][3];

		IssmDouble  xyz_list[NUMVERTICES][3];
		::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

		/*Recover parameters and values*/
		parameters->FindParam(&domaintype,DomainTypeEnum);
		Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

		/*Be sure that values are not zero*/
		if(gl[0]==0.) gl[0]=gl[0]+epsilon;
		if(gl[1]==0.) gl[1]=gl[1]+epsilon;
		if(gl[2]==0.) gl[2]=gl[2]+epsilon;

		if(domaintype==Domain2DverticalEnum){
			_error_("not implemented");
		}
		else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

		}
		else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

		/*Some checks in debugging mode*/
		_assert_(s1>=0 && s1<=1.);
		_assert_(s2>=0 && s2<=1.);

		/*Get normal vector*/
		IssmDouble normal[3];
		this->NormalSection(&normal[0],&xyz_front[0][0]);
		normal[0] = -normal[0];
		normal[1] = -normal[1];

		/*Get inputs*/
		IssmDouble flux = 0.;
		IssmDouble area = 0.;
		IssmDouble calvingratex,calvingratey,thickness,Jdet,flux_per_area;
		IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
		Input* thickness_input    = this->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* calvingratex_input = this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
		Input* calvingratey_input = this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);

		/*Start looping on Gaussian points*/
		Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
		while(gauss->next()){
			thickness_input->GetInputValue(&thickness,gauss);
			calvingratex_input->GetInputValue(&calvingratex,gauss);
			calvingratey_input->GetInputValue(&calvingratey,gauss);
			this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

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
void       Tria::CalvingMeltingFluxLevelset(){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)){
		IssmDouble flux_per_area=0;
		this->AddInput(CalvingMeltingFluxLevelsetEnum,&flux_per_area,P0Enum);
	}
	else{
		int               domaintype,index1,index2;
		const IssmPDouble epsilon = 1.e-15;
		IssmDouble        s1,s2;
		IssmDouble        gl[NUMVERTICES];
		IssmDouble        xyz_front[2][3];

		IssmDouble  xyz_list[NUMVERTICES][3];
      ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

		/*Recover parameters and values*/
		parameters->FindParam(&domaintype,DomainTypeEnum);
		Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

		/*Be sure that values are not zero*/
		if(gl[0]==0.) gl[0]=gl[0]+epsilon;
		if(gl[1]==0.) gl[1]=gl[1]+epsilon;
		if(gl[2]==0.) gl[2]=gl[2]+epsilon;

		if(domaintype==Domain2DverticalEnum){
			_error_("not implemented");
		}
		else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

		}
		else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

		/*Some checks in debugging mode*/
		_assert_(s1>=0 && s1<=1.);
		_assert_(s2>=0 && s2<=1.);

		/*Get normal vector*/
		IssmDouble normal[3];
		this->NormalSection(&normal[0],&xyz_front[0][0]);
		normal[0] = -normal[0];
		normal[1] = -normal[1];

		/*Get inputs*/
		IssmDouble flux = 0.;
		IssmDouble area = 0.;
		IssmDouble calvingratex,calvingratey,vx,vy,vel,meltingrate,meltingratex,meltingratey,thickness,Jdet,flux_per_area;
		IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
		Input* thickness_input    = this->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* calvingratex_input = this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
		Input* calvingratey_input = this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
		Input* vx_input           = this->GetInput(VxEnum); _assert_(vx_input);
		Input* vy_input           = this->GetInput(VyEnum); _assert_(vy_input);
		Input* meltingrate_input  = this->GetInput(CalvingMeltingrateEnum); _assert_(meltingrate_input);

		/*Start looping on Gaussian points*/
		Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
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
			this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

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
void       Tria::CalvingRateParameterization(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  calvingrate[NUMVERTICES];
	IssmDouble  lambda1,lambda2,ex,ey,vx,vy,vel;
	IssmDouble  sigma_vm[NUMVERTICES];
	IssmDouble  B, n;
	IssmDouble  epse_2,groundedice,bed,sealevel;
	IssmDouble  arate, rho_ice, rho_water, thickness;
	int			use_parameter=-1;
	IssmDouble  gamma, theta, alpha, xoffset, yoffset;
	IssmDouble  vel_lower, vel_upper, vel_threshold, vrate, truncateVrate, vel_max;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs and parameters we will need*/
	Input *vx_input      = this->GetInput(VxEnum);                                _assert_(vx_input);
	Input *vy_input      = this->GetInput(VyEnum);                                _assert_(vy_input);
	Input *B_input       = this->GetInput(MaterialsRheologyBbarEnum);             _assert_(B_input);
	Input *gr_input      = this->GetInput(MaskOceanLevelsetEnum);                 _assert_(gr_input);
	Input *bs_input      = this->GetInput(BedEnum);                               _assert_(bs_input);
	Input *H_input       = this->GetInput(ThicknessEnum);                         _assert_(H_input);
	Input *n_input       = this->GetInput(MaterialsRheologyNEnum);                _assert_(n_input);
	Input *sl_input      = this->GetInput(SealevelEnum);                          _assert_(sl_input);
	Input *arate_input   = this->GetInput(CalvingAblationrateEnum);               _assert_(arate_input);

	/* Ice and sea water density */
	this->FindParam(&rho_ice,MaterialsRhoIceEnum);
	this->FindParam(&rho_water,MaterialsRhoSeawaterEnum);

	/* Use which parameter  */
	this->FindParam(&use_parameter, CalvingUseParamEnum);
	this->FindParam(&theta, CalvingThetaEnum);
	this->FindParam(&alpha, CalvingAlphaEnum);
	this->FindParam(&xoffset, CalvingXoffsetEnum);
	this->FindParam(&yoffset, CalvingYoffsetEnum);
	this->FindParam(&vel_lower, CalvingVelLowerboundEnum);
	this->FindParam(&vel_upper, CalvingVelUpperboundEnum);
	this->FindParam(&vel_threshold, CalvingVelThresholdEnum);
	this->FindParam(&vel_max, CalvingVelMaxEnum);

	/* Start looping on the number of vertices: */
	GaussTria* gauss=new GaussTria();
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss->GaussVertex(iv);

		/*Get velocity components and thickness*/
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		gr_input->GetInputValue(&groundedice,gauss);
		bs_input->GetInputValue(&bed,gauss);
		H_input->GetInputValue(&thickness,gauss);
		vel=sqrt(vx*vx+vy*vy)+1.e-14;
		sl_input->GetInputValue(&sealevel,gauss);
		arate_input->GetInputValue(&arate,gauss);
		/* reduce the arate by a factor depends on a fixed velocity threshold */
		vrate = 1.0;
		if (vel < vel_threshold) vrate = vel / vel_threshold;
		/* parameter used in the calving law: real time max vel, or predefined vel range*/
		if (vel_upper <= vel_lower) {
			truncateVrate = vel / vel_max;
		}
		else {
			truncateVrate = (min(vel_upper, max(vel_lower, vel))-vel_lower) / (vel_upper - vel_lower);
		}

		/*Compute strain rate and viscosity: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],gauss,vx_input,vy_input);

		/*Get Eigen values*/
		Matrix2x2Eigen(&lambda1,&lambda2,&ex,&ey,epsilon[0],epsilon[2],epsilon[1]);
		_assert_(!xIsNan<IssmDouble>(lambda1));
		_assert_(!xIsNan<IssmDouble>(lambda2));

		/*Process Eigen values (only account for extension)*/
		lambda1 = max(lambda1,0.);
		lambda2 = max(lambda2,0.);

		/*Calculate sigma_vm*/
		epse_2    = 1./2. *(lambda1*lambda1 + lambda2*lambda2);
		sigma_vm[iv]  = sqrt(3.) * B * pow(epse_2,1./(2.*n));

		switch (use_parameter) { 
			case 0:
				/* 0 Linear: f(x) = y_{o} + \alpha (x+x_{o}) */
				gamma = yoffset + alpha * (bed+xoffset);
				break;
			case 1:
				/* 1 tanh: f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})) */
				gamma = yoffset -  0.5*theta*tanh(alpha*(bed+xoffset));
				break;
			case 2:
				/* 2 tanh(thicknes): f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})) */
				gamma = yoffset -  0.5*theta*tanh(alpha*(-thickness+xoffset));
				break;
			case 3:
				/* 3 tanh(normal vel): f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})) */
				gamma = yoffset -  0.5*theta*tanh(alpha*(truncateVrate+xoffset));
				break;
			case 4:
				/* 4 tanh(normal vel), cross (0,0): f(x)=y_{o}-\frac{\theta}{2}\tanh(\alpha(x+x_{o})) */
				gamma = 0.5*theta*(tanh(alpha*xoffset) - tanh(alpha*(truncateVrate+xoffset)));
				break;
			case 5:
				/* 5 Linear: f(x) = y_{o} + \alpha (x+x_{o}) */
				gamma = yoffset + alpha * (truncateVrate+xoffset);
				break;
			case -1:
				/* nothing, just the arate*/
				gamma = 1;
				break;
			default:
				_error_("The parameter is not supported yet!");
		}

		/* set upper and lower bounds */
		if (gamma > 1.0) gamma = 1.0;
		if (gamma < 0.0) gamma = 0.0;
		if (bed >= sealevel) gamma = 0.0;

		/*-------------------------------------------*/
		calvingrate[iv] = arate*gamma*vrate;
	}
	/*Add input*/
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	this->AddInput(SigmaVMEnum,&sigma_vm[0],P1DGEnum);
	this->CalvingRateToVector();

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void       Tria::CalvingRateCalvingMIP(){/*{{{*/

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

	/*Retrieve all inputs and parameters we will need*/
	Input *vx_input      = this->GetInput(VxEnum);                                _assert_(vx_input);
	Input *vy_input      = this->GetInput(VyEnum);                                _assert_(vy_input);
	Input *wrate_input   = this->GetInput(CalvingAblationrateEnum);               _assert_(wrate_input); 
	Input* gr_input      = this->GetInput(MaskOceanLevelsetEnum);						_assert_(gr_input);

	/* Use which experiment: use existing Enum */
	this->FindParam(&experiment, CalvingUseParamEnum);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
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
				/* Exp 2: set c=v-wrate(given)*/
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
	this->AddInput(CalvingCalvingrateEnum,&calvingrate[0],P1DGEnum);
	switch (experiment) {
		case 1:
		case 3:
			this->CalvingRateToVector(true);
			break;
		case 2:
		case 4:
			this->CalvingRateToVector(false);
			break;
		default:
			_error_("The experiment is not supported yet!");
	}
}
/*}}}*/
IssmDouble Tria::CharacteristicLength(void){/*{{{*/

	return sqrt(2*this->GetArea());
}
/*}}}*/
void       Tria::ComputeBasalStress(void){/*{{{*/
	_error_("Not Implemented yet");
}
/*}}}*/
void       Tria::ComputeDeviatoricStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  viscosity,lambda1,lambda2;
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  tau_xx[NUMVERTICES];
	IssmDouble	tau_yy[NUMVERTICES];
	IssmDouble	tau_zz[NUMVERTICES]={0,0,0};
	IssmDouble  tau_xy[NUMVERTICES];
	IssmDouble	tau_xz[NUMVERTICES]={0,0,0};
	IssmDouble	tau_yz[NUMVERTICES]={0,0,0};
	IssmDouble  tau_e[NUMVERTICES];
	IssmDouble  tau_1[NUMVERTICES];
	IssmDouble  tau_2[NUMVERTICES];
	int domaintype,dim=2;

	/*Get approximation*/
	int approximation;
	this->Element::GetInputValue(&approximation,ApproximationEnum);

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	this->FindParam(&domaintype,DomainTypeEnum);
	Input* vx_input=this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum); _assert_(vy_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate and viscosity: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);
		switch(approximation){
			case SSAApproximationEnum:
				this->material->ViscositySSA(&viscosity,dim,&xyz_list[0][0],&gauss,vx_input,vy_input);
				break;
			case HOApproximationEnum:
				this->material->ViscosityHO(&viscosity,dim,&xyz_list[0][0],&gauss,vx_input,vy_input);
				break;
			case FSApproximationEnum:
				this->material->ViscosityFS(&viscosity,dim,&xyz_list[0][0],&gauss,vx_input,vy_input,NULL);
				break;
			default:
				_error_("not supported yet");
		}

		/*Compute Stress*/
		tau_xx[iv]=2*viscosity*epsilon[0]; // tau = nu eps
		tau_yy[iv]=2*viscosity*epsilon[1];
		tau_xy[iv]=2*viscosity*epsilon[2];
		tau_e[iv]=1/sqrt(2)*sqrt(pow(tau_xx[iv],2)+pow(tau_yy[iv],2)+2*pow(tau_xy[iv],2));

		/*Get Eigen values*/
		Matrix2x2Eigen(&tau_2[iv],&tau_1[iv],NULL,NULL,tau_xx[iv],tau_xy[iv],tau_yy[iv]);
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(DeviatoricStressxxEnum,&tau_xx[0],P1DGEnum);
	this->AddInput(DeviatoricStressxyEnum,&tau_xy[0],P1DGEnum);
	this->AddInput(DeviatoricStressxzEnum,&tau_xz[0],P1DGEnum);
	this->AddInput(DeviatoricStressyyEnum,&tau_yy[0],P1DGEnum);
	this->AddInput(DeviatoricStressyzEnum,&tau_yz[0],P1DGEnum);
	this->AddInput(DeviatoricStresszzEnum,&tau_zz[0],P1DGEnum);
	this->AddInput(DeviatoricStresseffectiveEnum,&tau_e[0],P1DGEnum);
	this->AddInput(DeviatoricStress1Enum,&tau_1[0],P1DGEnum);
	this->AddInput(DeviatoricStress2Enum,&tau_2[0],P1DGEnum);
}
/*}}}*/
void	      Tria::ComputeEsaStrainAndVorticity(){ /*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[4]; /* epsilon=[exx,eyy,exy+ (shear),exy- (rotation)];*/
	IssmDouble  strain_xx[NUMVERTICES];
	IssmDouble  strain_yy[NUMVERTICES];
	IssmDouble  strain_xy[NUMVERTICES];
	IssmDouble  vorticity_xy[NUMVERTICES];

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input=this->GetInput(EsaXmotionEnum); _assert_(vx_input);
	Input* vy_input=this->GetInput(EsaYmotionEnum); _assert_(vy_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate and vorticity rate: */
		this->StrainRateESA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);

		/*Compute Stress*/
		strain_xx[iv]=epsilon[0];
		strain_yy[iv]=epsilon[1];
		strain_xy[iv]=epsilon[2];
		vorticity_xy[iv]=epsilon[3];
	}

	/*Add Stress tensor components into inputs*/
	this->AddInput(EsaStrainratexxEnum,&strain_xx[0],P1DGEnum);
	this->AddInput(EsaStrainrateyyEnum,&strain_yy[0],P1DGEnum);
	this->AddInput(EsaStrainratexyEnum,&strain_xy[0],P1DGEnum);
	this->AddInput(EsaRotationrateEnum,&vorticity_xy[0],P1DGEnum);
}
/*}}}*/
void       Tria::ComputeSigmaNN(){/*{{{*/

	if(!IsOnBase()){
		IssmDouble sigma_nn[3]={0.};
		this->AddInput(SigmaNNEnum,&sigma_nn[0],P1Enum);
		return;
	}
	else{
		IssmDouble* xyz_list=NULL;
		IssmDouble *xyz_list_base=NULL;
		IssmDouble  pressure,viscosity;
		IssmDouble  sigma_nn[3];
		IssmDouble  sigma_xx,sigma_xy,sigma_yy;
		IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
		IssmDouble  base_normal[2];
		int domaintype,dim=2;

		/* Get node coordinates and dof list: */
		GetVerticesCoordinates(&xyz_list);
	   GetVerticesCoordinatesBase(&xyz_list_base);

		/*Retrieve all inputs we will be needing: */
		this->FindParam(&domaintype,DomainTypeEnum);
		if(domaintype==Domain2DhorizontalEnum) _error_("stress tensor calculation not supported for mesh of type " <<EnumToStringx(domaintype)<<", extrude mesh or call ComputeDeviatoricStressTensor");
		Input* pressure_input=this->GetInput(PressureEnum); _assert_(pressure_input);
		Input* vx_input=this->GetInput(VxEnum);             _assert_(vx_input);
		Input* vy_input=this->GetInput(VyEnum);             _assert_(vy_input);

		/* Start looping on the number of vertices: */
		Gauss* gauss = this->NewGauss();
		for(int i=0;i<NUMVERTICES;i++){
			gauss->GaussNode(P1Enum,i);

			/*Compute strain rate viscosity and pressure: */
			this->StrainRateSSA(&epsilon[0],xyz_list,gauss,vx_input,vy_input);
			this->material->ViscosityFS(&viscosity,dim,xyz_list,gauss,vx_input,vy_input,NULL);
			pressure_input->GetInputValue(&pressure,gauss);

			/*Compute Stress*/
			sigma_xx=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
			sigma_yy=2*viscosity*epsilon[1]-pressure;
			sigma_xy=2*viscosity*epsilon[2];

			/*Get normal vector to the bed */
			NormalBase(&base_normal[0],xyz_list_base);

			/*Compute sigma_nn*/
			sigma_nn[i]=sigma_xx*base_normal[0]*base_normal[0] + 2*sigma_xy*base_normal[0]*base_normal[1] + sigma_yy*base_normal[1]*base_normal[1];
		}

		/*Add Stress tensor components into inputs*/
		this->AddInput(SigmaNNEnum,&sigma_nn[0],P1Enum);

		/*Clean up and return*/
		xDelete<IssmDouble>(xyz_list);
		xDelete<IssmDouble>(xyz_list_base);
		delete gauss;
	}
}
/*}}}*/
void       Tria::ComputeSigmaVM(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  lambda1,lambda2,ex,ey,vx,vy,vel;
	IssmDouble  sigma_vm[NUMVERTICES];
	IssmDouble  B,n;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs and parameters we will need*/
	Input* vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input = this->GetInput(VyEnum); _assert_(vy_input);
	Input* B_input  = this->GetInput(MaterialsRheologyBbarEnum);   _assert_(B_input);
	Input* n_input  = this->GetInput(MaterialsRheologyNEnum); _assert_(n_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
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
		IssmDouble epse_2 = 1./2. *(lambda1*lambda1 + lambda2*lambda2);
		sigma_vm[iv]  = sqrt(3.) * B * pow(epse_2,1./(2.*n));
	}

	/*Add input*/
	this->AddInput(SigmaVMEnum,&sigma_vm[0],P1DGEnum);
}
/*}}}*/
void       Tria::ComputeStressTensor(){/*{{{*/

	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  pressure,viscosity;
	IssmDouble  epsilon[3]; /* epsilon=[exx,eyy,exy];*/
	IssmDouble  sigma_xx[NUMVERTICES];
	IssmDouble	sigma_yy[NUMVERTICES];
	IssmDouble	sigma_zz[NUMVERTICES]={0,0,0};
	IssmDouble  sigma_xy[NUMVERTICES];
	IssmDouble	sigma_xz[NUMVERTICES]={0,0,0};
	IssmDouble	sigma_yz[NUMVERTICES]={0,0,0};
	int domaintype,dim=2;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype==Domain2DhorizontalEnum) _error_("stress tensor calculation not supported for mesh of type " <<EnumToStringx(domaintype)<<", extrude mesh or call ComputeDeviatoricStressTensor");
	Input* pressure_input=this->GetInput(PressureEnum); _assert_(pressure_input);
	Input* vx_input=this->GetInput(VxEnum);             _assert_(vx_input);
	Input* vy_input=this->GetInput(VyEnum);             _assert_(vy_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);
		this->material->ViscositySSA(&viscosity,dim,&xyz_list[0][0],&gauss,vx_input,vy_input);
		pressure_input->GetInputValue(&pressure,&gauss);

		/*Compute Stress*/
		sigma_xx[iv]=2*viscosity*epsilon[0]-pressure; // sigma = nu eps - pressure
		sigma_yy[iv]=2*viscosity*epsilon[1]-pressure;
		sigma_xy[iv]=2*viscosity*epsilon[2];
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
void       Tria::Configure(Elements* elementsin, Loads* loadsin,Nodes* nodesin,Vertices *verticesin,Materials* materialsin, Parameters* parametersin,Inputs* inputsin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	if (this->element_type_list) this->element_type=this->element_type_list[analysis_counter];

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective
	 * datasets, using internal ids and offsets hidden in hooks: */
	if(this->hnodes){
		if (this->hnodes[analysis_counter]) this->hnodes[analysis_counter]->configure(nodesin);
		else this->hnodes[analysis_counter] = NULL;
	}
	else this->hnodes = NULL;
	this->hvertices->configure(verticesin);
	if(this->hmaterial) this->hmaterial->configure(materialsin);

	/*Now, go pick up the objects inside the hooks: */
	if(this->hnodes && this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices = (Vertex**)this->hvertices->deliverp();
	if(this->hmaterial)this->material = (Material*)this->hmaterial->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
	this->inputs=inputsin;
}/*}}}*/
void       Tria::ControlInputSetGradient(IssmDouble* gradient,int control_enum,int control_index,int offset,int M,int N,int interp){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         lidlist[NUMVERTICES];

	/*Get list of ids for this element and this control*/
	int* idlist = xNew<int>(NUMVERTICES*N);
	GradientIndexing(&idlist[0],control_index);

	ControlInput* control_input=this->inputs->GetControlInput(control_enum); _assert_(control_input);
	this->GetVerticesLidList(&lidlist[0]);

	/*Get values on vertices*/
	if(control_input->layout_enum==TriaInputEnum){
		ElementInput* gradient_input = control_input->GetInput("gradient"); _assert_(gradient_input);
		if(gradient_input->GetInputInterpolationType()==P1Enum){
			_assert_(N==1);
			for(int i=0;i<NUMVERTICES;i++) values[i] = gradient[idlist[i]];
			gradient_input->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
		}
		else if(gradient_input->GetInputInterpolationType()==P0Enum){
			_assert_(N==1);
			gradient_input->SetInput(P0Enum,this->lid,gradient[idlist[0]]);
		}
		else{
			_error_("not implemented yet");
		}
	}
	else if(control_input->layout_enum==TransientInputEnum){
		_assert_(N>1);

		int* interp = NULL;
		parameters->FindParam(&interp,NULL,ControlInputInterpolationEnum);

		TransientInput* gradient_input = control_input->GetTransientInput("gradient"); _assert_(gradient_input);

		for(int n=0;n<N;n++){
			if(interp[control_index]==P1Enum){
				for(int i=0;i<NUMVERTICES;i++) values[i] = gradient[idlist[i]];
				gradient_input->AddTriaTimeInput(n,NUMVERTICES,&lidlist[0],&values[0],P1Enum);
			}
			else if(interp[control_index]==P0Enum){
				gradient_input->AddTriaTimeInput(n,1,&(this->lid),&gradient[idlist[n]],P0Enum);
			}
			else{
				_error_("not implemented yet");
			}
		}
		xDelete<int>(interp);
	}
	else _error_("Type not supported");

	/*Clean up*/
	xDelete<int>(idlist);

}/*}}}*/
void       Tria::ControlToVectors(Vector<IssmPDouble>* vector_control, Vector<IssmPDouble>* vector_gradient,int control_enum,int control_interp){/*{{{*/

	int         sidlist[NUMVERTICES];
	int         lidlist[NUMVERTICES];
	int         connectivity[NUMVERTICES];
	IssmPDouble values[NUMVERTICES];
	IssmPDouble gradients[NUMVERTICES];
	IssmDouble  value,gradient;

	/*Get relevant inputs*/
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

		GaussTria gauss;
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
void       Tria::CreateDistanceInputFromSegmentlist(IssmDouble* distances,int distanceenum){/*{{{*/

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
int        Tria::EdgeOnBaseIndex(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			return i;
		}
	}

	_printf_("list of vertices on bed: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on bed");
}
/*}}}*/
void       Tria::EdgeOnBaseIndices(int* pindex1,int* pindex2){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			return;
		}
	}

	_printf_("list of vertices on bed: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on bed");
}
/*}}}*/
int        Tria::EdgeOnSurfaceIndex(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			return i;
		}
	}

	_printf_("list of vertices on surface: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on surface");
}
/*}}}*/
void       Tria::EdgeOnSurfaceIndices(int* pindex1,int* pindex2){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[3][2] = {{1,2},{2,0},{0,1}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<3;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			return;
		}
	}

	_printf_("list of vertices on surface: "<<values[0]<<" "<<values[1]<<" "<<values[2]);
	_error_("Could not find 2 vertices on surface");
}
/*}}}*/
void       Tria::ElementCoordinates(Vector<IssmDouble>* vxe,Vector<IssmDouble>* vye,Vector<IssmDouble>* vze, Vector<IssmDouble>* vareae, bool spherical){/*{{{*/

	/*Look for x,y,z coordinates:*/
	IssmDouble xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Find centroid:*/
	IssmDouble xe=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	IssmDouble ye=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;
	IssmDouble ze=(xyz_list[0][2]+xyz_list[1][2]+xyz_list[2][2])/3.0;

	if(!spherical){
		IssmDouble area;
		vxe->SetValue(this->sid,xe,INS_VAL);
		vye->SetValue(this->sid,ye,INS_VAL);
		vze->SetValue(this->sid,ze,INS_VAL);
		area=this->GetAreaSpherical();
		vareae->SetValue(this->sid,area,INS_VAL);

		/*in addition, put in in the inputs:*/
		this->inputs->SetDoubleInput(AreaEnum,this->lid,area);
	}
	else _error_("spherical coordinates not supported yet!");
	return;
}
/*}}}*/
void       Tria::ElementCoordinates(Vector<IssmDouble>* vlonge,Vector<IssmDouble>* vlate,Vector<IssmDouble>* vareae){ /*{{{*/

	IssmDouble planetradius;

	/*Look for x,y,z coordinates:*/
	IssmDouble xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);
	IssmDouble area,xe,ye,ze,late,longe;

	/*Find centroid:*/
	xe=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	ye=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;
	ze=(xyz_list[0][2]+xyz_list[1][2]+xyz_list[2][2])/3.0;

	late= asin(ze/sqrt(pow(xe,2.0)+pow(ye,2.0)+pow(ze,2.0)));
	longe= atan2(ye,xe);
	area=this->GetAreaSpherical();

	vlonge->SetValue(this->sid,longe,INS_VAL);
	vlate->SetValue(this->sid,late,INS_VAL);
	vareae->SetValue(this->sid,area,INS_VAL);

	return;
}
/*}}}*/
void       Tria::ElementResponse(IssmDouble* presponse,int response_enum){/*{{{*/

	switch(response_enum){
		case MaterialsRheologyBbarEnum:
			*presponse=this->material->GetBbar(NULL);
			break;

		case VelEnum:{

			/*Get input:*/
			IssmDouble vel;
			Input* vel_input=this->GetInput(VelEnum); _assert_(vel_input);
			vel_input->GetInputAverage(&vel);

			/*Assign output pointers:*/
			*presponse=vel;}
			break;
		default:
			_error_("Response type " << EnumToStringx(response_enum) << " not supported yet!");
	}

}
/*}}}*/
void       Tria::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble xmin,ymin;
	IssmDouble xmax,ymax;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	xmin=xyz_list[0][0]; xmax=xyz_list[0][0];
	ymin=xyz_list[0][1]; ymax=xyz_list[0][1];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<xmin) xmin=xyz_list[i][0];
		if(xyz_list[i][0]>xmax) xmax=xyz_list[i][0];
		if(xyz_list[i][1]<ymin) ymin=xyz_list[i][1];
		if(xyz_list[i][1]>ymax) ymax=xyz_list[i][1];
	}

	*hx=xmax-xmin;
	*hy=ymax-ymin;
	*hz=0.;
}
/*}}}*/
int        Tria::FiniteElement(void){/*{{{*/
	return this->element_type;
}
/*}}}*/
IssmDouble Tria::FloatingArea(bool scaled){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,scalefactor,floatingarea;

	if(!IsIceInElement())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum && domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	floatingarea=(1-phi)*this->GetArea();
	if(scaled==true){
		Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
		scalefactor_input->GetInputAverage(&scalefactor);
		floatingarea=floatingarea*scalefactor;
	}

	/*Clean up and return*/
	return floatingarea;
}
/*}}}*/
void       Tria::FSContactMigration(Vector<IssmDouble>* vertex_sigmann,Vector<IssmDouble>* vertex_waterpressure){/*{{{*/

	if(!IsOnBase()) return;

	int approximation;
	this->Element::GetInputValue(&approximation,ApproximationEnum);

	if(approximation==HOApproximationEnum || approximation==SSAApproximationEnum || approximation==SSAHOApproximationEnum){
		_error_(" contact contiditon only works for FS elements");
	}
	/*Intermediaries*/
	IssmDouble  bed_normal[2],base[NUMVERTICES],bed[NUMVERTICES],surface[NUMVERTICES],phi[NUMVERTICES];
	IssmDouble  water_pressure[NUMVERTICES],pressureice[NUMVERTICES],pressure[NUMVERTICES];
	IssmDouble  sigmaxx[NUMVERTICES],sigmayy[NUMVERTICES],sigmaxy[NUMVERTICES],sigma_nn[NUMVERTICES];
	IssmDouble  viscosity,epsilon[NUMVERTICES];
	Element::GetInputListOnVertices(&base[0],BaseEnum);
	Element::GetInputListOnVertices(&bed[0],BedEnum);
	Element::GetInputListOnVertices(&surface[0],SurfaceEnum);
	Element::GetInputListOnVertices(&pressure[0],PressureEnum);
	Element::GetInputListOnVertices(&phi[0],MaskOceanLevelsetEnum);
	IssmDouble rho_ice   = FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = FindParam(MaterialsRhoSeawaterEnum);
	IssmDouble gravity   = FindParam(ConstantsGEnum);

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will be needing: */
	Input* vx_input       = this->GetInput(VxEnum);       _assert_(vx_input);
	Input* vy_input       = this->GetInput(VyEnum);       _assert_(vy_input);

	/*1. Recover stresses at the base*/
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);
		this->material->ViscosityFS(&viscosity,2,&xyz_list[0][0],&gauss,vx_input,vy_input,NULL);
		/*FIXME: this is for Hongju only*/
		//	pressureice[iv]=gravity*rho_ice*(surface[iv]-base[iv]);
		//	if (pressure[iv]/pressureice[iv]>1) pressure[iv]=pressureice[iv];

		/*Compute Stress*/
		sigmaxx[iv]=2*viscosity*epsilon[0]-pressure[iv];
		sigmayy[iv]=2*viscosity*epsilon[1]-pressure[iv];
		sigmaxy[iv]=2*viscosity*epsilon[2];
	}

	/*2. compute contact condition*/
	for(int i=0;i<NUMVERTICES;i++){
		/*If was grounded*/
		if (phi[i]>=0.){
			NormalBase(&bed_normal[0],&xyz_list[0][0]);
			sigma_nn[i]=-1*(sigmaxx[i]*bed_normal[0]*bed_normal[0] + sigmayy[i]*bed_normal[1]*bed_normal[1]+2*sigmaxy[i]*bed_normal[0]*bed_normal[1]);
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
			vertex_waterpressure->SetValue(vertices[i]->Pid(),+1.,ADD_VAL);
		}
	}
}
/*}}}*/
IssmDouble Tria::GetArea(void){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble x1,y1,x2,y2,x3,y3;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x1=xyz_list[0][0]; y1=xyz_list[0][1];
	x2=xyz_list[1][0]; y2=xyz_list[1][1];
	x3=xyz_list[2][0]; y3=xyz_list[2][1];

	_assert_(x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1>0);
	return (x2*y3 - y2*x3 + x1*y2 - y1*x2 + x3*y1 - y3*x1)/2;
}
/*}}}*/
IssmDouble Tria::GetHorizontalSurfaceArea(void){/*{{{*/

	return this->GetArea();
}
/*}}}*/
IssmDouble Tria::GetArea3D(void){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble x1,y1,z1,x2,y2,z2,x3,y3,z3;
	IssmDouble detm1,detm2,detm3;

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x1=xyz_list[0][0]; y1=xyz_list[0][1]; z1=xyz_list[0][2];
	x2=xyz_list[1][0]; y2=xyz_list[1][1]; z2=xyz_list[1][2];
	x3=xyz_list[2][0]; y3=xyz_list[2][1]; z3=xyz_list[2][2];

	detm1=x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2;
	detm2=y1*z2 - y2*z1 - y1*z3 + y3*z1 + y2*z3 - y3*z2;
	detm3=x2*z1 - x1*z2 + x1*z3 - x3*z1 - x2*z3 + x3*z2;

	return sqrt(pow(detm1,2) + pow(detm2,2) + pow(detm3,2))/2;
}
/*}}}*/
IssmDouble Tria::GetAreaIce(void){/*{{{*/

	/*return area of element covered by ice*/
	/*Intermediaries*/
	int numiceverts;
	IssmDouble area_fraction;
	IssmDouble s[2]; // s:fraction of intersected triangle edges that lie inside ice
	int* indices=NULL;

	this->GetLevelsetIntersection(&indices, &numiceverts, s, MaskIceLevelsetEnum, 0.);

	switch (numiceverts){
		case 0: // no vertex has ice: element is ice free
			area_fraction=0.;
			break;
		case 1: // one vertex has ice: get area of triangle
			area_fraction=s[0]*s[1];
			break;
		case 2: // two vertices have ice: get area of quadrangle
			area_fraction=s[0]+s[1]-s[0]*s[1];
			break;
		case NUMVERTICES: // all vertices have ice: return triangle area
			area_fraction=1.;
			break;
		default:
			_error_("Wrong number of ice vertices in Tria::GetAreaIce!");
			break;
	}
	_assert_((area_fraction>=0.) && (area_fraction<=1.));

	xDelete<int>(indices);
	return area_fraction*this->GetArea();
}/*}}}*/
IssmDouble Tria::GetAreaSpherical(void){/*{{{*/

	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble x1,y1,z1,x2,y2,z2,x3,y3,z3;
	IssmDouble arc12,arc23,arc31,semi_peri,excess;

	/*retrieve coordinates: lat,long,radius */
	::GetVerticesCoordinates(&llr_list[0][0],vertices,NUMVERTICES,spherical);
	x1=llr_list[0][0]/180.*M_PI; y1=llr_list[0][1]/180.*M_PI; z1=llr_list[0][2];
	x2=llr_list[1][0]/180.*M_PI; y2=llr_list[1][1]/180.*M_PI; z2=llr_list[1][2];
	x3=llr_list[2][0]/180.*M_PI; y3=llr_list[2][1]/180.*M_PI; z3=llr_list[2][2];

	/*compute great circle distance between vertices */
	arc12=2.*asin(sqrt(pow(sin(0.5*(x2-x1)),2)+cos(x1)*cos(x2)*pow(sin(0.5*(y2-y1)),2)));
	arc23=2.*asin(sqrt(pow(sin(0.5*(x3-x2)),2)+cos(x2)*cos(x3)*pow(sin(0.5*(y3-y2)),2)));
	arc31=2.*asin(sqrt(pow(sin(0.5*(x1-x3)),2)+cos(x3)*cos(x1)*pow(sin(0.5*(y1-y3)),2)));

	/*semi parameter */
	semi_peri=(arc12+arc23+arc31)/2;

	/*spherical excess */
	excess=4.*atan(sqrt(tan(semi_peri/2)*tan((semi_peri-arc12)/2)*tan((semi_peri-arc23)/2)*tan((semi_peri-arc31)/2)));

	/*area = excess*radius^2 */
	return excess*pow((z1+z2+z3)/3,2);
}
/*}}}*/
void       Tria::GetAreaCoordinates(IssmDouble* area_coordinates,IssmDouble* xyz_zero,IssmDouble* xyz_list,int numpoints){/*{{{*/
	/*Computeportion of the element that is grounded*/

	int         i,j,k;
	IssmDouble  area_init,area_portion;
	IssmDouble  xyz_bis[NUMVERTICES][3];

	area_init=GetArea();

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
int        Tria::GetElementType(){/*{{{*/

	/*return TriaRef field*/
	return this->element_type;

}
/*}}}*/
void       Tria::GetGroundedPart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* pmainlyfloating, int distance_enum, IssmDouble intrusion_distance){/*{{{*/
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
	*pmainlyfloating=floating;
}
/*}}}*/
IssmDouble Tria::GetGroundedPortion(IssmDouble* xyz_list){/*{{{*/
	/*Computeportion of the element that is grounded*/

	bool              mainlyfloating = true;
	int               domaintype,index1,index2;
	const IssmPDouble epsilon        = 1.e-15;
	IssmDouble        phi,s1,s2;
	IssmDouble        gl[NUMVERTICES];

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		this->EdgeOnBaseIndices(&index1,&index2);
		if(gl[index1]>0 && gl[index2]>0) phi=1; // All grounded
		else if(gl[index1]<0 && gl[index2]<0) phi=0; // All floating
		else if(gl[index1]<0 && gl[index2]>0){ //index2 grounded
			phi=1./(1.-gl[index1]/gl[index2]);
		}
		else if(gl[index2]<0 && gl[index1]>0){ //index1 grounded
			phi=1./(1.-gl[index2]/gl[index1]);
		}

	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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
	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	_assert_(phi<=1. && phi>=0.);
	return phi;
}
/*}}}*/
void       Tria::GetFractionGeometry(IssmDouble* weights, IssmDouble* pphi, int* ppoint1,IssmDouble* pfraction1,IssmDouble* pfraction2, bool* ptrapezeisnegative, IssmDouble* gl){/*{{{*/

	/*Computeportion of the element that is grounded*/
	bool               trapezeisnegative=true; //default value
	int                point;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         f1,f2,phi;

	/*Weights: */
	Gauss* gauss=NULL;
	IssmDouble loadweights_g[NUMVERTICES];
	IssmDouble total_weight=0;

	_assert_(!xIsNan<IssmDouble>(gl[0]));
	_assert_(!xIsNan<IssmDouble>(gl[1]));
	_assert_(!xIsNan<IssmDouble>(gl[2]));

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are positive or negative: */
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All positive
		point=0;
		f1=1.;
		f2=1.;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All negative
		point=0;
		f1=0.;
		f2=0.;
	}
	else{
		if(gl[0]*gl[1]*gl[2]<0) trapezeisnegative=false; //no matter what configuration, there has to be two positive vertices, which means the trapeze is positive.

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
	if(trapezeisnegative) phi=1-f1*f2;
	else phi=f1*f2;

	/*Compute weights:*/
	gauss = this->NewGauss(point,f1,f2,1-trapezeisnegative,2); //VV correction (16Nov2021)

	total_weight=0;
	for(int i=0;i<NUMVERTICES;i++)weights[i]=0;
	while(gauss->next()){
		TriaRef::GetNodalFunctions(&loadweights_g[0], gauss,P1Enum);
		for(int i=0;i<NUMVERTICES;i++)weights[i]+=loadweights_g[i]*gauss->weight;
		total_weight+=gauss->weight;
	}
	/*Normalizing to phi such that weights provide coefficients for integration over subelement (for averaging:phi*weights)*/  
	if(total_weight>0.) for(int i=0;i<NUMVERTICES;i++)weights[i]/=total_weight/phi; 
	else for(int i=0;i<NUMVERTICES;i++)weights[i]=0;

	/*Free resources:*/
	delete gauss;

	/*Assign output pointers:*/
	*pphi=phi;
	*ppoint1=point;
	*pfraction1=f1;
	*pfraction2=f2;
	*ptrapezeisnegative=trapezeisnegative;
}
/*}}}*/
IssmDouble Tria::GetTriangleAreaSpherical(IssmDouble xyz_list[3][3]){/*{{{*/

	IssmDouble x1,y1,z1,x2,y2,z2,x3,y3,z3;
	IssmDouble arc12,arc23,arc31,semi_peri,excess;
	IssmDouble lat1,lat2,lat3;
	IssmDouble long1,long2,long3;
	IssmDouble r1,r2,r3;

	/*retrieve x,y,z coordinates: */
	x1=xyz_list[0][0]; y1=xyz_list[0][1]; z1=xyz_list[0][2];
	x2=xyz_list[1][0]; y2=xyz_list[1][1]; z2=xyz_list[1][2];
	x3=xyz_list[2][0]; y3=xyz_list[2][1]; z3=xyz_list[2][2];

	/*Build lat,long, r:*/
	r1=sqrt(pow(x1,2.0)+pow(y1,2.0)+pow(z1,2.0));
	r2=sqrt(pow(x2,2.0)+pow(y2,2.0)+pow(z2,2.0));
	r3=sqrt(pow(x3,2.0)+pow(y3,2.0)+pow(z3,2.0));

	lat1=asin(z1/r1); long1=atan2(y1,x1);
	lat2=asin(z2/r2); long2=atan2(y2,x2);
	lat3=asin(z3/r3); long3=atan2(y3,x3);

	/*compute great circle distance between vertices */
	arc12=2.*asin(sqrt(pow(sin(0.5*(lat2-lat1)),2)+cos(lat1)*cos(lat2)*pow(sin(0.5*(long2-long1)),2)));
	arc23=2.*asin(sqrt(pow(sin(0.5*(lat3-lat2)),2)+cos(lat2)*cos(lat3)*pow(sin(0.5*(long3-long2)),2)));
	arc31=2.*asin(sqrt(pow(sin(0.5*(lat1-lat3)),2)+cos(lat3)*cos(lat1)*pow(sin(0.5*(long1-long3)),2)));

	/*semi parameter */
	semi_peri=(arc12+arc23+arc31)/2;

	/*spherical excess */
	excess=4.*atan(sqrt(tan(semi_peri/2)*tan((semi_peri-arc12)/2)*tan((semi_peri-arc23)/2)*tan((semi_peri-arc31)/2)));

	/*area = excess*radius^2 */
	return excess*pow((r1+r2+r3)/3,2);
}
/*}}}*/
void       Tria:: GetBarycenterFromLevelset(IssmDouble* platbar, IssmDouble* plongbar,IssmDouble phi,IssmDouble fraction1,IssmDouble fraction2,IssmDouble late, IssmDouble longe, int point1,int istrapeze1, IssmDouble planetradius){ /*{{{*/

	int i0,i1,i2;

	IssmDouble xyz0[3][3];
	IssmDouble barycenter[3]={0};
	IssmDouble centroid[3]={0};

	::GetVerticesCoordinates(&xyz0[0][0],vertices,NUMVERTICES); // initial triangle

	i0=point1;
	i1=(point1+1)%3;
	i2=(point1+2)%3;

	//Barycenter of the subelement triangle:
	for (int i=0;i<3;i++) barycenter[i]=xyz0[i0][i]*(3.0-fraction1-fraction2)/3.0 + fraction1/3.0*xyz0[i1][i] + fraction2/3.0*xyz0[i2][i];

	if (istrapeze1){
		centroid[0]=planetradius*cos(late*M_PI/180.0) * cos(longe*M_PI/180.0); //x
		centroid[1]=planetradius*cos(late*M_PI/180.0) * sin(longe*M_PI/180.0);  //y
		centroid[2]=planetradius*sin(late*M_PI/180.0);					//z

		// centroid_el *area_el= barycenter_triangle * area_triangle + barycenter_trapeze * area_trapeze
		// and phi_trapeze = area_trapeze/area_el = (1 - area_triangle/area_el)
		// => barycenter_trapeze = (centroid_el - barycenter_triangle * (1-phi_trapeze) )/phi_trapeze
		for (int i=0;i<3;i++) barycenter[i] =(centroid[i] -barycenter[i]*(1.0-phi))/phi;

	}

	//recompute planetradius from the barycenter onwards:
	planetradius=sqrt( pow(barycenter[0],2.0)+ pow(barycenter[1],2.0)+ pow(barycenter[2],2.0));

	*platbar=asin(barycenter[2]/planetradius)*180.0/M_PI;
	*plongbar=atan2(barycenter[1],barycenter[0])*180.0/M_PI;

} /*}}}*/
void       Tria::GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area,  int levelsetenum){ /*{{{*/

	IssmDouble phi;
	IssmDouble fraction1,fraction2;
	bool istrapeze1;  
	bool flip1=false;
	int  point1;
	IssmDouble levelset[NUMVERTICES];
	IssmDouble planetradius;

	this->parameters->FindParam(&planetradius,SolidearthPlanetRadiusEnum);

	//figure out if we are flipping the levelsets: 
	if(levelsetenum<0){
		levelsetenum=-levelsetenum;
		flip1=true;
	}
	//figure out area where we have loads
	Element::GetInputListOnVertices(&levelset[0],levelsetenum);
	if(flip1)for(int i=0;i<NUMVERTICES;i++)levelset[i]=-levelset[i];

	//compute sea level load weights
	this->GetFractionGeometry(loadweights,&phi,&point1,&fraction1,&fraction2,&istrapeze1,levelset);

	//failsafe for phi so small that GetFractionGeometry returns 0	
	if (phi==0) phi=1e-16;

	for (int i=0;i<NUMVERTICES;i++) loadweights[i]/=phi;
	this->GetBarycenterFromLevelset(platbar,plongbar, phi, fraction1, fraction2, late, longe, point1,istrapeze1,planetradius);

	/*assign output pointers:*/
	*ploadarea=phi*area;

} /*}}}*/
void       Tria::GetNodalWeightsAndAreaAndCentroidsFromLeveset(IssmDouble* loadweights, IssmDouble* ploadarea, IssmDouble* platbar, IssmDouble* plongbar, IssmDouble late, IssmDouble longe, IssmDouble area,  int levelset1enum, int levelset2enum){ /*{{{*/

	bool istrapeze1, istrapeze2;
	IssmDouble phi1,phi2, d,e,f,g,h1,h2;
	int point1, point2,  i0,i1,i2,j0,j1,j2;
	IssmDouble weights1[3],weights2[3];
	IssmDouble levelset1[3];
	IssmDouble levelset2[3];

	bool flip1=false;
	bool flip2=false;

	IssmDouble xyz0[3][3];
	IssmDouble xyz1[3][3]={0};
	IssmDouble xyz2[3][3]={0};
	IssmDouble xyz3[3][3]={0};
	IssmDouble xyz[8][3]={0};
	IssmDouble w[8][NUMVERTICES]={0};
	IssmDouble areasub=0;
	IssmDouble area1=0;
	IssmDouble area2=0;
	IssmDouble area3=0;

	int tria0[3]={0,1,2};
	int tria1[3]={-1};
	int tria2[3]={-1};
	int tria3[3]={-1};

	IssmDouble w1[3][3]={0};
	IssmDouble w2[3][3]={0};
	IssmDouble w3[3][3]={0};

	IssmDouble barycenter[3]={0};
	IssmDouble planetradius;

	//figure out if we are flipping the levelsets: 
	if(levelset1enum<0){
		levelset1enum=-levelset1enum;
		flip1=true;
	}
	if(levelset2enum<0){
		levelset2enum=-levelset2enum;
		flip2=true;
	}

	//recover levelsets: 
	Element::GetInputListOnVertices(&levelset1[0],levelset1enum);
	if(flip1)for(int i=0;i<NUMVERTICES;i++)levelset1[i]=-levelset1[i];
	Element::GetInputListOnVertices(&levelset2[0],levelset2enum);
	if(flip2)for(int i=0;i<NUMVERTICES;i++)levelset2[i]=-levelset2[i];

	//We want the fraction of the element where both levelsets are negative.
	//Early return if either of them is >=0 on all vertices
	if (levelset1[0]>=0 && levelset1[1]>=0 && levelset1[2]>=0) {
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]=0.0;
		*ploadarea= 0.0;
		*platbar=late; //just default to centroid of triangle
		*plongbar=longe; 
		return;
	}
	if (levelset2[0]>=0 && levelset2[1]>=0 && levelset2[2]>=0) {
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]=0.0;
		*ploadarea= 0.0;
		*platbar=late; //just default to centroid of triangle
		*plongbar=longe; 
		return;
	}

	//If everyone is negative, no need to calculate any fraction
	if (levelset1[0]<=0 && levelset1[1]<=0 && levelset1[2]<=0 && levelset2[0]<=0 && levelset2[1]<=0 && levelset2[2]<=0) {
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]=1.0/NUMVERTICES;
		*ploadarea= area;
		*platbar=late;
		*plongbar=longe;
		return;
	}

	/*recover planet radius:*/
	this->parameters->FindParam(&planetradius,SolidearthPlanetRadiusEnum);

	//If just one levelset is all negative, just take the partitioning of the other, no interaction between them
	if (levelset1[0]<=0 && levelset1[1]<=0 && levelset1[2]<=0) {
		this->GetFractionGeometry(loadweights,&phi2,&point2,&f,&g,&istrapeze2,levelset2);
		this->GetBarycenterFromLevelset(platbar,plongbar, phi2, f, g, late, longe, point2,istrapeze2,planetradius);
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]/=phi2;
		*ploadarea=area*phi2;
		return;
	}
	if (levelset2[0]<=0 && levelset2[1]<=0 && levelset2[2]<=0) {
		this->GetFractionGeometry(loadweights,&phi1,&point1,&d,&e,&istrapeze1,levelset1);
		this->GetBarycenterFromLevelset(platbar,plongbar, phi1, d, e, late, longe, point1,istrapeze1,planetradius);
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]/=phi1;
		*ploadarea=area*phi1;
		return;
	}

	this->GetFractionGeometry(&weights1[0],&phi1,&point1,&d,&e,&istrapeze1,levelset1);
	this->GetFractionGeometry(&weights2[0],&phi2,&point2,&f,&g,&istrapeze2,levelset2);

	//Early return if levelsets are not independent
	if (istrapeze1==istrapeze2 && point1==point2 && phi1==phi2){
		//the two levelsets are redundant: levelset1 = positivescalar * levelset2
		this->GetBarycenterFromLevelset(platbar,plongbar, phi1, d, e, late, longe, point1,istrapeze1,planetradius);
		*ploadarea=area*phi1;
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]=weights1[i]/phi1;
		return;
	}
	if (istrapeze1!=istrapeze2 && point1==point2 && phi1==(1.0-phi2)){
		//the two levelsets are incompatible: levelset1 = negativescalar * levelset2
		*plongbar=longe;
		*platbar=late;
		*ploadarea=0.0;
		for (int i=0;i<NUMVERTICES;i++) loadweights[i]=0.0;
		return;
	}

	::GetVerticesCoordinates(&xyz0[0][0],vertices,NUMVERTICES); // initial triangle

	//Let our element be triangle ABC with:
	i0=point1; //A
	i1=(point1+1)%3; //B
	i2=(point1+2)%3; //C

	j0=point2; //Can be A, B or C
	j1=(point2+1)%3; //anticlockwise point from j0
	j2=(point2+2)%3; //clockwise point from j0

	/* Below we define the relative fractional lengths of ABC where the zero-level contours of the two level sets intersect with ABC and each other. For example D is the intersection of level set 1 with side AB with fractional length d=[AD]/[AB].

	   levelset1 intersects ABC on D and E:
	   A------D---B
	   <--d--->		

	   A-----E----C
	   <--e-->

	   levelset2 intersects ABC on F and G:
	   j0---F------j1
	   <--f->

	   j0-------G--j2
	   <---g---->

	   levelset1 and 2 intersect on H (when that intersection exists inside the element)
	   D----H------E
	   <-h1->

	   F-----H----G
	   <--h2->
	*/

	if (point2==i0){
		h1= g*(d-f)/(d*g-e*f);
		h2= e/g * h1;
	}
        else if (point2==i1){
		h1=f*(1.0-g-d)/(f*(e-d)-e*g);
		h2= 1.0-e/f * h1;
	}
	else if (point2==i2){
		h1= (g*(1.0-f-d)+f*d)/(g*(e-d) +f*d); 
		h2= (d*(f+e-1))/(g*(e-d) +f*d);
	}

	//interpolant weights of each point. Any field F[0,1,2] provided at the original vertices [0,1,2] will be equal on point k to sum_i (F[i] * w[k][i])
	w[i0][i0]=1; //A
	w[i1][i1]=1; //B
	w[i2][i2]=1; //C
	w[3][i0]=1.0-d; w[3][i1]=d; //D
	w[4][i0]=1.0-e; w[4][i2]=e; //E
	w[5][j0]=1.0-f; w[5][j1]=f; //F
	w[6][j0]=1.0-g; w[6][j2]=g; //G
	for (int j=0;j<3;j++) w[7][j]=w[3][j]*(1.0-h1) + w[4][j]*h1; //H: we interpolate the intersection point H between D and E at fraction h1

	for (int k=0;k<8;k++){
		for (int i=0;i<NUMVERTICES;i++) {
			for (int j=0;j<3;j++) xyz[k][j]+=xyz0[i][j]*w[k][i];
		}
	}

		//point2 can be either i0,i1 or i2. We start the search with i1 and i2 as they have less computational cost in ifs
		if(point2==i2){ /*{{{*/
			if (e>1.0-f){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=5; tria1[1]= 7; tria1[2]= 4; // FHE
					area1=h2*(e+f-1.0)*g;
				}
				else if (!istrapeze1 && istrapeze2){
					tria1[0]=i0; tria1[1]= 3; tria1[2]= 5; //ADF
					area1=d*(1.0-f);
					tria2[0]=3; tria2[1]= 7; tria2[2]= 5; //DHF
					area2=d*h1*(e+f-1.0);
				}
				else if (istrapeze1 && !istrapeze2){
					tria1[0]=7; tria1[1]= 6; tria1[2]= 4; //HGE
					area1=g*(1.0-h2)*(e+f-1.0);
					tria2[0]=4; tria2[1]= 6; tria2[2]= i2; //EGC
					area2=g*(1.0-e);
				}
				else { //istrapeze1 && istrapeze2
					tria1[0]=3; tria1[1]= i1; tria1[2]= 6; //DBG
					area1=(1.0-d)*(1.0-g);
					tria2[0]=3; tria2[1]= 6; tria2[2]= 7; //DGH
					area2=g*((1.0-f)*(1.0-h2)+h2*e)+d*(1.0-e-g);
				}  /*}}}*/
			}
			else if (e<=1.0-f){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
				}
				else if (!istrapeze1 && istrapeze2){
					tria1[0]=i0; tria1[1]= 3; tria1[2]= 4; //ADE
					area1=d*e;
				}
				else if (istrapeze1 && !istrapeze2){
					tria1[0]=5; tria1[1]= 6; tria1[2]= i2; //FGC
					area1=f*g;
				}
				else { //istrapeze1 && istrapeze2
					tria1[0]=3; tria1[1]= i1; tria1[2]= 5; //DBF
			                area1=(1.0-d)*(1.0-f);
					tria2[0]=4; tria2[1]= 3; tria2[2]= 5;  //EDF
			                area2=d*(1.0-e-f);
					tria3[0]=5; tria3[1]= i1; tria3[2]= 6; //FBG
			                area3=f*(1.0-g);
				} /*}}}*/
			} 
		}/*}}}*/    
		else if(point2==i1){ /*{{{*/
			if (d>1.0-g){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=6; tria1[1]= 3; tria1[2]= 7; //GDH
					area1=(1.0-h2)*(d+g-1.0)*f;
				}
				else if (!istrapeze1 && istrapeze2){
					tria1[0]=i0; tria1[1]= 6; tria1[2]= 4; //AGE
			                area1=(1.0-g)*e;
					tria2[0]=6; tria2[1]= 7; tria2[2]= 4; //GHE
			                area2=e*(1.0-h1)*(g+d-1.0);
				}
				else if (istrapeze1 && !istrapeze2){
					tria1[0]=i1; tria1[1]= 5; tria1[2]= 3; //BFD
			                area1=(1.0-d)*f;
					tria2[0]=3; tria2[1]= 5; tria2[2]= 7; //DFH
					area2=f*h2*(d+g-1.0);
				}
				else { //istrapeze1 && istrapeze2
					tria1[0]=7; tria1[1]= 5; tria1[2]= 4; //HFE
			                area1=e*((1.0-d)*(1.0-h1)+h1*g)+f*(1.0-g-e);
					tria2[0]=4; tria2[1]= 5; tria2[2]= i2;//EFC
			                area2=(1.0-e)*(1.0-f);
				}  /*}}}*/
			}
			else if (d<=1.0-g){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
				}
				else if (!istrapeze1 && istrapeze2){
					tria1[0]=i0; tria1[1]= 3; tria1[2]= 4; //ADE
			                area1=d*e;
				}
				else if (istrapeze1 && !istrapeze2){
					tria1[0]=6; tria1[1]= i1; tria1[2]= 5; //GBF
			                area1=f*g;
				}
				else { //istrapeze1 && istrapeze2
					tria1[0]=3; tria1[1]= 6; tria1[2]= 4; //DGE
					area1=e*(1.0-d-g);
					tria2[0]=6; tria2[1]= i2; tria2[2]= 4; //GCE
					area2=(1.0-g)*(1.0-e);
					tria3[0]=6; tria3[1]= 5; tria3[2]= i2; //GFC
					area3=g*(1.0-f);
				} /*}}}*/
			}

		}/*}}}*/
		else{ /*{{{*/
			if (d<=f && e>=g){  /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=i0; tria1[1]= 3; tria1[2]= 7; //ADH
			                area1=h1*d*e;
					tria2[0]=i0; tria2[1]= 7; tria2[2]= 6; //AHG
			                area2=(1.0-h2)*f*g;
				}
				else if (!istrapeze1 && istrapeze2){
					tria1[0]=6; tria1[1]= 7; tria1[2]= 4; //GHE
			                area1=(e-g)*d*(1.0-h1);
				}
				else if (istrapeze1 && !istrapeze2){
					tria1[0]=3; tria1[1]= 5; tria1[2]= 7; //DFH
					area1=(f-d)*g*h2;
				}
				else { //istrapeze1 && istrapeze2
					tria1[0]=5; tria1[1]= i1; tria1[2]= 7; //FBH
			                area1=g*h2*(1.0-f);
					tria2[0]=7; tria2[1]= i1; tria2[2]= i2;//HBC
			                area2=1.0+d*(h1-e*h1-1.0)+g*h2*(d-1.0);
					tria3[0]=7; tria3[1]= i2; tria3[2]= 4; //HCE
			                area3=d*(1.0-h1)*(1.0-e);
				} /*}}}*/
			}
			else if (d>=f && e<=g){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=i0; tria1[1]= 5; tria1[2]= 7; //AFH
			                area1=h2*f*g;
					tria2[0]=i0; tria2[1]= 7; tria2[2]= 4; //AHE
			                area2=(1.0-h1)*d*e;
				}else if (!istrapeze1 && istrapeze2){
					tria1[0]=5; tria1[1]= 3; tria1[2]= 7; //FDH
			                area1=(d-f)*e*h1;
				}else if (istrapeze1 && !istrapeze2){
					tria1[0]=4; tria1[1]= 7; tria1[2]= 6; //EHG
			                area1=(g-e)*f*(1.0-h2);
				}else { //istrapeze1 && istrapeze2
					tria1[0]=3; tria1[1]= i1; tria1[2]= 7; //DBH
			                area1=e*h1*(1.0-d);
					tria2[0]=7; tria2[1]= i1; tria2[2]= 6; //HCG
			                area2=f*(1.0-h2)*(1.0-g);
					tria3[0]=6; tria3[1]= i1; tria3[2]= i2; //HBC
			                area3=1.0+f*(h2-g*h2-1.0)+e*h1*(f-1.0);
				}  /*}}}*/
			}
			else if (d<=f && e<=g){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=i0; tria1[1]= 3; tria1[2]= 4;//ADE
					area1=d*e;
				}else if (!istrapeze1 && istrapeze2){
				}else if (istrapeze1 && !istrapeze2){
					tria1[0]=3; tria1[1]= 5; tria1[2]= 6; //DFG
			                area1=g*(f-d);
					tria2[0]=3; tria2[1]= 6; tria2[2]= 4; //DGE
			                area2=d*(g-e);
				}else { //istrapeze1 && istrapeze2
					tria1[0]=5; tria1[1]= i2; tria1[2]= 6; //FCG
			                area1=f*(1.0-g);
					tria2[0]=5; tria2[1]= i1; tria2[2]= i2;//FBC
			                area2=1.0-f;
				}  /*}}}*/
			}
			else if (d>=f && e>=g){ /*{{{*/
				if (!istrapeze1 && !istrapeze2){
					tria1[0]=i0; tria1[1]= 5; tria1[2]= 6; //AFG
			                area1=f*g;
				}else if (!istrapeze1 && istrapeze2){
					tria1[0]=5; tria1[1]= 4; tria1[2]= 6; //FEG
			                area1=f*(e-g);
					tria2[0]=5; tria2[1]= 3; tria2[2]= 4; //FDE
			                area2=e*(d-f);
				}else if (istrapeze1 && !istrapeze2){
				}else { //istrapeze1 && istrapeze2
					tria1[0]=3; tria1[1]= i1; tria1[2]= i2; //DBC
			                area1=1.0-d;
					tria2[0]=3; tria2[1]= i2; tria2[2]= 4;//DCE
			                area2=d*(1.0-e);
				}  /*}}}*/
			} 
		} /*}}}*/

	if(tria1[0]>-1){ 
		for (int i=0;i<NUMVERTICES;i++){
			for (int j=0;j<3;j++) {
				xyz1[i][j]=xyz[tria1[i]][j];
				w1[i][j]=w[tria1[i]][j];
			}
		}
		area1*=area; //dimensionalize the fractional area from [0:1] to [0:area]
	}
	if(tria2[0]>-1){ 
		for (int i=0;i<NUMVERTICES;i++){
			for (int j=0;j<3;j++) {
				xyz2[i][j]=xyz[tria2[i]][j];
				w2[i][j]=w[tria2[i]][j];
			}
		}
		area2*=area;
	}
	if(tria3[0]>-1){ 
		for (int i=0;i<NUMVERTICES;i++){
			for (int j=0;j<3;j++) {
				xyz3[i][j]=xyz[tria3[i]][j];
				w3[i][j]=w[tria3[i]][j];
			}
		}
		area3*=area;
	}

	areasub=area1+area2+area3;

	if (areasub>0){
		for (int j=0;j<3;j++){
			for (int i=0;i<NUMVERTICES;i++) {
				loadweights[j]+=w1[i][j]*area1 + w2[i][j]*area2 + w3[i][j]*area3;
				barycenter[j]+=xyz1[i][j]*area1+xyz2[i][j]*area2+xyz3[i][j]*area3;
			}
			loadweights[j]/=areasub*3.0;
			barycenter[j]/=areasub *3.0;
		}
		*platbar=asin(barycenter[2]/sqrt(pow(barycenter[0],2.0)+pow(barycenter[1],2.0)+pow(barycenter[2],2.0)))*180.0/M_PI;
		*plongbar=atan2(barycenter[1],barycenter[0])*180.0/M_PI;
	} 
	else {
		for(int j=0;j<3;j++)loadweights[j]=0.0;
		*platbar=late;
		*plongbar=longe;
	}
	*ploadarea=areasub;

} /*}}}*/
IssmDouble Tria::GetIcefrontArea(){/*{{{*/

	IssmDouble  bed[NUMVERTICES];
	IssmDouble	Haverage,frontarea;
	IssmDouble  x1,y1,x2,y2,distance;
	IssmDouble lsf[NUMVERTICES], Haux[NUMVERTICES], surfaces[NUMVERTICES], bases[NUMVERTICES];
	int* indices=NULL;

	/*Return if no ice front present*/
	if(!IsZeroLevelset(MaskIceLevelsetEnum)) return 0;
	//if(!this->IsIcefront()) return 0.;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&bed[0],BedEnum);
	Element::GetInputListOnVertices(&surfaces[0],SurfaceEnum);
	Element::GetInputListOnVertices(&bases[0],BaseEnum);
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);

	/*Only continue if all 3 vertices are below sea level*/
	for(int i=0;i<NUMVERTICES;i++) if(bed[i]>=0.) return 0.;

	/*2. Find coordinates of where levelset crosses 0*/
	int         numiceverts;
	IssmDouble  s[2],x[2],y[2];
	this->GetLevelsetIntersection(&indices, &numiceverts, &s[0],MaskIceLevelsetEnum,0.);
	_assert_(numiceverts);
	if(numiceverts>2){
		Input* ls_input = this->GetInput(MaskIceLevelsetEnum);
		ls_input->Echo();
	}

	/*3 Write coordinates*/
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);
	int counter = 0;
	if((numiceverts>0) && (numiceverts<NUMVERTICES)){
		for(int i=0;i<numiceverts;i++){
			for(int n=numiceverts;n<NUMVERTICES;n++){ // iterate over no-ice vertices
				x[counter] = xyz_list[indices[i]][0]+s[counter]*(xyz_list[indices[n]][0]-xyz_list[indices[i]][0]);
				y[counter] = xyz_list[indices[i]][1]+s[counter]*(xyz_list[indices[n]][1]-xyz_list[indices[i]][1]);
				counter++;
			}
		}
	}
	else if(numiceverts==NUMVERTICES){ //NUMVERTICES ice vertices: calving front lies on element edge

		for(int i=0;i<NUMVERTICES;i++){
			if(lsf[indices[i]]==0.){
				x[counter]=xyz_list[indices[i]][0];
				y[counter]=xyz_list[indices[i]][1];
				counter++;
			}
			if(counter==2) break;
		}
		if(counter==1){
			/*We actually have only 1 vertex on levelset, write a single point as a segment*/
			x[counter]=x[0];
			y[counter]=y[0];
			counter++;
		}
	}
	else{
		_error_("not sure what's going on here...");
	}
	x1=x[0]; y1=y[0]; x2=x[1]; y2=y[1];
	distance=sqrt(pow((x1-x2),2)+pow((y1-y2),2));
	if(distance<1e-3) return 0.;

	IssmDouble H[4];
	for(int iv=0;iv<NUMVERTICES;iv++) Haux[iv]=-bed[indices[iv]]; //sort bed in ice/noice
	xDelete<int>(indices);

	switch(numiceverts){
		case 1: // average over triangle
			H[0]=Haux[0];
			H[1]=Haux[0]+s[0]*(Haux[1]-Haux[0]);
			H[2]=Haux[0]+s[1]*(Haux[2]-Haux[0]);
			Haverage=(H[1]+H[2])/2;
			break;
		case 2: // average over quadrangle
			H[0]=Haux[0];
			H[1]=Haux[1];
			H[2]=Haux[0]+s[0]*(Haux[2]-Haux[0]);
			H[3]=Haux[1]+s[1]*(Haux[2]-Haux[1]);
			Haverage=(H[2]+H[3])/2;
			break;
		case 3:
			if(counter==1) distance = 0; //front has 0 width on this element because levelset is 0 at a single vertex
			else if(counter==2){ //two vertices with levelset=0: averaging ice front depth over both
				Haverage = 0;
				for(int i=0;i<NUMVERTICES;i++){
					if(lsf[indices[i]]==0.) Haverage -= Haux[indices[i]]/2;
					if(Haverage<Haux[indices[i]]/2-1e-3) break; //done with the two vertices
				}
			}
			break;
		default:
			_error_("Number of ice covered vertices wrong in Tria::GetIceFrontArea(void)");
			break;
	}
	frontarea=distance*Haverage;

	_assert_(frontarea>0);
	return frontarea;
}
/*}}}*/
void       Tria::GetIcefrontCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum){/*{{{*/

	/* Intermediaries */
	IssmDouble  levelset[NUMVERTICES];
	int         indicesfront[NUMVERTICES];

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&levelset[0],levelsetenum);

	/* Get nodes where there is no ice */
	int num_frontnodes=0;
	for(int i=0;i<NUMVERTICES;i++){
		if(levelset[i]>=0.){
			indicesfront[num_frontnodes]=i;
			num_frontnodes++;
		}
	}
	_assert_(num_frontnodes==2);

	/* arrange order of frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES+indicesfront[0]-indicesfront[1])%NUMVERTICES!=NUMVERTICES-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}

	IssmDouble* xyz_front = xNew<IssmDouble>(3*2);
	/* Return nodes */
	for(int i=0;i<2;i++){
		for(int j=0;j<3;j++){
			xyz_front[3*i+j]=xyz_list[3*indicesfront[i]+j];
		}
	}

	*pxyz_front=xyz_front;
}/*}}}*/
Input*    Tria::GetInput(int inputenum){/*{{{*/

	/*Get Input from dataset*/
	if(this->iscollapsed){
		PentaInput* input = this->inputs->GetPentaInput(inputenum);
		if(!input) return input;

		this->InputServe(input);
		return input;
	}
	else{
		TriaInput* input = this->inputs->GetTriaInput(inputenum);
		if(!input) return input;
		this->InputServe(input);
		return input;
	}
}/*}}}*/
Input*    Tria::GetInput(int inputenum,IssmDouble time){/*{{{*/

	/*Get Input from dataset*/
	if(this->iscollapsed){
		PentaInput* input = this->inputs->GetPentaInput(inputenum,time);
		if(!input) return input;

		this->InputServe(input);
		return input;
	}
	else{
		TriaInput* input = this->inputs->GetTriaInput(inputenum,time);
		if(!input) return input;

		this->InputServe(input);
		return input;
	}
}/*}}}*/
Input*    Tria::GetInput(int inputenum,IssmDouble start_time, IssmDouble end_time, int averaging_method){/*{{{*/

	/*Get Input from dataset*/
	if(this->iscollapsed){
		PentaInput* input = this->inputs->GetPentaInput(inputenum,start_time,end_time,averaging_method);
		if(!input) return input;

		this->InputServe(input);
		return input;
	}
	else{
		TriaInput* input = this->inputs->GetTriaInput(inputenum,start_time,end_time,averaging_method);
		if(!input) return input;

		this->InputServe(input);
		return input;
	}
}/*}}}*/
void       Tria::GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/* Start looping on the number of vertices: */
	if(input){
		GaussTria gauss;
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
void       Tria::GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/*What type of finite element are we dealing with?*/
	int fe       = this->FiniteElement();
	int numnodes = this->GetNumberOfNodes();

	/* Start looping on the number of vertices: */
	if(input){
		GaussTria gauss;
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
void       Tria::InputServe(Input* input_in){/*{{{*/

	/*Return NULL pointer if input is NULL*/
	if(!input_in) return;

	/*Get Input from dataset*/
	if(this->iscollapsed){
		_assert_(input_in->ObjectEnum()==PentaInputEnum);
		PentaInput* input = xDynamicCast<PentaInput*>(input_in);

		/*Intermediaries*/
		int numindices;
		int indices[3];

		/*Check interpolation*/
		int interpolation = input->GetInterpolation();
		switch(interpolation){
			case P0Enum:
				numindices = 1;
				indices[0] = this->lid;
				input->Serve(numindices,&indices[0]);
				break;
			case P1Enum:
				numindices = 3;
				for(int i=0;i<3;i++) indices[i] = vertices[i]->lid;
				input->Serve(numindices,&indices[0]);
				break;
			case P1xP2Enum:
			case P1xP3Enum:
			case P1xP4Enum:
			case P1DGEnum:
			case P1bubbleEnum:
			case P2Enum:
				input->ServeCollapsed(this->lid,this->iscollapsed);
				break;
			default: _error_("interpolation "<<EnumToStringx(interpolation)<<" not supported");
		}

		/*Flag as collapsed for later use*/
		input->SetServeCollapsed(true);
		return;
	}
	else{
		_assert_(input_in->ObjectEnum()==TriaInputEnum);
		TriaInput* input = xDynamicCast<TriaInput*>(input_in);

		/*Intermediaries*/
		int numindices;
		int indices[7];

		/*Check interpolation*/
		int interpolation = input->GetInterpolation();
		switch(interpolation){
			case P0Enum:
				numindices = 1;
				indices[0] = this->lid;
				input->Serve(numindices,&indices[0]);
				break;
			case P1Enum:
				numindices = 3;
				for(int i=0;i<3;i++) indices[i] = vertices[i]->lid;
				input->Serve(numindices,&indices[0]);
				break;
			case P1DGEnum:
				numindices = 3;
				input->Serve(this->lid,numindices);
				break;
			default:
				input->Serve(this->lid,this->GetNumberOfNodes(interpolation));
				break;
		}
		return;
	}
}/*}}}*/
DatasetInput* Tria::GetDatasetInput(int inputenum){/*{{{*/

	DatasetInput* datasetinput = this->inputs->GetDatasetInput(inputenum);
	if(!datasetinput) return NULL;

	for(int i=0;i<datasetinput->GetNumIds();i++){

		/*Get Input from dataset*/
		if(this->iscollapsed){

			PentaInput* input = datasetinput->GetPentaInputByOffset(i); _assert_(input);

			/*Intermediaries*/
			int numindices;
			int indices[3];

			/*Check interpolation*/
			int interpolation = input->GetInterpolation();
			switch(interpolation){
				case P0Enum:
					numindices = 1;
					indices[0] = this->lid;
					input->Serve(numindices,&indices[0]);
					break;
				case P1Enum:
					numindices = 3;
					for(int i=0;i<3;i++) indices[i] = vertices[i]->lid;
					input->Serve(numindices,&indices[0]);
					break;
				case P1DGEnum:
				case P1bubbleEnum:
					input->ServeCollapsed(this->lid,this->iscollapsed);
					break;
				default: _error_("interpolation "<<EnumToStringx(interpolation)<<" not supported");
			}

			/*Flag as collapsed for later use*/
			input->SetServeCollapsed(true);
		}
		else{

			TriaInput* input = datasetinput->GetTriaInputByOffset(i); _assert_(input);

			/*Intermediaries*/
			int numindices;
			int indices[7];

			/*Check interpolation*/
			int interpolation = input->GetInterpolation();
			switch(interpolation){
				case P0Enum:
					numindices = 1;
					indices[0] = this->lid;
					input->Serve(numindices,&indices[0]);
					break;
				case P1Enum:
					numindices = 3;
					for(int i=0;i<3;i++) indices[i] = vertices[i]->lid;
					input->Serve(numindices,&indices[0]);
					break;
				case P1DGEnum:
					numindices = 3;
					input->Serve(this->lid,numindices);
					break;
				default: _error_("interpolation "<<EnumToStringx(interpolation)<<" not supported");
			}

		}
	}

	return datasetinput;
}/*}}}*/
void       Tria::CreateInputTimeAverage(int transientinput_enum,int averagedinput_enum,IssmDouble start_time,IssmDouble end_time, int averaging_method){/*{{{*/

	_assert_(end_time>start_time);

	/*Get transient input time steps*/
	TransientInput* transient_input  = this->inputs->GetTransientInput(transientinput_enum);
	TriaInput* averaged_input = transient_input->GetTriaInput(start_time,end_time,averaging_method);
	Input* averaged_copy = averaged_input->copy();

	averaged_copy->ChangeEnum(averagedinput_enum);
	this->inputs->AddInput(averaged_copy);
}
/*}}}*/
void       Tria::GetInputAveragesUpToCurrentTime(int input_enum,IssmDouble** pvalues, IssmDouble** ptimes, int* pnumtimes, IssmDouble currenttime){/*{{{*/

	/*Get transient input time steps*/
	int         numtimesteps;
	IssmDouble *timesteps    = NULL;
	TransientInput* transient_input  = this->inputs->GetTransientInput(input_enum);

	transient_input->GetAllTimes(&timesteps,&numtimesteps);

	/*Figure out how many time steps we are going to return: */
	int  numsteps               = 0;
	bool iscurrenttime_included = false;
	for(int i=0;i<numtimesteps;i++){
		if(timesteps[i]==currenttime) iscurrenttime_included=true;
		if(timesteps[i]>currenttime)  break;
		else numsteps++;
	}
	if(iscurrenttime_included==false)numsteps++;

	/*allocate: */
	IssmDouble* times=xNew<IssmDouble>(numsteps);
	IssmDouble* values=xNew<IssmDouble>(numsteps);

	for(int i=0;i<numsteps;i++){
		if((iscurrenttime_included==false) && (i==(numsteps-1))){
			Input* input = this->GetInput(input_enum,currenttime);
			input->GetInputAverage(&values[i]);
			times[i]=currenttime;
		}
		else{
			TriaInput* input = transient_input->GetTriaInput(i);
			this->InputServe(input);
			input->GetInputAverage(&values[i]);
			times[i]=timesteps[i];
		}
	}

	/*Assign output pointers*/
	xDelete<IssmDouble>(timesteps);
	*pvalues=values;
	*ptimes=times;
	*pnumtimes=numtimesteps;
}
/*}}}*/
void       Tria::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=this->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	int index = this->GetNodeIndex(node);

	GaussTria gauss;
	gauss.GaussNode(this->element_type,index);
	input->GetInputValue(pvalue,&gauss);
}
/*}}}*/
void       Tria::GetInputValue(IssmDouble* pvalue,Vertex* vertex,int enumtype){/*{{{*/

	Input* input=this->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	int index = this->GetVertexIndex(vertex);

	GaussTria gauss;
	gauss.GaussVertex(index);
	input->GetInputValue(pvalue,&gauss);
}
/*}}}*/
void       Tria::GetLevelCoordinates(IssmDouble** pxyz_front,IssmDouble* xyz_list,int levelsetenum,IssmDouble level){/*{{{*/

	/* Intermediaries */
	int i, dir,nrfrontnodes;
	IssmDouble  levelset[NUMVERTICES];
	int indicesfront[NUMVERTICES];

	/*Recover parameters and values*/
	Element::GetInputListOnVertices(&levelset[0],levelsetenum);

	/* Get nodes where there is no ice */
	nrfrontnodes=0;
	for(i=0;i<NUMVERTICES;i++){
		if(levelset[i]==level){
			indicesfront[nrfrontnodes]=i;
			nrfrontnodes++;
		}
	}

	_assert_(nrfrontnodes==2);

	/* arrange order of frontnodes such that they are oriented counterclockwise */
	if((NUMVERTICES+indicesfront[0]-indicesfront[1])%NUMVERTICES!=NUMVERTICES-1){
		int index=indicesfront[0];
		indicesfront[0]=indicesfront[1];
		indicesfront[1]=index;
	}

	IssmDouble* xyz_front = xNew<IssmDouble>(3*nrfrontnodes);
	/* Return nodes */
	for(i=0;i<nrfrontnodes;i++){
		for(dir=0;dir<3;dir++){
			xyz_front[3*i+dir]=xyz_list[3*indicesfront[i]+dir];
		}
	}

	*pxyz_front=xyz_front;

}/*}}}*/
void       Tria::GetLevelsetIntersection(int** pindices, int* pnumiceverts, IssmDouble* fraction, int levelset_enum, IssmDouble level){/*{{{*/

	/* GetLevelsetIntersection computes:
	 * 1. indices of element, sorted in [iceverts, noiceverts] in counterclockwise fashion,
	 * 2. fraction of intersected triangle edges intersected by levelset, lying below level*/

	/*Intermediaries*/
	int numiceverts, numnoiceverts;
	int ind0, ind1, lastindex;
	int indices_ice[NUMVERTICES],indices_noice[NUMVERTICES];
	IssmDouble lsf[NUMVERTICES];
	int* indices = xNew<int>(NUMVERTICES);

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&lsf[0],levelset_enum);

	/* Determine distribution of ice over element.
	 * Exploit: ice/no-ice parts are connected, so find starting vertex of segment*/
	lastindex=0;
	for(int i=0;i<NUMVERTICES;i++){ // go backwards along vertices, and check for sign change
		ind0=(NUMVERTICES-i)%NUMVERTICES;
		ind1=(ind0-1+NUMVERTICES)%NUMVERTICES;
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
	for(int i=0;i<NUMVERTICES;i++){
		ind0=(lastindex+i)%NUMVERTICES;
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
	for(int i=0;i<numiceverts;  i++){
		indices[i]=indices_ice[i];
	}
	for(int i=0;i<numnoiceverts;i++){
		indices[numiceverts+i]=indices_noice[i];
	}

	switch(numiceverts){
		case 0: // no vertex has ice: element is ice free, no intersection
			for(int i=0;i<2;i++)
				fraction[i]=0.;
			break;
		case 1: // one vertex has ice:
			for(int i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[0]])/(lsf[indices[numiceverts+i]]-lsf[indices[0]]);
			}
			break;
		case 2: // two vertices have ice: fraction is computed from first ice vertex to last in CCW fashion
			for(int i=0;i<2;i++){
				fraction[i]=(level-lsf[indices[i]])/(lsf[indices[numiceverts]]-lsf[indices[i]]);
			}
			break;
		case NUMVERTICES: // all vertices have ice: return triangle area
			for(int i=0;i<2;i++)
				fraction[i]=1.;
			break;
		default:
			_error_("Wrong number of ice vertices in Tria::GetLevelsetIntersection!");
			break;
	}

	*pindices=indices;
	*pnumiceverts=numiceverts;
}
/*}}}*/
void       Tria::GetLevelsetPositivePart(int* point1,IssmDouble* fraction1,IssmDouble* fraction2, bool* mainlynegative,IssmDouble* gl){/*{{{*/

	/*Computeportion of the element that has a positive levelset*/

	bool               negative=true;
	int                point=0;
	const IssmPDouble  epsilon= 1.e-15;
	IssmDouble         f1,f2;

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	/*Check that not all nodes are positive or negative*/
	if(gl[0]>0 && gl[1]>0 && gl[2]>0){ // All positive
		point=0;
		f1=1.;
		f2=1.;
	}
	else if(gl[0]<0 && gl[1]<0 && gl[2]<0){ //All negative
		point=0;
		f1=0.;
		f2=0.;
	}
	else{
		if(gl[0]*gl[1]*gl[2]<0) negative=false;

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
		else{
			_error_("This case should NOT be happening");
		}
	}
	*point1=point;
	*fraction1=f1;
	*fraction2=f2;
	*mainlynegative=negative;
}
/*}}}*/
int        Tria::GetVertexIndex(Vertex* vertex){/*{{{*/
	_assert_(vertices);
	for(int i=0;i<NUMVERTICES;i++){
		if(vertex==vertices[i])
		 return i;
	}
	_error_("Vertex provided not found among element vertices");
}
/*}}}*/
int        Tria::GetNumberOfNodes(void){/*{{{*/
	if (this->nodes) return this->NumberofNodes(this->element_type);
	else return 0;
}
/*}}}*/
int        Tria::GetNumberOfNodes(int enum_type){/*{{{*/
	return this->NumberofNodes(enum_type);
}
/*}}}*/
int        Tria::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
void       Tria::GetVectorFromControlInputs(Vector<IssmDouble>* vector,int control_enum,int control_index,int N,const char* data,int offset){/*{{{*/

	/*Get list of ids for this element and this control*/
	int* idlist = xNew<int>(NUMVERTICES*N);
	int  lidlist[NUMVERTICES];
	this->GradientIndexing(&idlist[0],control_index);
	this->GetVerticesLidList(&lidlist[0]);

	/*Get Control*/
	ControlInput* control_input=this->inputs->GetControlInput(control_enum); _assert_(control_input);

	/*Get values on vertices*/
	if(control_input->layout_enum==TriaInputEnum){
		_assert_(N==1);
		TriaInput* input = xDynamicCast<TriaInput*>(control_input->GetInput(data)); _assert_(input);

		if(input->GetInputInterpolationType()==P1Enum){
			IssmDouble values[NUMVERTICES];
			input->Serve(NUMVERTICES,&lidlist[0]);
			for(int i=0;i<NUMVERTICES;i++) values[i] = input->element_values[i];
			vector->SetValues(NUMVERTICES,idlist,values,INS_VAL);
		}
		else if(input->GetInputInterpolationType()==P0Enum){
			input->Serve(1,&this->lid);
			vector->SetValue(idlist[0],input->element_values[0],INS_VAL);
		}
		else{
			_error_("not supported yet");
		}
	}
	else if(control_input->layout_enum==TransientInputEnum){
		TransientInput* input = control_input->GetTransientInput(data); _assert_(input);
		_assert_(N==input->numtimesteps);
		IssmDouble* values = xNew<IssmDouble>(NUMVERTICES*N);

		int count = 0;
		for(int n=0;n<N;n++){
			TriaInput* input_n = input->GetTriaInput(n); _assert_(input_n);
			if(input_n->GetInputInterpolationType()==P1Enum){
				input_n->Serve(NUMVERTICES,&lidlist[0]);
				for(int i=0;i<NUMVERTICES;i++) values[count+i] = input_n->element_values[i];
				count=count+NUMVERTICES;
			}
			else if(input_n->GetInputInterpolationType()==P0Enum){
				input_n->Serve(1,&this->lid);
				values[n] = input_n->element_values[0];
				count++;
			}
			else{
				_error_("not implemented yet");
			}
		}
		vector->SetValues(count,idlist,values,INS_VAL);
		xDelete<IssmDouble>(values);
	}
	else _error_("Type not supported");

	/*Clean up*/
	xDelete<int>(idlist);
}
/*}}}*/
void       Tria::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	int        indices[2];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(2*3);
	this->EdgeOnBaseIndices(&indices[0],&indices[1]);
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
void       Tria::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	int        indices[2];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(2*3);
	this->EdgeOnSurfaceIndices(&indices[0],&indices[1]);
	for(int i=0;i<2;i++) for(int j=0;j<2;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
IssmDouble Tria::GroundedArea(bool scaled){/*{{{*/

	/*Intermediaries*/
	int         domaintype;
	IssmDouble  phi,scalefactor,groundedarea;

	if(!IsIceInElement())return 0.;

	/*Get problem dimension*/
	this->FindParam(&domaintype,DomainTypeEnum);
	if(domaintype!=Domain2DhorizontalEnum && domaintype!=Domain3DEnum) _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	phi=this->GetGroundedPortion(&xyz_list[0][0]);
	groundedarea=phi*this->GetArea();
	if(scaled==true){
		Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
		scalefactor_input->GetInputAverage(&scalefactor);
		groundedarea=groundedarea*scalefactor;
	}

	/*Clean up and return*/
	return groundedarea;
}
/*}}}*/
bool       Tria::HasEdgeOnBase(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on bed not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool       Tria::HasEdgeOnSurface(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on surface not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
IssmDouble Tria::IcefrontMassFlux(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsIcefront()) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	/*Get domain type*/
	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);

	/*Get ice front coordinates*/
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	IssmDouble* xyz_front = NULL;
	this->GetIcefrontCoordinates(&xyz_front,&xyz_list[0][0],MaskIceLevelsetEnum);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],xyz_front);
	//normal[0] = -normal[0];
	//normal[1] = -normal[1];

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],xyz_front,3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		this->JacobianDeterminantSurface(&Jdet,xyz_front,gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(vx*normal[0] + vy*normal[1]);
	}

	/*Cleanup and return*/
	xDelete<IssmDouble>(xyz_front);
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Tria::IcefrontMassFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               domaintype,index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		_error_("not implemented");
	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(vx*normal[0] + vy*normal[1]);
	}

	/*Cleanup and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Tria::GroundinglineMassFlux(bool scaled){/*{{{*/

	/*Make sure there is a grounding line here*/
	if(!IsIceInElement()) return 0;
	if(!IsZeroLevelset(MaskOceanLevelsetEnum)) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               domaintype,index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskOceanLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		_error_("not implemented");
	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],&xyz_front[0][0]);

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble vx,vy,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(vx*normal[0] + vy*normal[1]);
	}

	/*Cleanup and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Tria::IceVolume(bool scaled){/*{{{*/

	/*The volume of a truncated prism is area_base * 1/numedges sum(length of edges)*/

	/*Intermediaries*/
	int i, numiceverts;
	IssmDouble area_base,surface,base,Haverage,scalefactor;
	IssmDouble Haux[NUMVERTICES], surfaces[NUMVERTICES], bases[NUMVERTICES];
	IssmDouble SFaux[NUMVERTICES], scalefactors[NUMVERTICES];
	IssmDouble s[2]; // s:fraction of intersected triangle edges, that lies inside ice
	int* indices=NULL;
	IssmDouble* H=NULL;
	IssmDouble* SF=NULL;

	if(!IsIceInElement())return 0.;

	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);

	/*Relict code
	if(false && IsIcefront()){
		//Assumption: linear ice thickness profile on element.
		//Hence ice thickness at intersection of levelset function with triangle edge is linear interpolation of ice thickness at vertices.
		this->GetLevelsetIntersection(&indices, &numiceverts, s, MaskIceLevelsetEnum, 0.);
		int numthk=numiceverts+2;
		H=xNew<IssmDouble>(numthk);
		//Correct area distortion caused by projection if requestion
		area_base=this->GetAreaIce();
		if(scaled==true){
			Element::GetInputListOnVertices(&scalefactors[0],MeshScaleFactorEnum);
			for(i=0;i<NUMVERTICES;i++) SFaux[i]= scalefactors[indices[i]]; //sort thicknesses in ice/noice
			switch(numiceverts){
				case 1: // average over triangle
					SF[0]=SFaux[0];
					SF[1]=SFaux[0]+s[0]*(SFaux[1]-SFaux[0]);
					SF[2]=SFaux[0]+s[1]*(SFaux[2]-SFaux[0]);
					break;
				case 2: // average over quadrangle
					SF[0]=SFaux[0];
					SF[1]=SFaux[1];
					SF[2]=SFaux[0]+s[0]*(SFaux[2]-SFaux[0]);
					SF[3]=SFaux[1]+s[1]*(SFaux[2]-SFaux[1]);
					break;
				default:
					_error_("Number of ice covered vertices wrong in Tria::IceVolume()");
					break;
			}
			scalefactor=0.;
			for(i=0;i<numthk;i++)	scalefactor+=SF[i];
			scalefactor/=IssmDouble(numthk);
			area_base=area_base*scalefactor;
		}
		Element::GetInputListOnVertices(&surfaces[0],SurfaceEnum);
		Element::GetInputListOnVertices(&bases[0],BaseEnum);
		for(i=0;i<NUMVERTICES;i++) Haux[i]= surfaces[indices[i]]-bases[indices[i]]; //sort thicknesses in ice/noice
		switch(numiceverts){
			case 1: // average over triangle
				H[0]=Haux[0];
				H[1]=Haux[0]+s[0]*(Haux[1]-Haux[0]);
				H[2]=Haux[0]+s[1]*(Haux[2]-Haux[0]);
				break;
			case 2: // average over quadrangle
				H[0]=Haux[0];
				H[1]=Haux[1];
				H[2]=Haux[0]+s[0]*(Haux[2]-Haux[0]);
				H[3]=Haux[1]+s[1]*(Haux[2]-Haux[1]);
				break;
			default:
				_error_("Number of ice covered vertices wrong in Tria::IceVolume()");
				break;
		}
		Haverage=0.;
		for(i=0;i<numthk;i++)	Haverage+=H[i];
		Haverage/=IssmDouble(numthk);
	}*/

   IssmDouble lsf[NUMVERTICES];
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
   /*Deal with partially ice-covered elements*/
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		bool istrapneg;
      int point;
      IssmDouble weights[NUMVERTICES];
      IssmDouble surfaces[NUMVERTICES];
      IssmDouble bases[NUMVERTICES];
      IssmDouble Hice[NUMVERTICES];
      IssmDouble area_basetot,f1,f2,phi;
      /*Average thickness over subelement*/
		Element::GetInputListOnVertices(&surfaces[0],SurfaceEnum);
      Element::GetInputListOnVertices(&bases[0],BaseEnum);
      GetFractionGeometry(weights,&phi,&point,&f1,&f2,&istrapneg,lsf);
      for(int i=0;i<NUMVERTICES;i++) Hice[i] = surfaces[i]-bases[i];
      Haverage = 0.0;
		/*Use weights[i]/phi to get average thickness over subelement*/
      for(int i=0;i<NUMVERTICES;i++) Haverage += weights[i]/phi*Hice[i];
		/*Get back area of ice-covered base*/
		area_basetot = this->GetArea();
		area_base    = phi*area_basetot;

		/*Account for scaling factor averaged over subelement*/
		if(scaled==true){
			IssmDouble* scalefactor_vertices = xNew<IssmDouble>(NUMVERTICES);
			Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
			scalefactor = 0.0;
			for(int i=0;i<NUMVERTICES;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
			area_base = area_base*scalefactor;
			xDelete<IssmDouble>(scalefactor_vertices);
		}
   }
	else{
		/*First get back the area of the base*/
		area_base=this->GetArea();
		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);
			area_base=area_base*scalefactor;
		}

		/*Now get the average height*/
		Input* surface_input = this->GetInput(SurfaceEnum); _assert_(surface_input);
		Input* base_input    = this->GetInput(BaseEnum);    _assert_(base_input);
		surface_input->GetInputAverage(&surface);
		base_input->GetInputAverage(&base);
		Haverage=surface-base;
	}

	/*Cleanup & return: */
	xDelete<int>(indices);
	xDelete<IssmDouble>(H);
	xDelete<IssmDouble>(SF);

	if(domaintype==Domain2DverticalEnum){
	  return area_base;
	}
	else{
	  return area_base*Haverage;
	}
}
/*}}}*/
IssmDouble Tria::IceVolumeAboveFloatation(bool scaled){/*{{{*/

	/*The volume above floatation: H + rho_water/rho_ice * bathymetry */
	IssmDouble rho_ice,rho_water,vaf,HAFaverage;
	IssmDouble area_base,surface,bed,bathymetry,scalefactor;
	IssmDouble xyz_list[NUMVERTICES][3];
   IssmDouble lsf[NUMVERTICES];

	if(!IsIceInElement() || IsAllFloating())return 0;

	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);

	rho_ice=FindParam(MaterialsRhoIceEnum);
	rho_water=FindParam(MaterialsRhoSeawaterEnum);
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		bool istrapneg;
      int point;
      IssmDouble weights[NUMVERTICES];
      IssmDouble surfaces[NUMVERTICES];
      IssmDouble bases[NUMVERTICES];
      IssmDouble bathys[NUMVERTICES];
      IssmDouble HAF[NUMVERTICES];
      IssmDouble area_basetot,f1,f2,phi;
      /*Average thickness over subelement*/
		Element::GetInputListOnVertices(&surfaces[0],SurfaceEnum);
      Element::GetInputListOnVertices(&bases[0],BaseEnum);
      Element::GetInputListOnVertices(&bathys[0],BedEnum);
      GetFractionGeometry(weights,&phi,&point,&f1,&f2,&istrapneg,lsf);
      for(int i=0;i<NUMVERTICES;i++) HAF[i] = surfaces[i]-bases[i]+min(rho_water/rho_ice*bathys[i],0.);
      HAFaverage = 0.0;
		/*Use weights[i]/phi to get average thickness over subelement*/
      for(int i=0;i<NUMVERTICES;i++) HAFaverage += weights[i]/phi*HAF[i];
		/*Get back area of ice-covered base*/
		area_basetot = this->GetArea();
		area_base    = phi*area_basetot;

		/*Account for scaling factor averaged over subelement*/
		if(scaled==true){
			IssmDouble* scalefactor_vertices = xNew<IssmDouble>(NUMVERTICES);
			Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
			scalefactor = 0.0;
			for(int i=0;i<NUMVERTICES;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
			area_base = area_base*scalefactor;
			xDelete<IssmDouble>(scalefactor_vertices);
		}
		vaf = area_base*HAFaverage;
   }
	else{
		/*First calculate the area of the base (cross section triangle)
		 * http://en.wikipedia.org/wiki/Triangle
		 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
		area_base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));
		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);
			area_base=area_base*scalefactor;
		}

		/*Now get the average height and bathymetry*/
		Input* surface_input = this->GetInput(SurfaceEnum); _assert_(surface_input);
		Input* base_input    = this->GetInput(BaseEnum);    _assert_(base_input);
		Input* bed_input     = this->GetInput(BedEnum);     _assert_(bed_input);
		if(!bed_input) _error_("Could not find bed");
		surface_input->GetInputAverage(&surface);
		base_input->GetInputAverage(&bed);
		bed_input->GetInputAverage(&bathymetry);
		vaf = area_base*(surface-bed+min(rho_water/rho_ice*bathymetry,0.));
	}

	/*Return: */
	return vaf;
}
/*}}}*/
void       Tria::InputDepthAverageAtBase(int enum_type,int average_enum_type){/*{{{*/

	/*New input*/
	Input* oldinput=NULL;
	Input* newinput=NULL;

	/*copy input of enum_type*/
	oldinput=this->GetInput(enum_type);
	if(!oldinput)_error_("could not find old input with enum: " << EnumToStringx(enum_type));
	newinput=oldinput->copy();

	/*Assign new name (average)*/
	newinput->ChangeEnum(average_enum_type);

	/*Add new input to current element*/
	_error_("not implemented");
}
/*}}}*/
void       Tria::InputUpdateFromIoModel(int index, IoModel* iomodel){ //i is the element index/*{{{*/

	/*Intermediaries*/
	int        i,j;
	int        tria_vertex_ids[3];
	IssmDouble nodeinputs[3];
	IssmDouble cmmininputs[3];
	IssmDouble cmmaxinputs[3];
	bool       control_analysis,ad_analysis   = false;
	int        num_control_type,num_responses;
	char**     controls = NULL;
	IssmDouble yts;

	/*Get parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&control_analysis,"md.inversion.iscontrol");
	iomodel->FindConstant(&ad_analysis, "md.autodiff.isautodiff");
	if(control_analysis && !ad_analysis) iomodel->FindConstant(&num_control_type,"md.inversion.num_control_parameters");
	if(control_analysis && !ad_analysis) iomodel->FindConstant(&num_responses,"md.inversion.num_cost_functions");
	if(control_analysis && ad_analysis) iomodel->FindConstant(&num_control_type,"md.autodiff.num_independent_objects");
	if(control_analysis && ad_analysis) iomodel->FindConstant(&num_responses,"md.autodiff.num_dependent_objects");

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<3;i++){
		tria_vertex_ids[i]=reCast<int>(iomodel->elements[3*index+i]); //ids for vertices are in the elements array from Matlab
	}
}
/*}}}*/
void       Tria::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

	/*Intermediary*/
	int* doflist = NULL;

	/*Fetch number of nodes for this finite element*/
	int numnodes = this->NumberofNodes(this->element_type);

	/*Fetch dof list and allocate solution vector*/
	GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
	IssmDouble* values    = xNew<IssmDouble>(numnodes);

	/*Use the dof list to index into the solution vector: */
	for(int i=0;i<numnodes;i++){
		values[i]=solution[doflist[i]];
		if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in solution vector");
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector, SID = " << this->Sid());
	}

	/*Add input to the element: */
	this->AddInput(enum_type,values,this->element_type);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
void       Tria::InputUpdateFromVector(IssmDouble* vector, int name, int type){/*{{{*/

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	int         numnodes;
	IssmDouble  value;
	int         lidlist[NUMVERTICES];
	int        *doflist = NULL;
	IssmDouble *values  = NULL;

	GetVerticesLidList(&lidlist[0]);

	switch(type){
		case VertexLIdEnum:
			values = xNew<IssmDouble>(NUMVERTICES);
			for(int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Lid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
			}
			/*update input*/
			inputs->SetTriaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			break;

		case VertexPIdEnum:
			values = xNew<IssmDouble>(NUMVERTICES);
			for(int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Pid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
			}
			/*update input*/
			inputs->SetTriaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			break;

		case VertexSIdEnum:
			values = xNew<IssmDouble>(NUMVERTICES);
			for(int i=0;i<NUMVERTICES;i++){
				values[i]=vector[this->vertices[i]->Sid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
			}
			/*update input*/
			inputs->SetTriaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			break;

		case NodesEnum:
			/*Get number of nodes and dof list: */
			numnodes = this->NumberofNodes(this->element_type);
			values   = xNew<IssmDouble>(numnodes);
			GetDofList(&doflist,NoneApproximationEnum,GsetEnum);

			for(int i=0;i<numnodes;i++){
				values[i]=vector[doflist[i]];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
			}
			//this->inputs->AddInput(new TriaInput(name,values,this->element_type));
			_error_("not implemented");
			break;

		case NodeSIdEnum:
			/*Get number of nodes and dof list: */
			numnodes = this->NumberofNodes(this->element_type);
			values   = xNew<IssmDouble>(numnodes);

			for(int i=0;i<numnodes;i++){
				values[i]=vector[nodes[i]->Sid()];
				if(xIsNan<IssmDouble>(values[i])) _error_("NaN found in vector");
				if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in vector");
			}
			if(this->element_type==P1Enum){
				inputs->SetTriaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			}
			else{
				inputs->SetTriaInput(name,this->element_type,this->lid,numnodes,values);
			}
			break;

		case ElementEnum:
			value=vector[this->Sid()];
			if(xIsNan<IssmDouble>(value)) _error_("NaN found in vector");
			if(xIsInf<IssmDouble>(value)) _error_("Inf found in vector");
			/*update input*/
			//this->inputs->AddInput(new DoubleInput(name,value));
			//inputs->SetTriaInput(name,P1Enum,NUMVERTICES,lidlist,values);
			_error_("not implemented");
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}

	/*Clean-up*/
	xDelete<int>(doflist);
	xDelete<IssmDouble>(values);

}
/*}}}*/
bool       Tria::IsFaceOnBoundary(void){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonboundaryEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2.);

	if(sum==3.)  _error_("Two edges on boundary not supported yet...");

	if(sum>1.){
		return true;
	}
	else{
		return false;
	}
}/*}}}*/
bool       Tria::IsIcefront(void){/*{{{*/

	bool isicefront;
	int i,nrice;
   IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	isicefront=false;
	if(IsIceInElement()){
		nrice=0;
		for(i=0;i<NUMVERTICES;i++)
			if(ls[i]<0.) nrice++;
		if(nrice==1) isicefront= true;
	}
	return isicefront;
}/*}}}*/
bool       Tria::IsNodeOnShelfFromFlags(IssmDouble* flags){/*{{{*/

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
bool       Tria::IsZeroLevelset(int levelset_enum){/*{{{*/

	bool iszerols;
	IssmDouble ls[NUMVERTICES];

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&ls[0],levelset_enum);

	/*If the level set is awlays <0, there is no ice front here*/
	iszerols= false;
	if(IsIceInElement()){
		if(ls[0]*ls[1]<0. || ls[0]*ls[2]<0. || (ls[0]*ls[1]*ls[2]==0. && ls[0]*ls[1]+ls[0]*ls[2]+ls[1]*ls[2]<=0.)){
			iszerols = true;
		}
	}

	return iszerols;
}
/*}}}*/
void       Tria::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_base,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list,(GaussTria*)gauss);

}
/*}}}*/
void       Tria::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_top,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetSegmentJacobianDeterminant(pJdet,xyz_list_top,(GaussTria*)gauss);

}
/*}}}*/
IssmDouble Tria::Masscon(IssmDouble* levelset){ /*{{{*/

	/*intermediary: */
	IssmDouble  thickness;
	IssmDouble  Jdet;
	int         point1;
	IssmDouble  fraction1,fraction2;
	bool        mainlynegative=true;

	/* Get node coordinates and dof list: */
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve inputs required:*/
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);

	/*Retrieve material parameters: */
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);

	/*Retrieve values of the levelset defining the masscon: */
	IssmDouble values[NUMVERTICES];
	for(int i=0;i<NUMVERTICES;i++){
		values[i]=levelset[this->vertices[i]->Sid()];
	}

	/*Ok, use the level set values to figure out where we put our gaussian points:*/
	this->GetLevelsetPositivePart(&point1,&fraction1,&fraction2,&mainlynegative,&values[0]);
	Gauss* gauss = this->NewGauss(point1,fraction1,fraction2,mainlynegative,4);

	IssmDouble volume=0.;
	while(gauss->next()){
		this->JacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		thickness_input->GetInputValue(&thickness, gauss);
		volume+=thickness*gauss->weight*Jdet;
	}

	/* clean up and Return: */
	delete gauss;
	return rho_ice*volume;
}
/*}}}*/
IssmDouble Tria::MassFlux(IssmDouble x1, IssmDouble y1, IssmDouble x2, IssmDouble y2,int segment_id){/*{{{*/

	int        domaintype;
	IssmDouble mass_flux=0.;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble vx1,vx2,vy1,vy2,h1,h2;

	/*Get material parameters :*/
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);

	/*First off, check that this segment belongs to this element: */
	if (segment_id!=this->id)_error_("error message: segment with id " << segment_id << " does not belong to element with id:" << this->id);

	/*Get xyz list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*get area coordinates of 0 and 1 locations: */
	GaussTria gauss_1;
	gauss_1.GaussFromCoords(x1,y1,&xyz_list[0][0]);
	GaussTria gauss_2;
	gauss_2.GaussFromCoords(x2,y2,&xyz_list[0][0]);

	/*Get segment length and normal (needs to be properly oriented)*/
	IssmDouble nx=cos(atan2(x1-x2,y2-y1));
	IssmDouble ny=sin(atan2(x1-x2,y2-y1));
	IssmDouble length=sqrt(pow(x2-x1,2)+pow(y2-y1,2));

	/*Get velocity and thickness*/
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	Input* thickness_input=this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* vx_input=NULL;
	Input* vy_input=NULL;
	if(domaintype==Domain2DhorizontalEnum){
		vx_input=this->GetInput(VxEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyEnum); _assert_(vy_input);
	}
	else{
		vx_input=this->GetInput(VxAverageEnum); _assert_(vx_input);
		vy_input=this->GetInput(VyAverageEnum); _assert_(vy_input);
	}

	thickness_input->GetInputValue(&h1, &gauss_1);
	thickness_input->GetInputValue(&h2, &gauss_2);
	vx_input->GetInputValue(&vx1,&gauss_1);
	vx_input->GetInputValue(&vx2,&gauss_2);
	vy_input->GetInputValue(&vy1,&gauss_1);
	vy_input->GetInputValue(&vy2,&gauss_2);

	mass_flux= rho_ice*length*(
				(1./3.*(h1-h2)*(vx1-vx2)+0.5*h2*(vx1-vx2)+0.5*(h1-h2)*vx2+h2*vx2)*nx+
				(1./3.*(h1-h2)*(vy1-vy2)+0.5*h2*(vy1-vy2)+0.5*(h1-h2)*vy2+h2*vy2)*ny
				);

	/*clean up and return:*/
	return mass_flux;
}
/*}}}*/
IssmDouble Tria::MassFlux(IssmDouble* segment){/*{{{*/
	return this->MassFlux(segment[0],segment[1],segment[2],segment[3],reCast<int>(segment[4]));
}
/*}}}*/
IssmDouble Tria::Misfit(int modelenum,int observationenum,int weightsenum){/*{{{*/

	/*Intermediaries*/
	IssmDouble model,observation,weight;
	IssmDouble Jdet;
	IssmDouble Jelem = 0;
	IssmDouble xyz_list[NUMVERTICES][3];
	GaussTria *gauss = NULL;

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Retrieve all inputs we will be needing: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* model_input=this->GetInput(modelenum);   _assert_(model_input);
	Input* observation_input=this->GetInput(observationenum);_assert_(observation_input);
	Input* weights_input     =this->GetInput(weightsenum);     _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		GetJacobianDeterminant(&Jdet, &xyz_list[0][0],gauss);

		/*Get parameters at gauss point*/
		model_input->GetInputValue(&model,gauss);
		observation_input->GetInputValue(&observation,gauss);
		weights_input->GetInputValue(&weight,gauss);

		/*compute misfit between model and observation */
		Jelem+=0.5*(model-observation)*(model-observation)*Jdet*weight*gauss->weight;
	}

	/* clean up and Return: */
	delete gauss;
	return Jelem;
}
/*}}}*/
IssmDouble Tria::MisfitArea(int weightsenum){/*{{{*/

	/*Intermediaries*/
	IssmDouble weight;
	IssmDouble Jdet;
	IssmDouble Jelem = 0;
	IssmDouble xyz_list[NUMVERTICES][3];
	GaussTria *gauss = NULL;

	/*If on water, return 0: */
	if(!IsIceInElement())return 0;

	/*Retrieve all inputs we will be needing: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* weights_input     =this->GetInput(weightsenum);     _assert_(weights_input);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(2);
	while(gauss->next()){

		/* Get Jacobian determinant: */
		GetJacobianDeterminant(&Jdet, &xyz_list[0][0],gauss);

		/*Get parameters at gauss point*/
		weights_input->GetInputValue(&weight,gauss);

		/*compute misfit between model and observation */
		Jelem+=Jdet*weight*gauss->weight;
	}

	/* clean up and Return: */
	delete gauss;
	return Jelem;
}
/*}}}*/
void	      Tria::MovingFrontalVelocity(void){/*{{{*/

	int  dim, domaintype, calvinglaw, i;
	IssmDouble v[3],w[3],c[3],m[3],dlsf[3];
	IssmDouble norm_dlsf, norm_calving, calvingrate, meltingrate, groundedice;
	IssmDouble migrationmax, calvinghaf, heaviside, haf_eps;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble movingfrontvx[NUMVERTICES];
	IssmDouble movingfrontvy[NUMVERTICES];
	IssmDouble vel;

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	Input* calvingratex_input = NULL;
	Input* calvingratey_input = NULL;
	Input* calvingrate_input  = NULL;
	Input* meltingrate_input  = NULL;

	/*Get problem dimension and whether there is moving front or not*/
	this->FindParam(&domaintype,DomainTypeEnum);
	this->FindParam(&calvinglaw,CalvingLawEnum);

	switch(domaintype){
		case Domain2DverticalEnum:   dim = 1; break;
		case Domain2DhorizontalEnum: dim = 2; break;
		case Domain3DEnum:           dim = 2; break;
		default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
	}

	/*Load levelset function gradients*/
   Input *gr_input         = this->GetInput(MaskOceanLevelsetEnum); _assert_(gr_input);
	Input *lsf_slopex_input = this->GetInput(LevelsetfunctionSlopeXEnum); _assert_(lsf_slopex_input);
	Input *lsf_slopey_input = this->GetInput(LevelsetfunctionSlopeYEnum); _assert_(lsf_slopey_input);
	Input *vx_input         = this->GetInput(VxEnum);                     _assert_(vx_input);
	Input *vy_input         = this->GetInput(VyEnum); _assert_(vy_input);

	switch(calvinglaw){
		case DefaultCalvingEnum:
		case CalvingVonmisesEnum:
		case CalvingVonmisesADEnum:
		case CalvingLevermannEnum:
		case CalvingPollardEnum:
		case CalvingTestEnum:
		case CalvingParameterizationEnum:
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
	GaussTria gauss;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Advection */
		vx_input->GetInputValue(&v[0],&gauss);
		vy_input->GetInputValue(&v[1],&gauss);
		gr_input->GetInputValue(&groundedice,&gauss);
		lsf_slopex_input->GetInputValue(&dlsf[0],&gauss);
		if(dim==2) lsf_slopey_input->GetInputValue(&dlsf[1],&gauss);

		norm_dlsf=0.;
		for(i=0;i<dim;i++) norm_dlsf+=pow(dlsf[i],2);
		norm_dlsf=sqrt(norm_dlsf);

		/*Get calving speed*/
		switch(calvinglaw){

			/*RATE calving laws*/
			case DefaultCalvingEnum:
			case CalvingVonmisesEnum:
			case CalvingVonmisesADEnum:
			case CalvingTestEnum:
			case CalvingParameterizationEnum:
			case CalvingLevermannEnum:
			case CalvingPollardEnum:
			case CalvingCalvingMIPEnum:
				calvingratex_input->GetInputValue(&c[0],&gauss);
				calvingratey_input->GetInputValue(&c[1],&gauss);
				meltingrate_input->GetInputValue(&meltingrate,&gauss);
				if(groundedice<0) meltingrate = 0.;
				m[0]=meltingrate*dlsf[0]/norm_dlsf;
				m[1]=meltingrate*dlsf[1]/norm_dlsf;

				if(norm_dlsf<1.e-10){
					for(i=0;i<dim;i++){
						c[i]=0.; m[i]=0.;
					}
				}
				break;

			/*Discrete calving laws*/
			case CalvingMinthicknessEnum:
			case CalvingHabEnum:
			case CalvingCrevasseDepthEnum:
				meltingrate_input->GetInputValue(&meltingrate,&gauss);

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
						heaviside=(groundedice-calvinghaf+haf_eps)/(2.*haf_eps) + sin(PI*(groundedice-calvinghaf)/haf_eps)/(2.*M_PI);
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
				calvingratex_input->GetInputValue(&c[0],&gauss);
				if(dim==2) calvingratey_input->GetInputValue(&c[1],&gauss);
				for(i=0;i<dim;i++) m[i]=0.;
				break;

			default:
				_error_("Calving law "<<EnumToStringx(calvinglaw)<<" not supported yet");
		}
		for(i=0;i<dim;i++) w[i]=v[i]-c[i]-m[i];

		movingfrontvx[iv] = w[0];
		movingfrontvy[iv] = w[1];		
	}

	#ifdef MICI
	/**************************************  MICI  START ************************************/
	/*MICI from Crawford et al. 2021:
	 * Tice = 20C; Bf  ->  I = 3.7e-16; α = 6.9
	 * Tice = 20C; Bn  ->  I = 5.1e-14; α = 6.0
	 * Tice = 20C; Bh  ->  I = 3.2e-17; α = 7.2
	 * Tice = 10C; Bn  ->  I = 6.9e-17; α = 7.3
	 * Tice = 5C;  Bn  ->  I = 1.9e-16; α = 7.3
	 */
	Input* surface_input = this->GetInput(SurfaceEnum); _assert_(surface_input);
	Input* bed_input     = this->GetInput(BedEnum);     _assert_(bed_input);
	Input* ls_input      = this->GetInput(MaskIceLevelsetEnum);   _assert_(ls_input);
	IssmDouble Hc,bed,ls;
	for(int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		surface_input->GetInputValue(&Hc,&gauss);
		bed_input->GetInputValue(&bed,&gauss);
		ls_input->GetInputValue(&ls,&gauss);

		/*use vel direction instead of LSF*/
		vx_input->GetInputValue(&v[0],&gauss);
		vy_input->GetInputValue(&v[1],&gauss);
		vel=sqrt(v[0]*v[0] + v[1]*v[1]);
		norm_dlsf = max(vel,1.e-10);
		dlsf[0] = v[0];
		dlsf[1] = v[1];

		/*Do we assume that the calving front does not move if MICI is not engaged?*/
		bool regrowth = false;
		bool apply_as_retreat = true;
		if(!regrowth){
			movingfrontvx[iv] = 0.;
			movingfrontvy[iv] = 0.;
		}

		//movingfrontvx[iv] = -2000./(365*24*3600.)*dlsf[0]/norm_dlsf;
		//movingfrontvy[iv] = -2000./(365*24*3600.)*dlsf[1]/norm_dlsf;

		if(MICI==1 && Hc>80. && bed<0. && fabs(ls)<100.e3){ //Pollard & De Conto
			IssmDouble C = (min(max(Hc,80.),100.) - 80.)/20. * 10./(24*3600.); /*Original MICI! convert from m/day to m/s*/

			/*Front motion = ice speed (v) - calving rate*/
			movingfrontvx[iv] = v[0]-C*dlsf[0]/norm_dlsf;
			movingfrontvy[iv] = v[1]-C*dlsf[1]/norm_dlsf;
		}
		else if (MICI==2 && Hc>135. && bed<0. && fabs(ls)<100.e3){ // Crawford et all

			/*if 1: RETREAT rate
			 *if 0: calving rate*/
			if(0) v[0]=0.; v[1]=0.;

			/*5C Bn (worst case scenario)*/
			IssmDouble I     = 1.9e-16;
			IssmDouble alpha = 7.3;
			IssmDouble C = min(2000.,I*pow(Hc,alpha))/(24*3600.); /*convert from m/day to m/s*/

			/*Front motion = ice speed (v) - calving rate*/
			movingfrontvx[iv] = v[0] -C*dlsf[0]/norm_dlsf;
			movingfrontvy[iv] = v[1] -C*dlsf[1]/norm_dlsf;

			/*disable regrowth if calving rate is too low*/
			if(!regrowth && C<vel){
				movingfrontvx[iv] = 0.;
				movingfrontvy[iv] = 0.;
			}
		}
	}
	#endif
	/**************************************  END MICI  *************************************/

	/*Add input*/
	this->AddInput(MovingFrontalVxEnum,&movingfrontvx[0],P1DGEnum);
	this->AddInput(MovingFrontalVyEnum,&movingfrontvy[0],P1DGEnum);
}
/*}}}*/
Gauss*     Tria::NewGauss(void){/*{{{*/
	return new GaussTria();
}
/*}}}*/
Gauss*     Tria::NewGauss(int order){/*{{{*/
	return new GaussTria(order);
}
/*}}}*/
Gauss*     Tria::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order){/*{{{*/

	IssmDouble  area_coordinates[2][3];
	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,2);
	return new GaussTria(area_coordinates,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,bool mainlyfloating,int order){/*{{{*/

	return new GaussTria(point1,fraction1,fraction2,mainlyfloating,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(int point1,IssmDouble fraction1,IssmDouble fraction2,int order){/*{{{*/

	return new GaussTria(point1,fraction1,fraction2,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(IssmDouble fraction1,IssmDouble fraction2,int order){/*{{{*/

	return new GaussTria(fraction1,fraction2,order);
}
/*}}}*/
Gauss*     Tria::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/

	IssmDouble  area_coordinates[2][3];
	GetAreaCoordinates(&area_coordinates[0][0],xyz_list_front,xyz_list,2);
	return new GaussTria(area_coordinates,order_vert);
}
/*}}}*/
Gauss*     Tria::NewGaussBase(int order){/*{{{*/

	int indices[2];
	this->EdgeOnBaseIndices(&indices[0],&indices[1]);
	return new GaussTria(indices[0],indices[1],order);
}
/*}}}*/
Gauss*     Tria::NewGaussTop(int order){/*{{{*/

	int indices[2];
	this->EdgeOnSurfaceIndices(&indices[0],&indices[1]);
	return new GaussTria(indices[0],indices[1],order);
}
/*}}}*/
void       Tria::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->element_type);

}
/*}}}*/
void       Tria::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,this->element_type);

}
/*}}}*/
void       Tria::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->PressureInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsP1(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,P1Enum);

}
/*}}}*/
void       Tria::NodalFunctionsP1Derivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTria*)gauss,P1Enum);

}
/*}}}*/
void       Tria::NodalFunctionsP2(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,P2Enum);

}
/*}}}*/
void       Tria::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->TensorInterpolation());

}
/*}}}*/
void       Tria::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTriaEnum);
	this->GetNodalFunctions(basis,(GaussTria*)gauss,this->VelocityInterpolation());

}
/*}}}*/
int        Tria::NodalValue(IssmDouble* pvalue, int index, int natureofdataenum){/*{{{*/

	int         found = 0;
	IssmDouble  value;
	GaussTria   gauss;

	/*First, serarch the input: */
	Input* data=this->GetInput(natureofdataenum);

	/*figure out if we have the vertex id: */
	found=0;
	for(int i=0;i<NUMVERTICES;i++){
		if(index==vertices[i]->Sid()){
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
void       Tria::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/
	LineSectionNormal(bed_normal, xyz_list);
	_assert_(bed_normal[1]<0);
}
/*}}}*/
void       Tria::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/
	LineSectionNormal(normal, xyz_list);
}
/*}}}*/
void       Tria::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/
	LineSectionNormal(top_normal, xyz_list);
	_assert_(top_normal[1]>0);
}
/*}}}*/
int        Tria::ObjectEnum(void){/*{{{*/

	return TriaEnum;

}
/*}}}*/
int        Tria::NumberofNodesPressure(void){/*{{{*/
	return TriaRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int        Tria::NumberofNodesVelocity(void){/*{{{*/
	return TriaRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
void       Tria::PotentialUngrounding(Vector<IssmDouble>* potential_ungrounding){/*{{{*/

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

	/*go through vertices, and figure out which ones are grounded and want to unground: */
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
int        Tria::PressureInterpolation(void){/*{{{*/
	return TriaRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void       Tria::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	/*Static condensation if requested*/
	if(pe){
		if(this->element_type==MINIcondensedEnum){
			int indices[2]={6,7};
			pe->StaticCondensation(Ke,2,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[3]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<3;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			pe->StaticCondensation(Ke,size,indices);
			xDelete<int>(indices);
		}
	}

	if(Ke){
		if(this->element_type==MINIcondensedEnum){
			int indices[2]={6,7};
			Ke->StaticCondensation(2,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[3]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<3;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			Ke->StaticCondensation(size,indices);
			xDelete<int>(indices);
		}
	}

}
/*}}}*/
void       Tria::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int numnodes = this->NumberofNodesVelocity();

	int          approximation;
	IssmDouble*  vertexonbase= NULL;
	IssmDouble   slope,groundedice;
	IssmDouble   xz_plane[6];

	/*For FS only: we want the CS to be tangential to the bedrock*/
	this->Element::GetInputValue(&approximation,ApproximationEnum);
	if(!HasNodeOnBase() ||  approximation!=FSApproximationEnum) return;

	/*Get inputs*/
	Input* slope_input=this->GetInput(BedSlopeXEnum);                             _assert_(slope_input);
	Input* groundedicelevelset_input=this->GetInput(MaskOceanLevelsetEnum); _assert_(groundedicelevelset_input);
	vertexonbase = xNew<IssmDouble>(numnodes);
	this->GetInputListOnNodesVelocity(&vertexonbase[0],MeshVertexonbaseEnum);

	/*Loop over basal nodes and update their CS*/
	GaussTria gauss;
	for(int i=0;i<this->NumberofNodesVelocity();i++){

		if(vertexonbase[i]==1){
			gauss.GaussNode(this->VelocityInterpolation(),i);
			slope_input->GetInputValue(&slope,&gauss);
			groundedicelevelset_input->GetInputValue(&groundedice,&gauss);
			IssmDouble theta = atan(slope);

			/*New X axis                  New Z axis*/
			xz_plane[0]=cos(theta);       xz_plane[3]=0.;
			xz_plane[1]=sin(theta);       xz_plane[4]=0.;
			xz_plane[2]=0.;               xz_plane[5]=1.;

			if(groundedice>=0){
				this->nodes[i]->DofInSSet(1); //vy
			}
			else{
				this->nodes[i]->DofInFSet(1); //vy
			}

			XZvectorsToCoordinateSystem(&this->nodes[i]->coord_system[0][0],&xz_plane[0]);
			this->nodes[i]->isrotated = true;
		}
	}

	/*cleanup*/
	xDelete<IssmDouble>(vertexonbase);
}
/*}}}*/
void       Tria::ResetHooks(){/*{{{*/

	if(this->nodes) xDelete<Node*>(this->nodes);
	this->nodes=NULL;
	this->vertices=NULL;
	this->material=NULL;
	this->parameters=NULL;

	//deal with ElementHook mother class
	for(int i=0;i<this->numanalyses;i++) if(this->hnodes[i]) this->hnodes[i]->reset();
	this->hvertices->reset();
	if(this->hmaterial)this->hmaterial->reset();
	if(this->hneighbors) this->hneighbors->reset();

}
/*}}}*/
void       Tria::SetControlInputsFromVector(IssmDouble* vector,int control_enum,int control_index,int offset,int M, int N){/*{{{*/

	IssmDouble  values[NUMVERTICES];
	int         lidlist[NUMVERTICES];

	/*Get Domain type*/
	int domaintype;
	parameters->FindParam(&domaintype,DomainTypeEnum);

	/*Specific case for depth averaged quantities*/
	int control_init=control_enum;
	if(domaintype==Domain2DverticalEnum){
		if(control_enum==MaterialsRheologyBbarEnum){
			control_enum=MaterialsRheologyBEnum;
			if(!IsOnBase()) return;
		}
		if(control_enum==DamageDbarEnum){
			control_enum=DamageDEnum;
			if(!IsOnBase()) return;
		}
	}

	/*Get list of ids for this element and this control*/
	int* idlist = xNew<int>(NUMVERTICES*N);
	GradientIndexing(&idlist[0],control_index);

	ControlInput* control_input=this->inputs->GetControlInput(control_enum); _assert_(control_input);
	this->GetVerticesLidList(&lidlist[0]);

	if(control_input->layout_enum==TriaInputEnum){
		ElementInput* input = control_input->GetInput("value"); _assert_(input);
		if(input->GetInputInterpolationType()==P1Enum){
			_assert_(N==1);
			for(int i=0;i<NUMVERTICES;i++){
				values[i] = vector[idlist[i]];
			}
			input->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
		}
		else if(input->GetInputInterpolationType()==P0Enum){
			_assert_(N==1);
			input->SetInput(P0Enum,this->lid,vector[idlist[0]]);
		}
		else{
			_error_("not implemented yet");
		}
	}
	else if(control_input->layout_enum==TransientInputEnum){
		_assert_(N>1);
		TransientInput* input = control_input->GetTransientInput("value"); _assert_(input);

		int count = 0;
		for(int n=0;n<N;n++){
			TriaInput* input_n = input->GetTriaInput(n); _assert_(input_n);
			if(input_n->GetInputInterpolationType()==P1Enum){
				for(int i=0;i<NUMVERTICES;i++){
					values[i] = vector[idlist[count]];
					count++;
				}
				input_n->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
			}
			else if(input_n->GetInputInterpolationType()==P0Enum){
				input_n->SetInput(P0Enum,this->lid,vector[idlist[count]]);
				count++;
			}
			else{
				_error_("not implemented yet");
			}
		}
	}
	else _error_("Type not supported");

	/*Special handling of some inputs*/
	if(control_enum==BedEnum){
		IssmDouble bedlist[NUMVERTICES];
		IssmDouble thicknesslist[NUMVERTICES];
		IssmDouble surfacelist[NUMVERTICES];
		IssmDouble phi[NUMVERTICES];
		Element::GetInputListOnVertices(&bedlist[0], BedEnum);
		Element::GetInputListOnVertices(&thicknesslist[0], ThicknessEnum);
		Element::GetInputListOnVertices(&surfacelist[0], SurfaceEnum);
		Element::GetInputListOnVertices(&phi[0],MaskOceanLevelsetEnum);

		for(int i=0;i<NUMVERTICES;i++) {
			if(phi[i]>0.){
				thicknesslist[i] = surfacelist[i]-bedlist[i]; //surface = oldbase + newthickness
			}
			if(thicknesslist[i]<1.){
				thicknesslist[i]=1.;
			}
		}

		/*Add input to the element: */
		this->AddBasalInput(ThicknessEnum,thicknesslist,P1Enum);
	}

	/*Clean up*/
	xDelete<int>(idlist);
}
/*}}}*/
void       Tria::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

	/*go into parameters and get the analysis_counter: */
	int analysis_counter;
	parametersin->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Get Element type*/
	if(this->element_type_list) this->element_type=this->element_type_list[analysis_counter];

	/*Pick up nodes*/
	if(this->hnodes && this->hnodes[analysis_counter]){
		this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	}

}
/*}}}*/
void       Tria::SetElementInput(int enum_in,IssmDouble value){/*{{{*/

	this->SetElementInput(this->inputs,enum_in,value);

}
/*}}}*/
void       Tria::SetElementInput(Inputs* inputs,int enum_in,IssmDouble value){/*{{{*/

	_assert_(inputs);
	inputs->SetTriaInput(enum_in,P0Enum,this->lid,value);

}
/*}}}*/
void       Tria::SetElementInput(int enum_in,IssmDouble value,int type){/*{{{*/

	if(type==P0Enum){
		this->inputs->SetTriaInput(enum_in,P0Enum,this->lid,value);
	}
	else if(type==P1Enum){
		IssmDouble values[3]; 
		for(int i=0;i<3;i++)values[i]=value;
		int lidlist[3];
		this->GetVerticesLidList(&lidlist[0]);
		this->inputs->SetTriaInput(enum_in,P1Enum,3,&lidlist[0],&values[0]);
	}
	else _error_("interpolation type not supported yet");
}
/*}}}*/
void       Tria::SetElementInput(Inputs* inputs,int numindices,int* indices,IssmDouble* values,int enum_in){/*{{{*/

	_assert_(inputs);
	inputs->SetTriaInput(enum_in,P1Enum,numindices,indices,values);

}
/*}}}*/
Element*   Tria::SpawnBasalElement(bool depthaverage_materials){/*{{{*/

	int index1,index2;
	int domaintype;

	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
		case Domain3DsurfaceEnum:
			return this;
		case Domain2DverticalEnum:
			_assert_(HasEdgeOnBase());
			this->EdgeOnBaseIndices(&index1,&index2);
			return SpawnSeg(index1,index2);
		default:
			_error_("not implemented yet");
	}
}
/*}}}*/
Seg*       Tria::SpawnSeg(int index1,int index2){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	this->parameters->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Create Seg*/
	Seg* seg=new Seg();
	seg->id=this->id;
	seg->sid=this->sid;
	seg->lid=this->lid;
	seg->inputs=this->inputs;
	seg->parameters=this->parameters;
	seg->element_type=P1Enum; //Only P1 CG for now (TO BE CHANGED)
	this->SpawnSegHook(xDynamicCast<ElementHook*>(seg),index1,index2);

	seg->iscollapsed = 1;
	seg->collapsed_ids[0] = index1;
	seg->collapsed_ids[1] = index2;

	/*Spawn material*/
	seg->material=(Material*)this->material->copy2(seg);

	/*recover nodes, material*/
	seg->nodes    = (Node**)seg->hnodes[analysis_counter]->deliverp();
	seg->vertices = (Vertex**)seg->hvertices->deliverp();

	/*Return new Seg*/
	return seg;
}
/*}}}*/
Element*   Tria::SpawnTopElement(void){/*{{{*/

	int index1,index2;
	int domaintype;

	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	switch(domaintype){
		case Domain2DhorizontalEnum:
			return this;
		case Domain2DverticalEnum:
			_assert_(HasEdgeOnSurface());
			this->EdgeOnSurfaceIndices(&index1,&index2);
			return SpawnSeg(index2,index1); //reverse order
		default:
			_error_("not implemented yet");
	}
}
/*}}}*/
bool       Tria::IsSpawnedElement(void){/*{{{*/

	if(this->iscollapsed!=0){
		return true;
	}

	return false;

}/*}}}*/
void       Tria::StrainRateparallel(){/*{{{*/

	IssmDouble  epsilon[3];
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainparallel[NUMVERTICES];

	/* Get node coordinates and dof list: */
   IssmDouble  xyz_list[NUMVERTICES][3];
   ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will need*/
	Input *vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input *vy_input = this->GetInput(VyEnum); _assert_(vy_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[2];

		/*strainparallel= Strain rate along the ice flow direction */
		strainparallel[iv]=(vx*vx*(strainxx)+vy*vy*(strainyy)+2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->AddInput(StrainRateparallelEnum,&strainparallel[0],P1DGEnum);
}
/*}}}*/
void       Tria::StrainRateperpendicular(){/*{{{*/

	IssmDouble  epsilon[3];
	IssmDouble  vx,vy,vel;
	IssmDouble  strainxx;
	IssmDouble  strainxy;
	IssmDouble  strainyy;
	IssmDouble  strainperpendicular[NUMVERTICES];

	/* Get node coordinates and dof list: */
   IssmDouble  xyz_list[NUMVERTICES][3];
   ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Retrieve all inputs we will need*/
	Input *vx_input = this->GetInput(VxEnum); _assert_(vx_input);
	Input *vy_input = this->GetInput(VyEnum); _assert_(vy_input);

	/* Start looping on the number of vertices: */
	GaussTria gauss;
	for (int iv=0;iv<NUMVERTICES;iv++){
		gauss.GaussVertex(iv);

		/* Get the value we need*/
		vx_input->GetInputValue(&vx,&gauss);
		vy_input->GetInputValue(&vy,&gauss);
		vel=vx*vx+vy*vy;

		/*Compute strain rate viscosity and pressure: */
		this->StrainRateSSA(&epsilon[0],&xyz_list[0][0],&gauss,vx_input,vy_input);
		strainxx=epsilon[0];
		strainyy=epsilon[1];
		strainxy=epsilon[2];

		/*strainperpendicular= Strain rate perpendicular to the ice flow direction */
		strainperpendicular[iv]=(vx*vx*(strainyy)+vy*vy*(strainxx)-2*vy*vx*strainxy)/(vel+1.e-14);
	}

	/*Add input*/
	this->AddInput(StrainRateperpendicularEnum,&strainperpendicular[0],P1DGEnum);
}
/*}}}*/
IssmDouble Tria::SurfaceArea(void){/*{{{*/

	IssmDouble S;
	IssmDouble normal[3];
	IssmDouble v13[3],v23[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*If on water, return 0: */
	if(!IsIceInElement()) return 0.;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0][i]-xyz_list[2][i];
		v23[i]=xyz_list[1][i]-xyz_list[2][i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];

	S = 0.5 * sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

	/*Return: */
	return S;
}
/*}}}*/
int        Tria::TensorInterpolation(void){/*{{{*/
	return TriaRef::TensorInterpolation(this->element_type);
}
/*}}}*/
IssmDouble Tria::TimeAdapt(void){/*{{{*/

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
	IssmDouble maxabsvx = vx_input->GetInputMaxAbs();
	IssmDouble maxabsvy = vy_input->GetInputMaxAbs();

	/* Get node coordinates and dof list: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	IssmDouble minx=xyz_list[0][0];
	IssmDouble maxx=xyz_list[0][0];
	IssmDouble miny=xyz_list[0][1];
	IssmDouble maxy=xyz_list[0][1];

	for(int i=1;i<NUMVERTICES;i++){
		if(xyz_list[i][0]<minx) minx=xyz_list[i][0];
		if(xyz_list[i][0]>maxx) maxx=xyz_list[i][0];
		if(xyz_list[i][1]<miny) miny=xyz_list[i][1];
		if(xyz_list[i][1]>maxy) maxy=xyz_list[i][1];
	}
	IssmDouble dx=maxx-minx;
	IssmDouble dy=maxy-miny;

	/*CFL criterion: */
	IssmDouble dt = C/(maxabsvx/dx+maxabsvy/dy + 1.e-18);

	/*Check hydro timestep also and take the minimum*/
	if(ishydro){
		Input* vx_input = this->GetInput(HydrologyWaterVxEnum);
		Input* vy_input = this->GetInput(HydrologyWaterVyEnum);
		IssmDouble maxabsvx = 0.03;
		IssmDouble maxabsvy = 0.03;
		if(vx_input) maxabsvx = vx_input->GetInputMaxAbs();
		if(vy_input) maxabsvy = vy_input->GetInputMaxAbs();

		/*CFL criterion: */
		IssmDouble dt2 = C/(maxabsvx/dx+maxabsvy/dy + 1.e-18);
		if(dt2<dt) dt = dt2;
	}

	return dt;
}
/*}}}*/
IssmDouble Tria::TotalCalvingFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               domaintype,index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

   IssmDouble  xyz_list[NUMVERTICES][3];
   ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		_error_("not implemented");
	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble calvingratex,calvingratey,thickness,Jdet;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input    = this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input* calvingratex_input = this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
	Input* calvingratey_input = this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		calvingratex_input->GetInputValue(&calvingratex,gauss);
		calvingratey_input->GetInputValue(&calvingratey,gauss);
		this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*(calvingratex*normal[0] + calvingratey*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Tria::TotalCalvingMeltingFluxLevelset(bool scaled){/*{{{*/

	/*Make sure there is an ice front here*/
	if(!IsIceInElement() || !IsZeroLevelset(MaskIceLevelsetEnum)) return 0;

	/*Scaled not implemented yet...*/
	_assert_(!scaled);

	int               domaintype,index1,index2;
	const IssmPDouble epsilon = 1.e-15;
	IssmDouble        s1,s2;
	IssmDouble        gl[NUMVERTICES];
	IssmDouble        xyz_front[2][3];

   IssmDouble  xyz_list[NUMVERTICES][3];
   ::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*Recover parameters and values*/
	parameters->FindParam(&domaintype,DomainTypeEnum);
	Element::GetInputListOnVertices(&gl[0],MaskIceLevelsetEnum);

	/*Be sure that values are not zero*/
	if(gl[0]==0.) gl[0]=gl[0]+epsilon;
	if(gl[1]==0.) gl[1]=gl[1]+epsilon;
	if(gl[2]==0.) gl[2]=gl[2]+epsilon;

	if(domaintype==Domain2DverticalEnum){
		_error_("not implemented");
	}
	else if(domaintype==Domain2DhorizontalEnum || domaintype==Domain3DEnum || domaintype==Domain3DsurfaceEnum){
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

	}
	else _error_("mesh type "<<EnumToStringx(domaintype)<<"not supported yet ");

	/*Some checks in debugging mode*/
	_assert_(s1>=0 && s1<=1.);
	_assert_(s2>=0 && s2<=1.);

	/*Get normal vector*/
	IssmDouble normal[3];
	this->NormalSection(&normal[0],&xyz_front[0][0]);
	normal[0] = -normal[0];
	normal[1] = -normal[1];

	/*Get inputs*/
	IssmDouble flux = 0.;
	IssmDouble calvingratex,calvingratey,vx,vy,vel,meltingrate,meltingratex,meltingratey,thickness,Jdet,groundedice;
	IssmDouble rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* thickness_input    = this->GetInput(ThicknessEnum); _assert_(thickness_input);
	Input *gr_input           = this->GetInput(MaskOceanLevelsetEnum); _assert_(gr_input);
	Input* calvingratex_input = this->GetInput(CalvingratexEnum); _assert_(calvingratex_input);
	Input* calvingratey_input = this->GetInput(CalvingrateyEnum); _assert_(calvingratey_input);
	Input* vx_input           = this->GetInput(VxEnum); _assert_(vx_input);
	Input* vy_input           = this->GetInput(VyEnum); _assert_(vy_input);
	Input* meltingrate_input  = this->GetInput(CalvingMeltingrateEnum); _assert_(meltingrate_input);

	/*Start looping on Gaussian points*/
	Gauss* gauss=this->NewGauss(&xyz_list[0][0],&xyz_front[0][0],3);
	while(gauss->next()){
		thickness_input->GetInputValue(&thickness,gauss);
		calvingratex_input->GetInputValue(&calvingratex,gauss);
		calvingratey_input->GetInputValue(&calvingratey,gauss);
		gr_input->GetInputValue(&groundedice,gauss);
		vx_input->GetInputValue(&vx,gauss);
		vy_input->GetInputValue(&vy,gauss);
		vel=vx*vx+vy*vy;
		meltingrate_input->GetInputValue(&meltingrate,gauss);
		if(groundedice<0) meltingrate = 0.;
		meltingratex=meltingrate*vx/(sqrt(vel)+1.e-14);
		meltingratey=meltingrate*vy/(sqrt(vel)+1.e-14);
		this->JacobianDeterminantSurface(&Jdet,&xyz_front[0][0],gauss);

		flux += rho_ice*Jdet*gauss->weight*thickness*((calvingratex+meltingratex)*normal[0] + (calvingratey+meltingratey)*normal[1]);
	}

	/*Clean up and return*/
	delete gauss;
	return flux;
}
/*}}}*/
IssmDouble Tria::TotalFloatingBmb(bool scaled){/*{{{*/

	/*The fbmb[kg yr-1] of one element is area[m2] * melting_rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble fbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,floatingmelt,Jdet,scalefactor;
	IssmDouble Total_Fbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement())return 0;

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* floatingmelt_input = this->GetInput(BasalforcingsFloatingiceMeltingRateEnum); _assert_(floatingmelt_input);
	Input* gllevelset_input   = this->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
	Input* scalefactor_input  = NULL;
	if(scaled==true){
		scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
	}
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,1-mainlyfloating,3);
	while(gauss->next()){
		this->JacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
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
IssmDouble Tria::TotalGroundedBmb(bool scaled){/*{{{*/

	/*The gbmb[kg yr-1] of one element is area[m2] * gounded melting rate [kg m^-2 yr^-1]*/
	int        point1;
	bool       mainlyfloating;
	IssmDouble gbmb=0;
	IssmDouble rho_ice,fraction1,fraction2,groundedmelt,Jdet,scalefactor;
	IssmDouble Total_Gbmb=0;
	IssmDouble xyz_list[NUMVERTICES][3];
	Gauss*     gauss     = NULL;

   if(!IsIceInElement())return 0;

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);
	Input* groundedmelt_input = this->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedmelt_input);
	Input* gllevelset_input = this->GetInput(MaskOceanLevelsetEnum); _assert_(gllevelset_input);
	Input* scalefactor_input = NULL;
	if(scaled==true){
		scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
	}
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	this->GetGroundedPart(&point1,&fraction1,&fraction2,&mainlyfloating,MaskOceanLevelsetEnum,0);
	/* Start  looping on the number of gaussian points: */
	gauss = this->NewGauss(point1,fraction1,fraction2,mainlyfloating,2);
	while(gauss->next()){
		this->JacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
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
IssmDouble Tria::TotalSmb(bool scaled){/*{{{*/

	/*The smb[kg yr-1] of one element is area[m2] * smb [kg m^-2 yr^-1]*/
	IssmDouble base,smb,rho_ice,scalefactor;
	IssmDouble Total_Smb=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

   if(!IsIceInElement())return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));	// area of element in m2

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* weights       = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble* smb_vertices  = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble f1,f2,phi;

		Element::GetInputListOnVertices(&smb_vertices[0],SmbMassBalanceEnum);
		GetFractionGeometry(weights,&phi,&point,&f1,&f2,&mainlyice,lsf);
		smb = 0.0;
		for(int i=0;i<NUMVERTICES;i++) smb += weights[i]*smb_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
      }
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(weights);
      xDelete<IssmDouble>(smb_vertices);
	}
	else{
		/*Fully ice-covered element*/
		Input* smb_input = this->GetInput(SmbMassBalanceEnum); _assert_(smb_input);
		smb_input->GetInputAverage(&smb);   // average smb on element in m ice s-1

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

   Total_Smb=rho_ice*base*smb*scalefactor;	// smb on element in kg s-1

	/*Return: */
	return Total_Smb;
}
/*}}}*/
IssmDouble Tria::TotalSmbMelt(bool scaled){/*{{{*/

	/*The smbmelt[kg yr-1] of one element is area[m2] * smbmelt [kg m^-2 yr^-1]*/
	IssmDouble base,smbmelt,rho_ice,scalefactor;
	IssmDouble Total_Melt=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

   if(!IsIceInElement())return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));	// area of element in m2

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* weights       = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble* smbmelt_vertices  = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble f1,f2,phi;

		Element::GetInputListOnVertices(&smbmelt_vertices[0],SmbMeltEnum);
		GetFractionGeometry(weights,&phi,&point,&f1,&f2,&mainlyice,lsf);
		smbmelt = 0.0;
		for(int i=0;i<NUMVERTICES;i++) smbmelt += weights[i]*smbmelt_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
      }
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(weights);
      xDelete<IssmDouble>(smbmelt_vertices);
	}
	else{
		/*Fully ice-covered element*/
		Input* smbmelt_input = this->GetInput(SmbMeltEnum); _assert_(smbmelt_input);
		smbmelt_input->GetInputAverage(&smbmelt);   // average smbmelt on element in m ice s-1

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

   Total_Melt=rho_ice*base*smbmelt*scalefactor;	// smbmelt on element in kg s-1

	/*Return: */
	return Total_Melt;
}
/*}}}*/
IssmDouble Tria::TotalSmbRefreeze(bool scaled){/*{{{*/

	/*The smb[kg yr-1] of one element is area[m2] * smb [kg m^-2 yr^-1]*/
	IssmDouble base,smbrefreeze,rho_ice,scalefactor;
	IssmDouble Total_Refreeze=0;
	IssmDouble lsf[NUMVERTICES];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Get material parameters :*/
	rho_ice=FindParam(MaterialsRhoIceEnum);

   if(!IsIceInElement())return 0;

	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*First calculate the area of the base (cross section triangle)
	 * http://en.wikipedia.org/wiki/Triangle
	 * base = 1/2 abs((xA-xC)(yB-yA)-(xA-xB)(yC-yA))*/
	base = 1./2. * fabs((xyz_list[0][0]-xyz_list[2][0])*(xyz_list[1][1]-xyz_list[0][1]) - (xyz_list[0][0]-xyz_list[1][0])*(xyz_list[2][1]-xyz_list[0][1]));	// area of element in m2

	/*Now get the average SMB over the element*/
	Element::GetInputListOnVertices(&lsf[0],MaskIceLevelsetEnum);
	if(lsf[0]*lsf[1]<=0 || lsf[0]*lsf[2]<=0 || lsf[1]*lsf[2]<=0){
		/*Partially ice-covered element*/
		bool mainlyice;
      int point;
      IssmDouble* weights       = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble* smbrefreeze_vertices  = xNew<IssmDouble>(NUMVERTICES);
      IssmDouble f1,f2,phi;

		Element::GetInputListOnVertices(&smbrefreeze_vertices[0],SmbRefreezeEnum);
		GetFractionGeometry(weights,&phi,&point,&f1,&f2,&mainlyice,lsf);
		smbrefreeze = 0.0;
		for(int i=0;i<NUMVERTICES;i++) smbrefreeze += weights[i]*smbrefreeze_vertices[i];

		if(scaled==true){
         IssmDouble* scalefactor_vertices = xNew<IssmDouble>(NUMVERTICES);
         Element::GetInputListOnVertices(&scalefactor_vertices[0],MeshScaleFactorEnum);
         scalefactor = 0.0;
         for(int i=0;i<NUMVERTICES;i++) scalefactor += weights[i]/phi*scalefactor_vertices[i];
         xDelete<IssmDouble>(scalefactor_vertices);
      }
		else scalefactor = 1.0;

		/*Cleanup*/
      xDelete<IssmDouble>(weights);
      xDelete<IssmDouble>(smbrefreeze_vertices);
	}
	else{
		/*Fully ice-covered element*/
		Input* smbrefreeze_input = this->GetInput(SmbRefreezeEnum); _assert_(smbrefreeze_input);
		smbrefreeze_input->GetInputAverage(&smbrefreeze);   // average smbrefreeze on element in m ice s-1

		if(scaled==true){
			Input* scalefactor_input = this->GetInput(MeshScaleFactorEnum); _assert_(scalefactor_input);
			scalefactor_input->GetInputAverage(&scalefactor);// average scalefactor on element
		}
		else scalefactor=1.0;
	}

   Total_Refreeze=rho_ice*base*smbrefreeze*scalefactor;	// smbrefreeze on element in kg s-1

	/*Return: */
	return Total_Refreeze;
}
/*}}}*/
void       Tria::Update(Inputs* inputs,int index, IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){/*{{{*/

	/*Intermediaries*/
	int  numnodes;
	int* tria_node_ids = NULL;

	/*Checks if debuging*/
	_assert_(iomodel->elements);
	_assert_(index==this->sid);

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P0DGEnum:
			numnodes        = 1;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]= index + 1;
			break;
		case P1Enum:
			numnodes        = 3;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			break;
		case P1DGEnum:
			numnodes        = 3;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=3*index+1;
			tria_node_ids[1]=3*index+2;
			tria_node_ids[2]=3*index+3;
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes        = 4;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+index+1;
			break;
		case P2Enum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			break;
		case P2bubbleEnum: case P2bubblecondensedEnum:
			numnodes        = 7;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+index+1;
			break;
		case P1P1Enum: case P1P1GLSEnum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];

			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elements[3*index+0];
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elements[3*index+1];
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elements[3*index+2];
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes       = 7;
			tria_node_ids  = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+index+1;

			tria_node_ids[4]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+0];
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+1];
			tria_node_ids[6]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[3*index+2];
			break;
		case TaylorHoodEnum:
		case XTaylorHoodEnum:
			numnodes        = 9;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;

			tria_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+0];
			tria_node_ids[7]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+1];
			tria_node_ids[8]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[3*index+2];
			break;
		case LATaylorHoodEnum:
			numnodes        = 6;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			break;
		case CrouzeixRaviartEnum:
			numnodes        = 10;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+index+1;

			tria_node_ids[7]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+1;
			tria_node_ids[8]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+2;
			tria_node_ids[9]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->numberofelements+3*index+3;
			break;
		case LACrouzeixRaviartEnum:
			numnodes        = 7;
			tria_node_ids   = xNew<int>(numnodes);
			tria_node_ids[0]=iomodel->elements[3*index+0];
			tria_node_ids[1]=iomodel->elements[3*index+1];
			tria_node_ids[2]=iomodel->elements[3*index+2];
			tria_node_ids[3]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+0]+1;
			tria_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+1]+1;
			tria_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[3*index+2]+1;
			tria_node_ids[6]=iomodel->numberofvertices+iomodel->numberofedges+index+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(tria_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(tria_node_ids);
}
/*}}}*/
void       Tria::UpdateConstraintsExtrudeFromBase(void){/*{{{*/

	if(!HasNodeOnBase()) return;

	int        extrusioninput;
	IssmDouble value,isonbase;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = this->GetInput(extrusioninput);      _assert_(input);
	Input* onbase = this->GetInput(MeshVertexonbaseEnum); _assert_(onbase);

	GaussTria gauss;
	for(int iv=0;iv<this->NumberofNodes(this->element_type);iv++){
		gauss.GaussNode(this->element_type,iv);
		onbase->GetInputValue(&isonbase,&gauss);
		if(isonbase==1.){
			input->GetInputValue(&value,&gauss);
			this->nodes[iv]->ApplyConstraint(0,value);
		}
	}

}
/*}}}*/
void       Tria::UpdateConstraintsExtrudeFromTop(void){/*{{{*/

	if(!HasNodeOnSurface()) return;

	int        extrusioninput;
	IssmDouble value,isonsurface;

	this->parameters->FindParam(&extrusioninput,InputToExtrudeEnum);
	Input* input = this->GetInput(extrusioninput); _assert_(input);
	Input* onsurf = this->GetInput(MeshVertexonsurfaceEnum); _assert_(onsurf);

	GaussTria gauss;
	for(int iv=0;iv<this->NumberofNodes(this->element_type);iv++){
		gauss.GaussNode(this->element_type,iv);
		onsurf->GetInputValue(&isonsurface,&gauss);
		if(isonsurface==1.){
			input->GetInputValue(&value,&gauss);
			this->nodes[iv]->ApplyConstraint(0,value);
		}
	}
}
/*}}}*/
int        Tria::UpdatePotentialUngrounding(IssmDouble* vertices_potentially_ungrounding,Vector<IssmDouble>* vec_nodes_on_iceshelf,IssmDouble* nodes_on_iceshelf){/*{{{*/

	int i;
	int nflipped=0;

	/*Go through nodes, and whoever is on the potential_ungrounding, ends up in nodes_on_iceshelf: */
	for(i=0;i<3;i++){
		if (reCast<bool>(vertices_potentially_ungrounding[vertices[i]->Pid()])){
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
void       Tria::ValueP1DerivativesOnGauss(IssmDouble* dvalue,IssmDouble* values,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/
	TriaRef::GetInputDerivativeValue(dvalue,values,xyz_list,gauss,P1Enum);
}
/*}}}*/
void       Tria::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	TriaRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int        Tria::VelocityInterpolation(void){/*{{{*/
	return TriaRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
int        Tria::VertexConnectivity(int vertexindex){/*{{{*/
	_assert_(this->vertices);
	return this->vertices[vertexindex]->Connectivity();
}
/*}}}*/
void       Tria::WriteFieldIsovalueSegment(DataSet* segments,int fieldenum,IssmDouble fieldvalue){/*{{{*/

	_assert_(fieldvalue==0.); //field value != 0 not implemented yet

	/*Get field on vertices (we do not allow for higher order elements!!)*/
	IssmDouble lsf[NUMVERTICES];
	Element::GetInputListOnVertices(&lsf[0],fieldenum);

	/*1. check that we do cross fieldvalue in this element*/
	IssmDouble minvalue = lsf[0];
	IssmDouble maxvalue = lsf[0];
	for(int i=1;i<NUMVERTICES;i++){
		if(lsf[i]>maxvalue) maxvalue = lsf[i];
		if(lsf[i]<minvalue) minvalue = lsf[i];
	}
	if(minvalue>fieldvalue) return;
	if(maxvalue<fieldvalue) return;

	/*2. Find coordinates of where levelset crosses 0*/
	int         numiceverts;
	IssmDouble  s[2],x[2],y[2];
	int        *indices = NULL;
	this->GetLevelsetIntersection(&indices, &numiceverts,&s[0],fieldenum,fieldvalue);
	_assert_(numiceverts);

	/*3 Write coordinates*/
	IssmDouble  xyz_list[NUMVERTICES][3];
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);
	int counter = 0;
	if((numiceverts>0) && (numiceverts<NUMVERTICES)){
		for(int i=0;i<numiceverts;i++){
			for(int n=numiceverts;n<NUMVERTICES;n++){ // iterate over no-ice vertices
				x[counter] = xyz_list[indices[i]][0]+s[counter]*(xyz_list[indices[n]][0]-xyz_list[indices[i]][0]);
				y[counter] = xyz_list[indices[i]][1]+s[counter]*(xyz_list[indices[n]][1]-xyz_list[indices[i]][1]);
				counter++;
			}
		}
	}
	else if(numiceverts==NUMVERTICES){ //NUMVERTICES ice vertices: calving front lies on element edge

		for(int i=0;i<NUMVERTICES;i++){
			if(lsf[indices[i]]==0.){
				x[counter]=xyz_list[indices[i]][0];
				y[counter]=xyz_list[indices[i]][1];
				counter++;
			}
			if(counter==2) break;
		}
		if(counter==1){
			/*We actually have only 1 vertex on levelset, write a single point as a segment*/
			x[counter]=x[0];
			y[counter]=y[0];
			counter++;
		}
	}
	else{
		_error_("not sure what's going on here...");
	}

	/*4. Write segment*/
	_assert_(counter==2);
	segments->AddObject(new Contour<IssmDouble>(segments->Size()+1,2,&x[0],&y[0],false));

	/*Cleanup and return*/
	xDelete<int>(indices);
}
/*}}}*/

#ifdef _HAVE_ESA_
void    Tria::EsaGeodetic2D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,Vector<IssmDouble>* pX,Vector<IssmDouble>* pY,IssmDouble* xx,IssmDouble* yy){ /*{{{*/

	/*diverse:*/
	int gsize;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble area;
	IssmDouble earth_radius = 6371012.0;	// Earth's radius [m]
	IssmDouble g_earth = 9.81;	// Gravitational acceleration on Earth's surface [m/s2]
	IssmDouble I;		//ice/water loading
	IssmDouble rho_ice, rho_earth;

	/*precomputed elastic green functions:*/
	IssmDouble* U_elastic_precomputed = NULL;
	IssmDouble* H_elastic_precomputed = NULL;
	IssmDouble* G_elastic_precomputed = NULL;
	int         M, hemi;

	/*computation of Green functions:*/
	IssmDouble* U_elastic= NULL;
	IssmDouble* N_elastic= NULL;
	IssmDouble* E_elastic= NULL;
	IssmDouble* X_elastic= NULL;
	IssmDouble* Y_elastic= NULL;
	IssmDouble* G_elastic= NULL;

	/*optimization:*/
	bool store_green_functions=false;

	/*Compute ice thickness change: */
	Input* deltathickness_input=this->GetInput(DeltaIceThicknessEnum);
	if (!deltathickness_input)_error_("delta thickness input needed to compute elastic adjustment!");
	deltathickness_input->GetInputAverage(&I);

	/*early return if we are not on the (ice) loading point: */
	if(I==0) return;

	/*recover material parameters: */
	rho_ice=FindParam(MaterialsRhoIceEnum);
	rho_earth=FindParam(MaterialsEarthDensityEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*which hemisphere? for north-south, east-west components*/
	this->parameters->FindParam(&hemi,EsaHemisphereEnum);

	/*compute area of element:*/
	area=GetArea();

	/*figure out gravity center of our element (Cartesian): */
	IssmDouble x_element, y_element;
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x_element=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y_element=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;

	/*recover elastic Green's functions for displacement:*/
	DoubleVecParam* U_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaUElasticEnum)); _assert_(U_parameter);
	DoubleVecParam* H_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaHElasticEnum)); _assert_(H_parameter);
	DoubleVecParam* G_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaGElasticEnum)); _assert_(G_parameter);
	U_parameter->GetParameterValueByPointer(&U_elastic_precomputed,&M);
	H_parameter->GetParameterValueByPointer(&H_elastic_precomputed,&M);
	G_parameter->GetParameterValueByPointer(&G_elastic_precomputed,&M);

	/*initialize: */
	U_elastic=xNewZeroInit<IssmDouble>(gsize);
	N_elastic=xNewZeroInit<IssmDouble>(gsize);
	E_elastic=xNewZeroInit<IssmDouble>(gsize);
	X_elastic=xNewZeroInit<IssmDouble>(gsize);
	Y_elastic=xNewZeroInit<IssmDouble>(gsize);
	G_elastic=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* U_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* N_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* E_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* G_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* X_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* Y_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble dx, dy, dist, alpha, ang, ang2;
	IssmDouble N_azim, E_azim, X_azim, Y_azim;

	for(int i=0;i<gsize;i++){

		indices[i]=i;

		IssmDouble N_azim=0;
		IssmDouble E_azim=0;

		/*Compute alpha angle between centroid and current vertex: */
		dx = x_element - xx[i];		dy = y_element - yy[i];
		dist = sqrt(pow(dx,2)+pow(dy,2));						// distance between vertex and elemental centroid [m]
		alpha = dist*360.0/(2*M_PI*earth_radius) * M_PI/180.0;	// [in radians] 360 degree = 2*pi*earth_radius

		/*Compute azimuths, both north and east components: */
		ang = M_PI/2 - atan2(dy,dx);		// this is bearing angle!
		Y_azim = -cos(ang);
		X_azim = -sin(ang);

		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		int index=reCast<int,IssmDouble>(alpha/M_PI*(M-1));
		U_elastic[i] += U_elastic_precomputed[index];
		Y_elastic[i] += H_elastic_precomputed[index]*Y_azim;
		X_elastic[i] += H_elastic_precomputed[index]*X_azim;
		G_elastic[i] += G_elastic_precomputed[index];

		/*Add all components to the pUp solution vectors:*/
		U_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*U_elastic[i];
		Y_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*Y_elastic[i];
		X_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*X_elastic[i];
		G_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*G_elastic[i]*g_earth/earth_radius;

		/*North-south, East-west components */
		if (hemi == -1) {
			ang2 = M_PI/2 - atan2(yy[i],xx[i]);
		}
		else if (hemi == 1) {
			ang2 = M_PI/2 - atan2(-yy[i],-xx[i]);
		}
		if (hemi != 0){
			N_azim = Y_azim*cos(ang2) + X_azim*sin(ang2);
			E_azim = X_azim*cos(ang2) - Y_azim*sin(ang2);
			N_elastic[i] += H_elastic_precomputed[index]*N_azim;
			E_elastic[i] += H_elastic_precomputed[index]*E_azim;
			N_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*N_elastic[i];
			E_values[i]+=3*rho_ice/rho_earth*area/(4*M_PI*pow(earth_radius,2))*I*E_elastic[i];
		}
	}

	pUp->SetValues(gsize,indices,U_values,ADD_VAL);
	pNorth->SetValues(gsize,indices,N_values,ADD_VAL);
	pEast->SetValues(gsize,indices,E_values,ADD_VAL);
	pGravity->SetValues(gsize,indices,G_values,ADD_VAL);
	pX->SetValues(gsize,indices,X_values,ADD_VAL);
	pY->SetValues(gsize,indices,Y_values,ADD_VAL);

	/*Free resources:*/
	xDelete<int>(indices);
	xDelete<IssmDouble>(U_values); xDelete<IssmDouble>(N_values); xDelete<IssmDouble>(E_values); xDelete<IssmDouble>(G_values);
	xDelete<IssmDouble>(U_elastic); xDelete<IssmDouble>(N_elastic); xDelete<IssmDouble>(E_elastic); xDelete<IssmDouble>(G_elastic);
	xDelete<IssmDouble>(X_values); xDelete<IssmDouble>(Y_values);
	xDelete<IssmDouble>(X_elastic); xDelete<IssmDouble>(Y_elastic);

	return;
}
/*}}}*/
void    Tria::EsaGeodetic3D(Vector<IssmDouble>* pUp,Vector<IssmDouble>* pNorth,Vector<IssmDouble>* pEast,Vector<IssmDouble>* pGravity,IssmDouble* latitude,IssmDouble* longitude,IssmDouble* radius,IssmDouble* xx,IssmDouble* yy,IssmDouble* zz){ /*{{{*/

	/*diverse:*/
	int gsize;
	bool spherical=true;
	IssmDouble llr_list[NUMVERTICES][3];
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble area,planetarea;
	IssmDouble earth_radius = 6371012.0;	// Earth's radius [m]
	IssmDouble g_earth = 9.81;	// Gravitational acceleration on Earth's surface [m/s2]
	IssmDouble I;		//ice/water loading
	IssmDouble late,longe,re;
	IssmDouble lati,longi,ri;
	IssmDouble rho_ice,rho_earth;
	IssmDouble minlong=400;
	IssmDouble maxlong=-20;

	/*precomputed elastic green functions:*/
	IssmDouble* U_elastic_precomputed = NULL;
	IssmDouble* H_elastic_precomputed = NULL;
	IssmDouble* G_elastic_precomputed = NULL;
	int         M;

	/*computation of Green functions:*/
	IssmDouble* U_elastic= NULL;
	IssmDouble* N_elastic= NULL;
	IssmDouble* E_elastic= NULL;
	IssmDouble* G_elastic= NULL;

	/*optimization:*/
	bool store_green_functions=false;

	/*Compute ice thickness change: */
	Input* deltathickness_input=this->GetInput(DeltaIceThicknessEnum);
	if (!deltathickness_input)_error_("delta thickness input needed to compute elastic adjustment!");
	deltathickness_input->GetInputAverage(&I);

	/*early return if we are not on the (ice) loading point: */
	if(I==0) return;

	/*recover material parameters: */
	rho_ice=FindParam(MaterialsRhoIceEnum);
	rho_earth=FindParam(MaterialsEarthDensityEnum);

	/*recover earth area and radius: */
	this->parameters->FindParam(&planetarea,SolidearthPlanetAreaEnum);

	/*how many dofs are we working with here? */
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);

	/*Get area of element: precomputed in the sealevelchange_geometry:*/
	area=GetAreaSpherical();

	/*element centroid (spherical): */
	/* Where is the centroid of this element?:{{{*/
	::GetVerticesCoordinates(&llr_list[0][0],this->vertices,NUMVERTICES,spherical);

	minlong=400; maxlong=-20;
	for (int i=0;i<NUMVERTICES;i++){
		llr_list[i][0]=(90-llr_list[i][0]);
		if(llr_list[i][1]<0)llr_list[i][1]=180+(180+llr_list[i][1]);
		if(llr_list[i][1]>maxlong)maxlong=llr_list[i][1];
		if(llr_list[i][1]<minlong)minlong=llr_list[i][1];
	}
	if(minlong==0 && maxlong>180){
		if (llr_list[0][1]==0)llr_list[0][1]=360;
		if (llr_list[1][1]==0)llr_list[1][1]=360;
		if (llr_list[2][1]==0)llr_list[2][1]=360;
	}

	// correction at the north pole: given longitude of the North pole a definition
	// closer to the other two vertices.
	if(llr_list[0][0]==0)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==0)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==0)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	// correction at the north pole: given longitude of the North pole a definition
	// closer to the other two vertices.
	if(llr_list[0][0]==180)llr_list[0][1]=(llr_list[1][1]+llr_list[2][1])/2.0;
	if(llr_list[1][0]==180)llr_list[1][1]=(llr_list[0][1]+llr_list[2][1])/2.0;
	if(llr_list[2][0]==180)llr_list[2][1]=(llr_list[0][1]+llr_list[1][1])/2.0;

	late=(llr_list[0][0]+llr_list[1][0]+llr_list[2][0])/3.0;
	longe=(llr_list[0][1]+llr_list[1][1]+llr_list[2][1])/3.0;

	late=90-late;
	if(longe>180)longe=longe-360;

	late=late/180.*M_PI;
	longe=longe/180.*M_PI;
	/*}}}*/

	/*figure out gravity center of our element (Cartesian): */
	IssmDouble x_element, y_element, z_element;
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	x_element=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	y_element=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;
	z_element=(xyz_list[0][2]+xyz_list[1][2]+xyz_list[2][2])/3.0;

	/*recover elastic Green's functions for displacement:*/
	DoubleVecParam* U_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaUElasticEnum)); _assert_(U_parameter);
	DoubleVecParam* H_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaHElasticEnum)); _assert_(H_parameter);
	DoubleVecParam* G_parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(EsaGElasticEnum)); _assert_(G_parameter);
	U_parameter->GetParameterValueByPointer(&U_elastic_precomputed,&M);
	H_parameter->GetParameterValueByPointer(&H_elastic_precomputed,&M);
	G_parameter->GetParameterValueByPointer(&G_elastic_precomputed,&M);

	/*initialize: */
	U_elastic=xNewZeroInit<IssmDouble>(gsize);
	N_elastic=xNewZeroInit<IssmDouble>(gsize);
	E_elastic=xNewZeroInit<IssmDouble>(gsize);
	G_elastic=xNewZeroInit<IssmDouble>(gsize);

	int* indices=xNew<int>(gsize);
	IssmDouble* U_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* N_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* E_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble* G_values=xNewZeroInit<IssmDouble>(gsize);
	IssmDouble alpha;
	IssmDouble delPhi,delLambda;
	IssmDouble dx, dy, dz, x, y, z;
	IssmDouble N_azim, E_azim;

	for(int i=0;i<gsize;i++){

		indices[i]=i;

		/*Compute alpha angle between centroid and current vertex: */
		lati=latitude[i]/180.*M_PI; longi=longitude[i]/180.*M_PI;

		delPhi=fabs(lati-late); delLambda=fabs(longi-longe);
		alpha=2.*asin(sqrt(pow(sin(delPhi/2),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));

		/*Compute azimuths, both north and east components: */
		x = xx[i]; y = yy[i]; z = zz[i];
		if(latitude[i]==90){
			x=1e-12; y=1e-12;
		}
		if(latitude[i]==-90){
			x=1e-12; y=1e-12;
		}
		dx = x_element-x; dy = y_element-y; dz = z_element-z;
		N_azim = -(-z*x*dx-z*y*dy+(pow(x,2)+pow(y,2))*dz) /pow((pow(x,2)+pow(y,2))*(pow(x,2)+pow(y,2)+pow(z,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);
		E_azim = -(-y*dx+x*dy) /pow((pow(x,2)+pow(y,2))*(pow(dx,2)+pow(dy,2)+pow(dz,2)),0.5);

		/*Elastic component  (from Eq 17 in Adhikari et al, GMD 2015): */
		int index=reCast<int,IssmDouble>(alpha/M_PI*(M-1));
		U_elastic[i] += U_elastic_precomputed[index];
		N_elastic[i] += H_elastic_precomputed[index]*N_azim;
		E_elastic[i] += H_elastic_precomputed[index]*E_azim;
		G_elastic[i] += G_elastic_precomputed[index];

		/*Add all components to the pUp solution vectors:*/
		U_values[i]+=3*rho_ice/rho_earth*area/planetarea*I*U_elastic[i];
		N_values[i]+=3*rho_ice/rho_earth*area/planetarea*I*N_elastic[i];
		E_values[i]+=3*rho_ice/rho_earth*area/planetarea*I*E_elastic[i];
		G_values[i]+=3*rho_ice/rho_earth*area/planetarea*I*G_elastic[i]*g_earth/earth_radius;
	}
	pUp->SetValues(gsize,indices,U_values,ADD_VAL);
	pNorth->SetValues(gsize,indices,N_values,ADD_VAL);
	pEast->SetValues(gsize,indices,E_values,ADD_VAL);
	pGravity->SetValues(gsize,indices,G_values,ADD_VAL);

	/*Free resources:*/
	xDelete<int>(indices);
	xDelete<IssmDouble>(U_values); xDelete<IssmDouble>(N_values); xDelete<IssmDouble>(E_values); xDelete<IssmDouble>(G_values);
	xDelete<IssmDouble>(U_elastic); xDelete<IssmDouble>(N_elastic); xDelete<IssmDouble>(E_elastic); xDelete<IssmDouble>(G_elastic);

	return;
}
/*}}}*/
#endif
#ifdef _HAVE_SEALEVELCHANGE_
void       Tria::GiaDeflection(Vector<IssmDouble>* wg,Vector<IssmDouble>* dwgdt, Matlitho* litho, IssmDouble* x, IssmDouble* y){/*{{{*/

	IssmDouble xyz_list[NUMVERTICES][3];

	/*gia solution parameters:*/
	IssmDouble ice_mask;

	/*output: */
	IssmDouble  wi;
	IssmDouble  dwidt;

	/*arguments to GiaDeflectionCorex: */
	GiaDeflectionCoreArgs arguments;

	/*how many dofs are we working with here? */
	int gsize;
	IssmDouble yts;
	this->parameters->FindParam(&gsize,MeshNumberofverticesEnum);
	this->parameters->FindParam(&yts,ConstantsYtsEnum);

	/*recover gia solution parameters: */
	int cross_section_shape;
	this->parameters->FindParam(&cross_section_shape,SolidearthSettingsCrossSectionShapeEnum);

	/*what time is it? :*/
	IssmDouble currenttime;
	this->parameters->FindParam(&currenttime,TimeEnum);

	/*recover material parameters: */
	IssmDouble rho_ice                   = FindParam(MaterialsRhoIceEnum);

	/*recover mantle and lithosphere material properties:*/
	int numlayers=litho->numlayers;

	/*lithosphere is the last layer, mantle is the penultimate layer. Watch out, radius represents the layers 
	 *from center to surface of the Earth:*/
	IssmDouble lithosphere_thickness = litho->radius[numlayers] - litho->radius[numlayers-1];
	IssmDouble lithosphere_shear_modulus = litho->lame_mu[numlayers-1];
	IssmDouble lithosphere_density = litho->density[numlayers-1];
	IssmDouble mantle_shear_modulus = litho->lame_mu[numlayers-2];
	IssmDouble mantle_density = litho->density[numlayers-2];
	IssmDouble mantle_viscosity = litho->viscosity[numlayers-2];

	/*early return if we are NOT on an icy element:*/
	if(!IsIceInElement()) return;

	/*pull thickness averages! */
	IssmDouble *hes      = NULL;
	IssmDouble *times    = NULL;
	int         numtimes;
	this->GetInputAveragesUpToCurrentTime(TransientAccumulatedDeltaIceThicknessEnum,&hes,&times,&numtimes,currenttime);

	/*pull area of this Tria: */
	IssmDouble area=this->GetArea();

	/*element radius: */
	IssmDouble re=sqrt(area/M_PI);

	/*figure out gravity center of our element: */
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	IssmDouble x0=(xyz_list[0][0]+xyz_list[1][0]+xyz_list[2][0])/3.0;
	IssmDouble y0=(xyz_list[0][1]+xyz_list[1][1]+xyz_list[2][1])/3.0;

	/*start loading GiaDeflectionCore arguments: */
	arguments.re=re;
	arguments.hes=hes;
	arguments.times=times;
	arguments.numtimes=numtimes;
	arguments.currenttime=currenttime;
	arguments.lithosphere_shear_modulus=lithosphere_shear_modulus;
	arguments.lithosphere_density=lithosphere_density;
	arguments.mantle_shear_modulus=mantle_shear_modulus;
	arguments.mantle_viscosity=mantle_viscosity;
	arguments.mantle_density=mantle_density;
	arguments.lithosphere_thickness=lithosphere_thickness;
	arguments.rho_ice=rho_ice;
	arguments.idisk=this->id;
	arguments.iedge=cross_section_shape;
	arguments.yts=yts;

	for(int i=0;i<gsize;i++){
		/*compute distance from the center of the tria to the vertex i: */
		IssmDouble xi=x[i];
		IssmDouble yi=y[i];
		IssmDouble ri=sqrt(pow(xi-x0,2)+pow(yi-y0,2));

		/*load ri onto arguments for this vertex i: */
		arguments.ri=ri;

		/*for this Tria, compute contribution to rebound at vertex i: */
		GiaDeflectionCorex(&wi,&dwidt,&arguments);

		/*plug value into solution vector: */
		wg->SetValue(i,wi,ADD_VAL);
		dwgdt->SetValue(i,dwidt,ADD_VAL);
	}

	/*Free resources: */
	xDelete<IssmDouble>(hes);
	xDelete<IssmDouble>(times);

	return;
}
/*}}}*/
void       Tria::SealevelchangeGeometryInitial(IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae, int* lids, int* n_activevertices){ /*{{{*/

	/*Declarations:{{{*/
	int nel;
	IssmDouble area,planetarea,planetradius;
	IssmDouble constant,ratioe;
	IssmDouble rho_earth;
	IssmDouble NewtonG;
	IssmDouble g, cent_scaling;
	IssmDouble lati,longi;
	IssmDouble latitude[NUMVERTICES];
	IssmDouble longitude[NUMVERTICES];
	IssmDouble x,y,z,dx,dy,dz,N_azim,E_azim;
	IssmDouble xyz_list[NUMVERTICES][3];

	/*viscous stacks:*/
	IssmDouble* viscousRSL = NULL;
	IssmDouble* viscousU = NULL;
	IssmDouble* viscousN = NULL;
	IssmDouble* viscousE = NULL;
	IssmDouble* G_gravi_precomputed=NULL;

	/*viscoelastic green function:*/
	int index;
	int M;
	IssmDouble degacc;
	IssmDouble doubleindex,lincoef;

	/*Computational flags:*/
	bool computeselfattraction = false;
	bool computeelastic = false;
	bool computerotation = false;
	bool computeviscous = false;
	int  horiz;
	bool istime=true;
	IssmDouble timeacc=0.;
	IssmDouble start_time,final_time;
	int  nt,precomputednt;
	int grd, grdmodel;

	/*Rotational:*/
	IssmDouble* Grot=NULL;
	IssmDouble* GUrot=NULL;
	IssmDouble* GNrot=NULL;
	IssmDouble* GErot=NULL;
	IssmDouble* tide_love_h  = NULL;
	IssmDouble* tide_love_k  = NULL;
	IssmDouble* tide_love_l  = NULL;
	IssmDouble* LoveRotRSL   = NULL;
	IssmDouble* LoveRotU     = NULL;
	IssmDouble* LoveRothoriz = NULL;
	int* AlphaIndex   = NULL;
	int* AzimuthIndex = NULL;

	IssmDouble  moi_e, moi_p, omega;
	IssmDouble  Y21cos     , Y21sin     , Y20;
	IssmDouble dY21cos_dlat,dY21sin_dlat,dY20_dlat;
	IssmDouble dY21cos_dlon,dY21sin_dlon;
	IssmDouble polenudge;
	/*}}}*/

	/*Recover parameters:{{{ */
	rho_earth=FindParam(MaterialsEarthDensityEnum);
	this->parameters->FindParam(&computeselfattraction,SolidearthSettingsSelfAttractionEnum);
	this->parameters->FindParam(&computeelastic,SolidearthSettingsElasticEnum);
	this->parameters->FindParam(&computerotation,SolidearthSettingsRotationEnum);
	this->parameters->FindParam(&computeviscous,SolidearthSettingsViscousEnum);
	this->parameters->FindParam(&nel,MeshNumberofelementsEnum);
	this->parameters->FindParam(&planetarea,SolidearthPlanetAreaEnum);
	this->parameters->FindParam(&planetradius,SolidearthPlanetRadiusEnum);
	this->parameters->FindParam(&NewtonG,ConstantsNewtonGravityEnum);
	this->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	this->parameters->FindParam(&grd,SolidearthSettingsGRDEnum); 
	this->parameters->FindParam(&grdmodel,GrdModelEnum);

	/*early return:*/
	if (!grd || grdmodel!=ElasticEnum) return; //Love numbers won't be found in this case, return before loading them
	if(!computeselfattraction)return;

	if(computerotation){
		parameters->FindParam(&moi_e,RotationalEquatorialMoiEnum);
		parameters->FindParam(&moi_p,RotationalPolarMoiEnum);
		parameters->FindParam(&omega,RotationalAngularVelocityEnum);
		//parameters->FindParam(&tide_love_h,NULL,NULL,SealevelchangeTidalH2Enum);
		//parameters->FindParam(&tide_love_k,NULL,NULL,SealevelchangeTidalK2Enum);
		//if(horiz) parameters->FindParam(&tide_love_l,NULL,NULL,SealevelchangeTidalL2Enum);
	}
	/*}}}*/

	/*Recover precomputed green function kernels:{{{*/
	parameters->FindParam(&degacc,SolidearthSettingsDegreeAccuracyEnum);
	M=reCast<int,IssmDouble>(180.0/degacc+1.);

	DoubleVecParam* parameter;
	if(computeelastic){
		if(computerotation){
			parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeTidalH2Enum)); _assert_(parameter);
			parameter->GetParameterValueByPointer((IssmDouble**)&tide_love_h,NULL);

			parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeTidalK2Enum)); _assert_(parameter);
			parameter->GetParameterValueByPointer((IssmDouble**)&tide_love_k,NULL);

			if (horiz) {
				parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeTidalL2Enum)); _assert_(parameter);
				parameter->GetParameterValueByPointer((IssmDouble**)&tide_love_l,NULL);
			}
		}
	}
	/*}}}*/
	/*Compute lat long of all vertices in the element:{{{*/
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	for(int i=0;i<NUMVERTICES;i++){
		latitude[i]= asin(xyz_list[i][2]/planetradius);
		if((xyz_list[i][2]/planetradius)==1.0)latitude[i]=M_PI/2;
		longitude[i]= atan2(xyz_list[i][1],xyz_list[i][0]);
	}
	/*}}}*/
	/*Compute green functions:{{{ */
	if(computeviscous){
		this->parameters->FindParam(&istime,LoveIsTimeEnum);
		if(!istime)_error_("Frequency love numbers not supported yet!");
		this->parameters->FindParam(&timeacc,SolidearthSettingsTimeAccEnum);
		this->parameters->FindParam(&start_time,TimesteppingStartTimeEnum);
		this->parameters->FindParam(&final_time,TimesteppingFinalTimeEnum);
		nt=reCast<int,IssmDouble>((final_time-start_time)/timeacc)+1;
	}
	else{
		nt=1; //in elastic, or if we run only selfattraction, we need only one step
	}

	AlphaIndex=xNewZeroInit<int>(n_activevertices[this->lid]*nel);
	if (horiz) AzimuthIndex=xNewZeroInit<int>(n_activevertices[this->lid]*nel);
	int intmax=pow(2,16)-1;

	int* activevertices=xNew<int>(n_activevertices[this->lid]);

	int av=0;

	for (int i=0;i<3;i++){
		if(lids[this->vertices[i]->lid]==this->lid){
			activevertices[av]=i;
			for(int e=0;e<nel;e++){
				IssmDouble alpha;
				IssmDouble delPhi,delLambda;
				/*recovers info for this element and vertex:*/
				IssmDouble late= asin(zze[e]/sqrt( pow(xxe[e],2.0)+ pow(yye[e],2.0)+ pow(zze[e],2.0)));
				IssmDouble longe= atan2(yye[e],xxe[e]);

				lati=latitude[i];
				longi=longitude[i];

				/*Computes alpha angle between centroid and current vertex, and indexes alpha in precomputed tables: */
				delPhi=fabs(lati-late); delLambda=fabs(longi-longe); if (delLambda>M_PI)delLambda=2*M_PI-delLambda;
				alpha=2.*asin(sqrt(pow(sin(delPhi/2),2)+cos(lati)*cos(late)*pow(sin(delLambda/2),2)));
				doubleindex=alpha/M_PI*reCast<IssmDouble,int>(M-1); //maps 0<alpha<PI on [0:M-1]
				index=reCast<int,IssmDouble>(doubleindex); //truncates doubleindex to integer part
				if ((doubleindex-index)>=0.5) index+=1; //nearest neighbour
				_assert_(index>=0 && index<M);

				if(horiz){
					/*Compute azimuths*/
					dx=cos(lati)*sin(late)-sin(lati)*cos(late)*cos(longe-longi);
					dy=sin(longe-longi)*cos(late);
					//angle between horiz motion and North, remapped from a double on [0,2*pi] to a int [0,intmax]
					AzimuthIndex[av*nel+e]=reCast<int,IssmDouble>(intmax*(atan2(dy,dx)/2/M_PI));
				}
				AlphaIndex[av*nel+e]=index;
			}
			av+=1;			
		} //for (int i=0;i<3;i++)
	} //for(int e=0;e<nel;e++)

	/*Add in inputs:*/
	this->inputs->SetIntArrayInput(SealevelchangeConvolutionVerticesEnum,this->lid,activevertices,n_activevertices[this->lid]);
	this->inputs->SetIntArrayInput(SealevelchangeAlphaIndexEnum,this->lid,AlphaIndex,nel*n_activevertices[this->lid]);
	if(horiz) this->inputs->SetIntArrayInput(SealevelchangeAzimuthIndexEnum,this->lid,AzimuthIndex,nel*n_activevertices[this->lid]);

	/*}}}*/
	/*Compute rotation kernel:{{{*/
	if(computerotation){
		//initialization
		LoveRotRSL  = xNewZeroInit<IssmDouble>(nt);
		LoveRotU    = xNewZeroInit<IssmDouble>(nt);
		if(horiz)LoveRothoriz= xNewZeroInit<IssmDouble>(nt);
		Grot        = xNewZeroInit<IssmDouble>(3*3*nt); //3 polar motion components * 3 vertices * number of time steps
		GUrot       = xNewZeroInit<IssmDouble>(3*3*nt);

		if (horiz){
			GErot=xNewZeroInit<IssmDouble>(3*3*nt);
			GNrot=xNewZeroInit<IssmDouble>(3*3*nt);
		}

		/*What is the gravity at planet's surface: */
		g=4.0/3.0*M_PI*rho_earth*NewtonG*planetradius;
		cent_scaling=pow(omega*planetradius,2.0); //centrifugal potential dimensioning constant
		for(int t=0;t<nt;t++){
			//Amplitude of the rotational feedback
			//to speed up calculation we include the dimension constant r^2*Omega^2/g, so all 3 of these are in meters
			LoveRotRSL[t]=((1.0+tide_love_k[t]-tide_love_h[t])/g)*cent_scaling;
			LoveRotU[t]=(tide_love_h[t]/g)*cent_scaling;
			if (horiz) LoveRothoriz[t]=(tide_love_l[t]/g)*cent_scaling;
		}
		for(int i=0;i<3;i++){

			//Avoid singularity of the poles
			if (latitude[i]==M_PI/2.){
				//North pole: nudge south
				polenudge=-1.e-12;
			}
			else if (latitude[i]==-M_PI/2.){
				//South pole: nudge north
				polenudge=1.e-12;
			}
			else {
				polenudge=0.0;
			}

			lati=latitude[i]+polenudge;
			longi=longitude[i];

			//Spherical harmonic functions of degree 2 (spatial pattern of rotation)
			Y21cos= -0.5*sin(2.*lati)*cos(longi);
			Y21sin= -0.5*sin(2.*lati)*sin(longi);
			Y20   = -(1.0/6.0 - 0.5*cos(2.0*lati));

			if (computeelastic && horiz){
				//bed_N = Love_l * d(Y)/dlat ;
				dY21cos_dlat=-cos(2.*lati)*cos(longi);
				dY21sin_dlat=-cos(2.*lati)*sin(longi);
				dY20_dlat=   -sin(2.*lati);

				//bed_E = Love_l * 1/cos(lat) * d(Y)/dlon ;
				dY21cos_dlon=-Y21sin/cos(lati);
				dY21sin_dlon=Y21cos/cos(lati);
				//dY20_dlon=0.;
			}

			for(int t=0;t<nt;t++){

				Grot[0*3*nt+i*nt+t]= LoveRotRSL[t]*Y21cos; //x component of polar motion
				Grot[1*3*nt+i*nt+t]= LoveRotRSL[t]*Y21sin; //y
				Grot[2*3*nt+i*nt+t]= LoveRotRSL[t]*Y20;    //z

				if (computeelastic){
					GUrot[0*3*nt+i*nt+t]= LoveRotU[t]*Y21cos;
					GUrot[1*3*nt+i*nt+t]= LoveRotU[t]*Y21sin;
					GUrot[2*3*nt+i*nt+t]= LoveRotU[t]*Y20;
					if (horiz){
						//bed_N = Love_l * d(Y)/dlat ;
						GNrot[0*3*nt+i*nt+t]= LoveRothoriz[t]*dY21cos_dlat;
						GNrot[1*3*nt+i*nt+t]= LoveRothoriz[t]*dY21sin_dlat;
						GNrot[2*3*nt+i*nt+t]= LoveRothoriz[t]*dY20_dlat;

						//bed_E = Love_l * 1/cos(lat) * d(Y)/dlon ;
						GErot[0*3*nt+i*nt+t]= LoveRothoriz[t]*dY21cos_dlon;
						GErot[1*3*nt+i*nt+t]= LoveRothoriz[t]*dY21sin_dlon;
						GErot[2*3*nt+i*nt+t]= 0.0;
					}
				}
			}
		}
		this->inputs->SetArrayInput(SealevelchangeGrotEnum,this->lid,Grot,3*3*nt);
		if (computeelastic){
			this->inputs->SetArrayInput(SealevelchangeGUrotEnum,this->lid,GUrot,3*3*nt);
			if(horiz){
				this->inputs->SetArrayInput(SealevelchangeGNrotEnum,this->lid,GNrot,3*3*nt);
				this->inputs->SetArrayInput(SealevelchangeGErotEnum,this->lid,GErot,3*3*nt);
			}
		}
		/*Free resources:*/
		xDelete<IssmDouble>(LoveRotRSL);
		xDelete<IssmDouble>(LoveRotU);
		if(horiz)xDelete<IssmDouble>(LoveRothoriz);
	}
	/*}}}*/
	/*Initialize viscous stacks: {{{*/
	if(computeviscous){
		viscousRSL=xNewZeroInit<IssmDouble>(3*nt);
		viscousU=xNewZeroInit<IssmDouble>(3*nt);

		this->inputs->SetArrayInput(SealevelchangeViscousRSLEnum,this->lid,viscousRSL,3*nt);
		this->inputs->SetArrayInput(SealevelchangeViscousUEnum,this->lid,viscousU,3*nt);
		this->parameters->SetParam(0,SealevelchangeViscousIndexEnum);
		if(horiz){
			viscousN=xNewZeroInit<IssmDouble>(3*nt);
			viscousE=xNewZeroInit<IssmDouble>(3*nt);
			this->inputs->SetArrayInput(SealevelchangeViscousNEnum,this->lid,viscousN,3*nt);
			this->inputs->SetArrayInput(SealevelchangeViscousEEnum,this->lid,viscousE,3*nt);
		}
	}
	/*}}}*/

	/*Free allocations:{{{*/
	xDelete<int>(activevertices);
	xDelete<int>(AlphaIndex);
	if(horiz){
		xDelete<int>(AzimuthIndex);
	}
	if(computerotation){
		xDelete<IssmDouble>(Grot);
		xDelete<IssmDouble>(GUrot);
		if (horiz){
			xDelete<IssmDouble>(GNrot);
			xDelete<IssmDouble>(GErot);
		}
	}
	/*}}}*/
	return;

}
/*}}}*/
void       Tria::SealevelchangeGeometrySubElementKernel(SealevelGeometry* slgeom){ /*{{{*/

	/*Declarations:{{{*/
	int nel;
	IssmDouble planetarea,planetradius;
	IssmDouble constant,ratioe;
	IssmDouble rho_earth;
	IssmDouble lati,longi;
	IssmDouble latitude[NUMVERTICES];
	IssmDouble longitude[NUMVERTICES];
	IssmDouble x,y,z,dx,dy,dz,N_azim,E_azim;
	IssmDouble xyz_list[NUMVERTICES][3];
	int* activevertices = NULL;
	int n_activevertices, av;
	int** AlphaIndex=NULL;
	int** AzimIndex=NULL;

	/*viscoelastic green function:*/
	int index;
	int M;
	IssmDouble doubleindex,lincoef, degacc;

	/*Computational flags:*/
	bool computeselfattraction = false;
	bool computeelastic = false;
	bool computeviscous = false;
	int  horiz;
	int grd, grdmodel;

	bool istime=true;
	IssmDouble timeacc=0;
	IssmDouble start_time,final_time;
	int  nt,precomputednt;
	int intmax=pow(2,16)-1;

	/*}}}*/
	/*Recover parameters:{{{ */
	rho_earth=FindParam(MaterialsEarthDensityEnum);
	this->parameters->FindParam(&computeselfattraction,SolidearthSettingsSelfAttractionEnum);
	this->parameters->FindParam(&computeelastic,SolidearthSettingsElasticEnum);
	this->parameters->FindParam(&computeviscous,SolidearthSettingsViscousEnum);
	this->parameters->FindParam(&nel,MeshNumberofelementsEnum);
	this->parameters->FindParam(&planetarea,SolidearthPlanetAreaEnum);
	this->parameters->FindParam(&planetradius,SolidearthPlanetRadiusEnum);
	this->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	this->parameters->FindParam(&grd,SolidearthSettingsGRDEnum); 
	this->parameters->FindParam(&grdmodel,GrdModelEnum);
	/*}}}*/

	/*early return:*/
	if (!grd || grdmodel!=ElasticEnum) return; //Love numbers won't be found in this case, return before loading them
	if(!computeselfattraction)return;

	/*Recover precomputed green function kernels:{{{*/
	parameters->FindParam(&degacc,SolidearthSettingsDegreeAccuracyEnum);
	M=reCast<int,IssmDouble>(180.0/degacc+1.);

	/*}}}*/
	/*Compute lat long of all vertices in the element:{{{*/
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	for(int i=0;i<NUMVERTICES;i++){
		latitude[i]= asin(xyz_list[i][2]/planetradius);
		longitude[i]= atan2(xyz_list[i][1],xyz_list[i][0]);
	}
	/*}}}*/
	/*Compute green functions:{{{ */

	if(computeviscous){
		this->parameters->FindParam(&istime,LoveIsTimeEnum);
		if(!istime)_error_("Frequency love numbers not supported yet!");
		this->parameters->FindParam(&timeacc,SolidearthSettingsTimeAccEnum);
		this->parameters->FindParam(&start_time,TimesteppingStartTimeEnum);
		this->parameters->FindParam(&final_time,TimesteppingFinalTimeEnum);
		nt=reCast<int,IssmDouble>((final_time-start_time)/timeacc)+1;
	}
	else{
		nt=1; //in elastic, or if we run only selfattraction, we need only one step
	}
	AlphaIndex=xNew<int*>(SLGEOM_NUMLOADS);
	if(horiz) AzimIndex=xNew<int*>(SLGEOM_NUMLOADS);

	this->inputs->GetIntArrayPtr(SealevelchangeConvolutionVerticesEnum,this->lid,&activevertices,&n_activevertices);
	// 0<=n_activevertices<=3 is the number of vertices this element is in charge of computing fields in during the sea level convolutions
	// activevertices contains the vertex indices (1,2 and/or 3) in case debugging is required, they are supposed to appear in the same order as slgeom->lids

	//Allocate: 
	for(int l=0;l<SLGEOM_NUMLOADS;l++){
		int nbar=slgeom->nbar[l];
		AlphaIndex[l]=xNewZeroInit<int>(n_activevertices*nbar);
		if(horiz) AzimIndex[l]=xNewZeroInit<int>(n_activevertices*nbar);

		//av=0;
		//for (int i=0;i<3;i++){
		for (int av=0;av<n_activevertices;av++){
			//if(slgeom->lids[this->vertices[i]->lid]==this->lid){
			int i=activevertices[av];
			for(int e=0;e<nbar;e++){
				IssmDouble alpha;
				IssmDouble delPhi,delLambda;
				/*recover info for this element and vertex:*/
				IssmDouble late= slgeom->latbarycentre[l][e]; 
				IssmDouble longe= slgeom->longbarycentre[l][e]; 
				late=late/180*M_PI;
				longe=longe/180*M_PI;
				lati=latitude[i];
				longi=longitude[i];
					if(horiz){
					/*Compute azimuths*/
						dx=cos(lati)*sin(late)-sin(lati)*cos(late)*cos(longe-longi);
						dy=sin(longe-longi)*cos(late);
						//angle between horiz motion and North, remapped from a double on [0,2*pi] to a int [0,intmax]
						AzimIndex[l][av*nbar+e]=reCast<int,IssmDouble>(intmax*(atan2(dy,dx)/2/M_PI));
					}

				/*Compute alpha angle between centroid and current vertex and index into precomputed tables: */
				delPhi=fabs(lati-late); delLambda=fabs(longi-longe); if (delLambda>M_PI)delLambda=2*M_PI-delLambda;
				alpha=2.*asin(sqrt(pow(sin(delPhi/2.0),2.0)+cos(lati)*cos(late)*pow(sin(delLambda/2.0),2.0)));
				doubleindex=alpha/M_PI*reCast<IssmDouble,int>(M-1); //maps 0<alpha<PI on [0:M-1]
				index=reCast<int,IssmDouble>(doubleindex); //truncates doubleindex to integer part

				if ((doubleindex-index)>=0.5) index+=1; //nearest neighbour
				if (index==M-1){ //avoids out of bound case
					index-=1;
					lincoef=1;
				}
				AlphaIndex[l][av*nbar+e]=index;
			//}
			//av+=1;
			}

		}
	}

	/*Save all these arrayins for each element:*/
	for (int l=0;l<SLGEOM_NUMLOADS;l++){
		this->inputs->SetIntArrayInput(slgeom->AlphaIndexEnum(l),this->lid,AlphaIndex[l],slgeom->nbar[l]*n_activevertices);
		if(horiz) this->inputs->SetIntArrayInput(slgeom->AzimuthIndexEnum(l),this->lid,AzimIndex[l],slgeom->nbar[l]*n_activevertices);
	}
	/*}}}*/
	/*Free memory:{{{*/
	for (int l=0;l<SLGEOM_NUMLOADS;l++){
		xDelete<int>(AlphaIndex[l]);
		if(horiz) xDelete<int>(AzimIndex[l]);
	}
	xDelete<int*>(AlphaIndex);
	if(horiz) xDelete<int*>(AzimIndex); 

	/*}}}*/
	return;

}
/*}}}*/
void       Tria::SealevelchangeGeometryCentroidLoads(SealevelGeometry* slgeom, IssmDouble* xxe, IssmDouble* yye, IssmDouble* zze, IssmDouble* areae){ /*{{{*/

	/* Classic buildup of load weights, centroids and areas *for elements which are fully inside a mask. 
	 * At the same time, we'll tag the elements that are fractionally only inside a mask*/

	IssmDouble loadweights[3]={0};
	IssmDouble area;
	IssmDouble loadweightsocean[3]; //to keep memory of these loads, no need to recompute for bottom pressure.
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble planetradius;
	IssmDouble late,longe;

	/*flags:*/
	bool isocean=false;
	bool isoceanonly=false;
	bool isice=false;
	bool isiceonly=false;
	bool computeice=false;
	bool computebp=false;
	bool computehydro=false;
	bool ismasstransport=false;
	bool ismmemasstransport=false;

	/*constants:*/
	IssmDouble constant=0;

	/*recover parameters:*/
	this->parameters->FindParam(&ismasstransport,TransientIsmasstransportEnum);
	this->parameters->FindParam(&ismmemasstransport,TransientIsmmemasstransportEnum);
	if(ismasstransport || ismmemasstransport)computeice=true;
	this->parameters->FindParam(&computebp,TransientIsoceantransportEnum);
	this->parameters->FindParam(&computehydro,TransientIshydrologyEnum);
	this->parameters->FindParam(&planetradius,SolidearthPlanetRadiusEnum);

	/*get vertex information:*/
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);

	/*answer mask questions:*/
	isiceonly=this->IsIceOnlyInElement();
	isice=this->IsIceInElement();
	isoceanonly=this->IsOceanOnlyInElement();
	isocean=this->IsOceanInElement();
	slgeom->isoceanin[this->lid]=isocean; //keep track for later.
	area=areae[this->sid];

	/*Compute element ids, used to speed up computations in convolution phase:{{{*/
	for(int i=0;i<NUMVERTICES;i++){
		slgeom->lids[this->vertices[i]->lid]=this->lid;
	}
	/*}}}*/

	/*set barycentre for all elements, to be updated for fractional loads in the next routine: */
	//late= asin(zze[this->sid]/planetradius)*180.0/M_PI;
	late= asin(zze[this->sid]/sqrt( pow(xxe[this->sid],2.0)+ pow(yye[this->sid],2.0)+ pow(zze[this->sid],2.0)))*180.0/M_PI;
	longe= atan2(yye[this->sid],xxe[this->sid])*180.0/M_PI;
	slgeom->longe[this->lid]=longe;
	slgeom->late[this->lid]=late;

	/*compute areas and load weights for ocean and flag elements only partially in the ocean:*/
	if(isoceanonly){ 
		slgeom->LoadArea[SLGEOM_OCEAN][this->lid]=area;
		for(int i=0;i<NUMVERTICES;i++) slgeom->LoadWeigths[SLGEOM_OCEAN][i][this->lid]=1.0/3.0;

		#ifdef _ISSM_DEBUG_ /*{{{*/
		/*Inform mask: */
		constant=1.0;
		for(int i=0;i<NUMVERTICES;i++) loadweightsocean[i]=1.0/3.0;
		this->AddInput(SealevelBarystaticOceanMaskEnum,&constant,P0Enum); 
		this->AddInput(SealevelBarystaticOceanWeightsEnum,loadweightsocean,P1DGEnum);
		this->AddInput(SealevelBarystaticOceanAreaEnum,&area,P0Enum);
		#endif /*}}}*/
	}
	else if(!isocean){
		slgeom->LoadArea[SLGEOM_OCEAN][this->lid]=0;
		for(int i=0;i<NUMVERTICES;i++) slgeom->LoadWeigths[SLGEOM_OCEAN][i][this->lid]=0.0;
		#ifdef _ISSM_DEBUG_ /*{{{*/
		/*Inform mask: */
		constant=0.0;
		for(int i=0;i<NUMVERTICES;i++) loadweightsocean[i]=0.0;
		this->AddInput(SealevelBarystaticOceanMaskEnum,&constant,P0Enum); 
		this->AddInput(SealevelBarystaticOceanWeightsEnum,loadweightsocean,P1DGEnum);
		this->AddInput(SealevelBarystaticOceanAreaEnum,&constant,P0Enum);
		#endif /*}}}*/
	}
	else{
		slgeom->issubelement[SLGEOM_OCEAN][this->lid]=true;
		slgeom->nsubel[SLGEOM_OCEAN]++;
	}

	/*early return if we are not on an ice sheet , and we are not requesting 
	 *hydrology or bottom pressure loads :*/
	if(!computebp && !computehydro){
		if(!isice  || isoceanonly) {
			#ifdef _ISSM_DEBUG_
			constant=0; 
			this->AddInput(SealevelBarystaticIceMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticIceAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticIceWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticHydroMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticHydroWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticHydroAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticBpAreaEnum,&constant,P0Enum);
			#endif
			for(int i=0;i<NUMVERTICES;i++){
				slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]=0;
				slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]=0;
			}
			slgeom->LoadArea[SLGEOM_ICE][this->lid]=0;
			slgeom->LoadArea[SLGEOM_WATER][this->lid]=0;
			return;
		}
	}

	/*early return if we are fully floating and we are not doing bottom pressure loads:*/
	if(!computebp){
		if (isoceanonly){
			#ifdef _ISSM_DEBUG_
			constant=0;
			this->AddInput(SealevelBarystaticIceMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticIceWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticIceAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticHydroMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticHydroWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticHydroAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticBpAreaEnum,&constant,P0Enum);
			#endif
			for(int i=0;i<NUMVERTICES;i++){
				slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]=0;
				slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]=0;
			}
			slgeom->LoadArea[SLGEOM_ICE][this->lid]=0;
			slgeom->LoadArea[SLGEOM_WATER][this->lid]=0;
			return;
		}
	}

	/*early return if we are not on the ocean and we are not doing ice mass transport of 
	 * hydrology:*/
	if(!computeice  && !computehydro){
		if(!isocean){
			#ifdef _ISSM_DEBUG_
			constant=0;
			this->AddInput(SealevelBarystaticIceMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticIceWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticIceAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticHydroMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticHydroWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticHydroAreaEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpMaskEnum,&constant,P0Enum);
			this->AddInput(SealevelBarystaticBpWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticBpAreaEnum,&constant,P0Enum);
			#endif
			for(int i=0;i<NUMVERTICES;i++){
				slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]=0;
				slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]=0;
			}
			slgeom->LoadArea[SLGEOM_ICE][this->lid]=0;
			slgeom->LoadArea[SLGEOM_WATER][this->lid]=0;
			return;
		}
	}

	/*Deal with ice loads if we are on grounded ice:*/
	if(isice && !isoceanonly && computeice){
		if(isiceonly && !isocean){
			slgeom->LoadArea[SLGEOM_ICE][this->lid]=area;
			for(int i=0;i<NUMVERTICES;i++) slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]=1.0/3.0;

			#ifdef _ISSM_DEBUG_ /*{{{*/
			/*Inform mask: */
			constant=1.0;
			for(int i=0;i<NUMVERTICES;i++) loadweights[i]=1.0/3.0;
			this->AddInput(SealevelBarystaticIceMaskEnum,&constant,P0Enum); 
			this->AddInput(SealevelBarystaticIceWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticIceAreaEnum,&area,P0Enum);
			#endif /*}}}*/
		}
		else{
			slgeom->issubelement[SLGEOM_ICE][this->lid]=true;
			slgeom->nsubel[SLGEOM_ICE]++;
		}
	} 

	/*Deal with water loads if we are on ground:*/
	if(!isoceanonly && computehydro){

		if(!isocean){
			slgeom->LoadArea[SLGEOM_WATER][this->lid]=area;
			for(int i=0;i<NUMVERTICES;i++) slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]=1.0/3.0;

			#ifdef _ISSM_DEBUG_ /*{{{*/
			/*Inform mask: */
			constant=1.0;
			for(int i=0;i<NUMVERTICES;i++) loadweights[i]=1.0/3.0;
			this->AddInput(SealevelBarystaticHydroMaskEnum,&constant,P0Enum); 
			this->AddInput(SealevelBarystaticHydroWeightsEnum,loadweights,P1DGEnum);
			this->AddInput(SealevelBarystaticHydroAreaEnum,&area,P0Enum);
			#endif /*}}}*/
		}
		else{
			slgeom->issubelement[SLGEOM_WATER][this->lid]=true;
			slgeom->nsubel[SLGEOM_WATER]++;
		}
	}

}
/*}}}*/
void       Tria::SealevelchangeBarystaticLoads(GrdLoads* loads,  BarystaticContributions* barycontrib, SealevelGeometry* slgeom){ /*{{{*/

	int nel;

	/*Inputs:*/
	IssmDouble I[NUMVERTICES]; 
	IssmDouble W[NUMVERTICES];
	IssmDouble BP[NUMVERTICES];
	IssmDouble* areae=NULL;

	/*output: */
	IssmDouble bslcice=0;
	IssmDouble bslchydro=0;
	IssmDouble bslcbp=0;
	IssmDouble BPavg=0;
	IssmDouble Iavg=0;
	IssmDouble Wavg=0;

	/*ice properties: */
	IssmDouble rho_ice,rho_water,rho_freshwater;

	/*recover some parameters:*/
	this->parameters->FindParam(&rho_ice,MaterialsRhoIceEnum);
	this->parameters->FindParam(&rho_water,MaterialsRhoSeawaterEnum);
	this->parameters->FindParam(&rho_freshwater,MaterialsRhoFreshwaterEnum);
	this->parameters->FindParam(&areae,&nel,AreaeEnum);

	/*Retrieve inputs:*/
	Element::GetInputListOnVertices(&I[0],DeltaIceThicknessEnum);
	Element::GetInputListOnVertices(&W[0],DeltaTwsEnum);
	Element::GetInputListOnVertices(&BP[0],DeltaBottomPressureEnum);

	for(int i=0;i<NUMVERTICES;i++){
		Iavg+=I[i]*slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]*slgeom->LoadArea[SLGEOM_ICE][this->lid];
		Wavg+=W[i]*slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]*slgeom->LoadArea[SLGEOM_WATER][this->lid];
		BPavg+=BP[i]*slgeom->LoadWeigths[SLGEOM_OCEAN][i][this->lid]*slgeom->LoadArea[SLGEOM_OCEAN][this->lid];
	}

	/*convert from m^3 to kg:*/
	Iavg*=rho_ice;
	Wavg*=rho_freshwater;
	BPavg*=rho_water;

	#ifdef _ISSM_DEBUG_ 
	this->AddInput(SealevelBarystaticIceLoadEnum,&Iavg,P0Enum);
	this->AddInput(SealevelBarystaticHydroLoadEnum,&Wavg,P0Enum);
	this->AddInput(SealevelBarystaticBpLoadEnum,&BPavg,P0Enum);
	#endif

	/*Compute barystatic component in kg:*/
	// Note: Iavg, etc, already include partial area factor phi for subelement loading
	bslcice =   -Iavg;
	bslchydro = -Wavg;
	bslcbp =    -BPavg;

	_assert_(!xIsNan<IssmDouble>(bslcice));
	_assert_(!xIsNan<IssmDouble>(bslchydro));
	_assert_(!xIsNan<IssmDouble>(bslcbp));

	/*Plug values into subelement load vector:*/
	if(slgeom->issubelement[SLGEOM_ICE][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_ICE][this->lid];
		loads->vsubloads[SLGEOM_ICE]->SetValue(intj,Iavg,INS_VAL);
		Iavg=0; //avoid double counting centroid loads and subelement loads!
	}
	if(slgeom->issubelement[SLGEOM_WATER][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_WATER][this->lid];
		loads->vsubloads[SLGEOM_WATER]->SetValue(intj,Wavg,INS_VAL);
		Wavg=0;
	}
	if(slgeom->issubelement[SLGEOM_OCEAN][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_OCEAN][this->lid];
		loads->vsubloads[SLGEOM_OCEAN]->SetValue(intj,BPavg,INS_VAL); 
		BPavg=0;
	}
	/*Plug remaining values into centroid load vector:*/
	loads->vloads->SetValue(this->sid,Iavg+Wavg+BPavg,INS_VAL);

	/*Keep track of barystatic contributions:*/
	barycontrib->Set(this->Sid(),bslcice,bslchydro,bslcbp);

	/*Free resources*/
	xDelete<IssmDouble>(areae);

}/*}}}*/
void       Tria::SealevelchangeGeometrySubElementLoads(SealevelGeometry* slgeom, IssmDouble* areae){ /*{{{*/

	/* Classic buildup of load weights, centroids and areas *for elements which are fully inside a mask. 
	 * At the same time, we'll tag the elements that are fractionally only inside a mask*/

	IssmDouble loadweights[3]={0};
	IssmDouble area,loadarea;
	IssmDouble loadareaocean;
	IssmDouble loadweightsocean[3]; //to keep memory of these loads, no need to recompute for bottom pressure.
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble latbar=slgeom->late[this->lid];
	IssmDouble longbar=slgeom->longe[this->lid];
	IssmDouble constant;
	IssmDouble nanconstant=NAN;

	/*get vertex and area information:*/
	::GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	area=areae[this->sid];

	#ifdef _ISSM_DEBUG_
	this->AddInput(SealevelBarystaticIceLatbarEnum,&latbar,P0Enum); 
	this->AddInput(SealevelBarystaticIceLongbarEnum,&longbar,P0Enum); 
	this->AddInput(SealevelBarystaticHydroLatbarEnum,&latbar,P0Enum); 
	this->AddInput(SealevelBarystaticHydroLongbarEnum,&longbar,P0Enum); 
	this->AddInput(SealevelBarystaticOceanLatbarEnum,&latbar,P0Enum); 
	this->AddInput(SealevelBarystaticOceanLongbarEnum,&longbar,P0Enum); 
	#endif

	if(slgeom->issubelement[SLGEOM_OCEAN][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_OCEAN][this->lid];

		this->GetNodalWeightsAndAreaAndCentroidsFromLeveset(&loadweightsocean[0],&loadareaocean,&latbar, &longbar, slgeom->late[this->lid], slgeom->longe[this->lid], area, MaskOceanLevelsetEnum);
		slgeom->LoadArea[SLGEOM_OCEAN][this->lid]=loadareaocean;
		slgeom->vareae_subel[SLGEOM_OCEAN]->SetValue(intj,loadareaocean,INS_VAL);
		slgeom->vlatbarycentre[SLGEOM_OCEAN]->SetValue(intj,latbar,INS_VAL);
		slgeom->vlongbarycentre[SLGEOM_OCEAN]->SetValue(intj,longbar,INS_VAL);

		for(int i=0;i<NUMVERTICES;i++) slgeom->LoadWeigths[SLGEOM_OCEAN][i][this->lid]=loadweightsocean[i];

		#ifdef _ISSM_DEBUG_ /*{{{*/
		/*Inform mask: */
		constant=loadareaocean/area;
		this->AddInput(SealevelBarystaticOceanMaskEnum,&constant,P0Enum); 
		this->AddInput(SealevelBarystaticOceanWeightsEnum,loadweightsocean,P1DGEnum);
		this->AddInput(SealevelBarystaticOceanAreaEnum,&loadareaocean,P0Enum);

		this->AddInput(SealevelBarystaticOceanLatbarEnum,&latbar,P0Enum); 
		this->AddInput(SealevelBarystaticOceanLongbarEnum,&longbar,P0Enum); 
		#endif /*}}}*/
	}
	if(slgeom->issubelement[SLGEOM_ICE][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_ICE][this->lid];

		this->GetNodalWeightsAndAreaAndCentroidsFromLeveset(&loadweights[0],&loadarea,&latbar, &longbar, slgeom->late[this->lid], slgeom->longe[this->lid], area, -MaskOceanLevelsetEnum,MaskIceLevelsetEnum);

		slgeom->LoadArea[SLGEOM_ICE][this->lid]=loadarea;
		slgeom->vareae_subel[SLGEOM_ICE]->SetValue(intj,loadarea,INS_VAL);
		slgeom->vlatbarycentre[SLGEOM_ICE]->SetValue(intj,latbar,INS_VAL);
		slgeom->vlongbarycentre[SLGEOM_ICE]->SetValue(intj,longbar,INS_VAL);

		for(int i=0;i<NUMVERTICES;i++)slgeom->LoadWeigths[SLGEOM_ICE][i][this->lid]=loadweights[i];

		#ifdef _ISSM_DEBUG_
		/*Inform mask: */
		constant=loadarea/area; 
		this->AddInput(SealevelBarystaticIceMaskEnum,&constant,P0Enum);
		this->AddInput(SealevelBarystaticIceWeightsEnum,loadweights,P1DGEnum);
		this->AddInput(SealevelBarystaticIceAreaEnum,&loadarea,P0Enum);

		this->AddInput(SealevelBarystaticIceLatbarEnum,&latbar,P0Enum); 
		this->AddInput(SealevelBarystaticIceLongbarEnum,&longbar,P0Enum); 

		#endif
	}
	if(slgeom->issubelement[SLGEOM_WATER][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_WATER][this->lid];

		this->GetNodalWeightsAndAreaAndCentroidsFromLeveset(&loadweights[0],&loadarea,&latbar, &longbar, slgeom->late[this->lid], slgeom->longe[this->lid], area, -MaskOceanLevelsetEnum);

		slgeom->LoadArea[SLGEOM_WATER][this->lid]=loadarea;
		slgeom->vareae_subel[SLGEOM_WATER]->SetValue(intj,loadarea,INS_VAL);
		slgeom->vlatbarycentre[SLGEOM_WATER]->SetValue(intj,latbar,INS_VAL);
		slgeom->vlongbarycentre[SLGEOM_WATER]->SetValue(intj,longbar,INS_VAL);

		for(int i=0;i<NUMVERTICES;i++)slgeom->LoadWeigths[SLGEOM_WATER][i][this->lid]=loadweights[i];

		#ifdef _ISSM_DEBUG_
		/*Inform mask: */
		constant=loadarea/area; 
		this->AddInput(SealevelBarystaticHydroMaskEnum,&constant,P0Enum);
		this->AddInput(SealevelBarystaticHydroWeightsEnum,loadweights,P1DGEnum);
		this->AddInput(SealevelBarystaticHydroAreaEnum,&loadarea,P0Enum);

		this->AddInput(SealevelBarystaticHydroLatbarEnum,&latbar,P0Enum); 
		this->AddInput(SealevelBarystaticHydroLongbarEnum,&longbar,P0Enum); 

		#endif
	}

}
/*}}}*/
void       Tria::SealevelchangeUpdateViscousFields(IssmDouble lincoeff, int newindex, int offset){ /*{{{*/

	/*Inputs:*/
	IssmDouble* viscousRSL=NULL;
	IssmDouble* viscousU=NULL;
	IssmDouble* viscousN=NULL;
	IssmDouble* viscousE=NULL;
	int         viscousnumsteps;
	int         size;
	bool        viscous=false;
	int	    horiz=0;

	this->parameters->FindParam(&viscous,SolidearthSettingsViscousEnum);

	if(viscous){
		this->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
		this->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);

		this->inputs->GetArrayPtr(SealevelchangeViscousRSLEnum,this->lid,&viscousRSL,&size);
		this->inputs->GetArrayPtr(SealevelchangeViscousUEnum,this->lid,&viscousU,&size);
		if(horiz){
			this->inputs->GetArrayPtr(SealevelchangeViscousNEnum,this->lid,&viscousN,&size);
			this->inputs->GetArrayPtr(SealevelchangeViscousEEnum,this->lid,&viscousE,&size);
		}

		for(int i=0;i<NUMVERTICES;i++){
			viscousRSL[i*viscousnumsteps+newindex+offset]=(1-lincoeff)*viscousRSL[i*viscousnumsteps+newindex]+lincoeff*viscousRSL[i*viscousnumsteps+newindex+1];
			viscousU[i*viscousnumsteps+newindex+offset]=(1-lincoeff)*viscousU[i*viscousnumsteps+newindex]+lincoeff*viscousU[i*viscousnumsteps+newindex+1];
			if(horiz){
				viscousN[i*viscousnumsteps+newindex+offset]=(1-lincoeff)*viscousN[i*viscousnumsteps+newindex]+lincoeff*viscousN[i*viscousnumsteps+newindex+1];
				viscousE[i*viscousnumsteps+newindex+offset]=(1-lincoeff)*viscousE[i*viscousnumsteps+newindex]+lincoeff*viscousE[i*viscousnumsteps+newindex+1];
			}
		}
	}
}
/*}}}*/
void       Tria::SealevelchangeOceanAverage(GrdLoads* loads, Vector<IssmDouble>* oceanareas, Vector<IssmDouble>* subelementoceanareas, IssmDouble* sealevelpercpu, SealevelGeometry* slgeom){ /*{{{*/

	IssmDouble oceanaverage=0;
	IssmDouble oceanarea=0;
	IssmDouble rho_water;

	this->parameters->FindParam(&rho_water,MaterialsRhoSeawaterEnum);

	/*retrieve ocean average and area:*/
	for(int i=0;i<NUMVERTICES;i++){
		oceanaverage+=sealevelpercpu[this->vertices[i]->lid]*slgeom->LoadWeigths[SLGEOM_OCEAN][i][this->lid];
	}

	oceanarea=slgeom->LoadArea[SLGEOM_OCEAN][this->lid];
	oceanaverage*=rho_water*oceanarea;

	/*add ocean average in the global sealevelloads vector:*/
	if(slgeom->issubelement[SLGEOM_OCEAN][this->lid]){
		int intj=slgeom->subelementmapping[SLGEOM_OCEAN][this->lid];
		loads->vsubsealevelloads->SetValue(intj,oceanaverage,INS_VAL);
		loads->vsealevelloads->SetValue(this->sid,0.,INS_VAL);
	}
	else loads->vsealevelloads->SetValue(this->sid,oceanaverage,INS_VAL);

	#ifdef _ISSM_DEBUG_ 
	this->AddInput(SealevelBarystaticOceanLoadEnum,&oceanaverage,P0Enum);
	#endif

	/*add ocean area into a global oceanareas vector:*/
	if(!loads->sealevelloads){
		oceanareas->SetValue(this->sid,oceanarea,INS_VAL);
		if(slgeom->issubelement[SLGEOM_OCEAN][this->lid]){
			int intj=slgeom->subelementmapping[SLGEOM_OCEAN][this->lid];
			subelementoceanareas->SetValue(intj,oceanarea,INS_VAL);
		}
	}
}
/*}}}*/
void       Tria::SealevelchangeConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* polarmotionvector,SealevelGeometry* slgeom){ /*{{{*/

	/*sal green function:*/
	int* AlphaIndex=NULL;
	int* AlphaIndexsub[SLGEOM_NUMLOADS];
	IssmDouble* G=NULL;
	IssmDouble* Grot=NULL;
	IssmDouble* rslfield=NULL;
	DoubleVecParam* parameter;
	bool computefuture=false;

	bool sal = false;
	bool viscous = false;
	bool rotation= false;
	bool percpu= false;
	int  size;
	int  nel,nbar;

	this->parameters->FindParam(&sal,SolidearthSettingsSelfAttractionEnum);
	this->parameters->FindParam(&viscous,SolidearthSettingsViscousEnum);
	this->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);
	this->parameters->FindParam(&nel,MeshNumberofelementsEnum);

	if(sal){
		parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeGViscoElasticEnum)); _assert_(parameter);
		parameter->GetParameterValueByPointer((IssmDouble**)&G,NULL);

		if (rotation)	this->inputs->GetArrayPtr(SealevelchangeGrotEnum,this->lid,&Grot,&size);

		rslfield=this->SealevelchangeGxL(G,Grot,loads,polarmotionvector,slgeom,nel,computefuture=false);
		this->SealevelchangeCollectGrdfield(sealevelpercpu,rslfield,slgeom,nel,percpu=true,SealevelchangeViscousRSLEnum,computefuture=false);

	}

	return;
} /*}}}*/
void       Tria::SealevelchangeDeformationConvolution(IssmDouble* sealevelpercpu, GrdLoads* loads, IssmDouble* polarmotionvector,SealevelGeometry* slgeom){ /*{{{*/

	IssmDouble SealevelGrd[3]={0,0,0};
	IssmDouble RSLGrd[3]={0,0,0};
	IssmDouble UGrd[3]={0,0,0};
	IssmDouble NGrd[3]={0,0,0};
	IssmDouble EGrd[3]={0,0,0};
	int nel,nbar;
	bool sal = false;
	int spatial_component=0;
	IssmDouble* G=NULL;
	IssmDouble* GU=NULL;
	IssmDouble* GH=NULL;
	IssmDouble* Grot=NULL;
	IssmDouble* GUrot=NULL;
	IssmDouble* GNrot=NULL;
	IssmDouble* GErot=NULL;
	IssmDouble* grdfield=NULL;

	DoubleVecParam* parameter;
	bool computefuture=false;

	int horiz;
	int size;

	bool rotation= false;
	bool elastic=false;
	bool percpu=false;
	bool planethasocean=false;

	this->parameters->FindParam(&nel,MeshNumberofelementsEnum);
	this->parameters->FindParam(&sal,SolidearthSettingsSelfAttractionEnum);
	this->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);
	this->parameters->FindParam(&elastic,SolidearthSettingsElasticEnum);
	this->parameters->FindParam(&horiz,SolidearthSettingsHorizEnum);
	this->parameters->FindParam(&planethasocean,SolidearthSettingsGrdOceanEnum);

	if(sal){
		parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeGViscoElasticEnum)); _assert_(parameter);
		parameter->GetParameterValueByPointer((IssmDouble**)&G,NULL);

		if(elastic){
			parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeUViscoElasticEnum)); _assert_(parameter);
			parameter->GetParameterValueByPointer((IssmDouble**)&GU,NULL);

			if(horiz){
				parameter = static_cast<DoubleVecParam*>(this->parameters->FindParamObject(SealevelchangeHViscoElasticEnum)); _assert_(parameter);
				parameter->GetParameterValueByPointer((IssmDouble**)&GH,NULL);
			}
			if (rotation) {
				this->inputs->GetArrayPtr(SealevelchangeGrotEnum,this->lid,&Grot,&size);
				this->inputs->GetArrayPtr(SealevelchangeGUrotEnum,this->lid,&GUrot,&size);
				if (horiz){
					this->inputs->GetArrayPtr(SealevelchangeGErotEnum,this->lid,&GErot,&size);
					this->inputs->GetArrayPtr(SealevelchangeGNrotEnum,this->lid,&GNrot,&size);
				}
			}
		}
		//Relative sea level convolution
		grdfield=this->SealevelchangeGxL(G,Grot,loads,polarmotionvector,slgeom,nel,computefuture=true);
		this->SealevelchangeCollectGrdfield(&RSLGrd[0],grdfield,slgeom,nel,percpu=false,SealevelchangeViscousRSLEnum,computefuture=true);

		if(elastic){
			//Bedrock Uplift
			grdfield=this->SealevelchangeGxL(GU,GUrot,loads,polarmotionvector,slgeom,nel,computefuture=true);
			this->SealevelchangeCollectGrdfield(&UGrd[0],grdfield,slgeom,nel,percpu=false,SealevelchangeViscousUEnum,computefuture=true);

			if(horiz){
				//Bedrock North displacement
				grdfield=this->SealevelchangeHorizGxL(spatial_component=1,GH,GNrot,loads,polarmotionvector,slgeom,nel,computefuture=true);
				this->SealevelchangeCollectGrdfield(&NGrd[0],grdfield,slgeom,nel,percpu=false,SealevelchangeViscousNEnum,computefuture=true);

				//Bedrock East displacement
				grdfield=this->SealevelchangeHorizGxL(spatial_component=2,GH,GErot,loads,polarmotionvector,slgeom,nel,computefuture=true);
				this->SealevelchangeCollectGrdfield(&EGrd[0],grdfield,slgeom,nel,percpu=false,SealevelchangeViscousEEnum,computefuture=true);
			}
		}
	}

	if (planethasocean){ //We must also output the RSL on vertices to compute the ocean mass conservation
		for(int i=0;i<NUMVERTICES;i++){
			if(slgeom->lids[this->vertices[i]->lid]==this->lid){
				sealevelpercpu[this->vertices[i]->lid]=RSLGrd[i];
			}
		}
	}

	/*Create geoid: */
	for(int i=0;i<NUMVERTICES;i++)SealevelGrd[i]=UGrd[i]+RSLGrd[i];

	/*Create inputs*/
	this->AddInput(SealevelGRDEnum,SealevelGrd,P1Enum);
	this->AddInput(BedGRDEnum,UGrd,P1Enum);
	if(horiz){
		this->AddInput(BedNorthGRDEnum,NGrd,P1Enum);
		this->AddInput(BedEastGRDEnum,EGrd,P1Enum);
	}

} /*}}}*/
IssmDouble*       Tria::SealevelchangeGxL(IssmDouble* G, IssmDouble* Grot, GrdLoads* loads, IssmDouble* polarmotionvector, SealevelGeometry* slgeom, int nel, bool computefuture) { /*{{{*/

	//This function performs the actual convolution between Green functions and surface Loads for a particular grd field
	int* AlphaIndex=NULL;
	int* AlphaIndexsub[SLGEOM_NUMLOADS];
	int* activevertices=NULL;
	IssmDouble* grdfield=NULL;
	int i,e,l,t,a, index, nbar, size, av,ae,b,c;
	bool rotation=false;
	int nt=1; //important, ensures there is a defined value if computeviscous is false
	int n_activevertices=0;

	//viscous
	bool computeviscous=false;
	int viscousindex=0; //important
	int viscousnumsteps=1; //important

	this->parameters->FindParam(&computeviscous,SolidearthSettingsViscousEnum);
	this->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);

	//Get green functions indexing & geometry
	this->inputs->GetIntArrayPtr(SealevelchangeConvolutionVerticesEnum,this->lid,&activevertices,&n_activevertices); //the order in which the vertices appear here should be the same as in slgeom->lids
	this->inputs->GetIntArrayPtr(SealevelchangeAlphaIndexEnum,this->lid,&AlphaIndex,&size);
	for (int l=0;l<SLGEOM_NUMLOADS;l++) this->inputs->GetIntArrayPtr(slgeom->AlphaIndexEnum(l),this->lid,&AlphaIndexsub[l],&size);

	if(computeviscous){
		this->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);
		this->parameters->FindParam(&viscousindex,SealevelchangeViscousIndexEnum);
		if(computefuture) {
			nt=viscousnumsteps-viscousindex+2; //number of time steps remaining to reach final_time, +1 is sufficient with no adaptative time stepping, +2 necessary otherwise; we assume the safe choice here for the sake of simplicity
			if (nt>viscousnumsteps) nt=viscousnumsteps;
		}
		else nt=1;
	}
	//allocate
	grdfield=xNewZeroInit<IssmDouble>(3*nt);

	//early return
	if (n_activevertices==0) return grdfield;

	if(rotation){ //add rotational feedback
		for(av=0;av<n_activevertices;av++) { //vertices
			i=activevertices[av];
			//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
			b=i*nt;
			for (int m=0;m<3;m++){ //polar motion components
				for(t=0;t<nt;t++){ //time
					int index=m*3*viscousnumsteps+i*viscousnumsteps+t;
					grdfield[b+t]+=Grot[index]*polarmotionvector[m];
				}
			}
		}
	}

	//Convolution
	for(av=0;av<n_activevertices;av++) { /*{{{*/
		i=activevertices[av];
		//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
		b=i*nt;
		c=av*nel;
		for(ae=0;ae<loads->nactiveloads;ae++){
			e=loads->combined_loads_index[ae];
			a=AlphaIndex[c+e]*viscousnumsteps;
			for(t=0;t<nt;t++){
				grdfield[b+t]+=G[a+t]*loads->combined_loads[ae];
			}
		}
		for(l=0;l<SLGEOM_NUMLOADS;l++){
			nbar=slgeom->nbar[l];
			c=av*nbar;
			for (ae=0;ae<loads->nactivesubloads[l];ae++){
				e=loads->combined_subloads_index[l][ae];
				a=AlphaIndexsub[l][c+e]*viscousnumsteps;
				for(t=0;t<nt;t++){
					grdfield[b+t]+=G[a+t]*loads->combined_subloads[l][ae];
				}
			}
		}
		//av+=1;
	} /*}}}*/

	return grdfield;

} /*}}}*/
IssmDouble*       Tria::SealevelchangeHorizGxL(int spatial_component, IssmDouble* G, IssmDouble* Grot, GrdLoads* loads, IssmDouble* polarmotionvector, SealevelGeometry* slgeom, int nel, bool computefuture) { /*{{{*/

	//This function performs the actual convolution between Green functions and surface Loads for a particular grd field
	int* AlphaIndex=NULL;
	int* AzimIndex=NULL;
	int* AlphaIndexsub[SLGEOM_NUMLOADS];
	int* AzimIndexsub[SLGEOM_NUMLOADS];
	int* activevertices = NULL;
	IssmDouble* grdfield=NULL;
	int i,e,l,t,a,b,c, index, nbar, av, ae,n_activevertices, size;
	bool rotation=false;
	IssmDouble* projected_loads=NULL;
	IssmDouble* projected_subloads[SLGEOM_NUMLOADS];
	IssmDouble* horiz_projection=NULL;
	IssmDouble* horiz_projectionsub[SLGEOM_NUMLOADS];
	int nt=1; //important, ensures there is a defined value if computeviscous is false

	//viscous
	bool computeviscous=false;
	int viscousindex=0; //important
	int viscousnumsteps=1; //important

	//Get green functions indexing & geometry
	this->inputs->GetIntArrayPtr(SealevelchangeConvolutionVerticesEnum,this->lid,&activevertices,&n_activevertices);
	this->inputs->GetIntArrayPtr(SealevelchangeAlphaIndexEnum,this->lid,&AlphaIndex,&size);
	for (int l=0;l<SLGEOM_NUMLOADS;l++) this->inputs->GetIntArrayPtr(slgeom->AlphaIndexEnum(l),this->lid,&AlphaIndexsub[l],&size);
	this->inputs->GetIntArrayPtr(SealevelchangeAzimuthIndexEnum,this->lid,&AzimIndex,&size);
	for (int l=0;l<SLGEOM_NUMLOADS;l++) this->inputs->GetIntArrayPtr(slgeom->AzimuthIndexEnum(l),this->lid,&AzimIndexsub[l],&size);

	//First, figure out how many time steps to compute grdfield for
	this->parameters->FindParam(&computeviscous,SolidearthSettingsViscousEnum);
	this->parameters->FindParam(&rotation,SolidearthSettingsRotationEnum);
	if(computeviscous){
		this->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);
		this->parameters->FindParam(&viscousindex,SealevelchangeViscousIndexEnum);
		if(computefuture) {
			nt=viscousnumsteps-viscousindex+2; //number of time steps remaining to reach final_time, +1 is sufficient with no adaptative time stepping, +2 necessary otherwise; we assume the safe choice here for the sake of simplicity
			if (nt>viscousnumsteps) nt=viscousnumsteps;
		}
		else nt=1;
	}
	//allocate
	grdfield=xNewZeroInit<IssmDouble>(3*nt);
	if (n_activevertices==0) return grdfield;

	if(rotation){ //add rotational feedback
		for(av=0;av<n_activevertices;av++) { //vertices
			i=activevertices[av];
			//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
			for (int m=0;m<3;m++){ //polar motion components
				for(t=0;t<nt;t++){ //time
					int index=m*3*viscousnumsteps+i*viscousnumsteps+t;
					grdfield[i*nt+t]+=Grot[index]*polarmotionvector[m];
				}
			}
			//}
		}
	}

	//Initialize projection vectors
	horiz_projection=xNewZeroInit<IssmDouble>(loads->nactiveloads);
	projected_loads=xNewZeroInit<IssmDouble>(loads->nactiveloads);
	for(l=0;l<SLGEOM_NUMLOADS;l++){
		//nbar=slgeom->nbar[l];
		projected_subloads[l]=xNewZeroInit<IssmDouble>(loads->nactivesubloads[l]);
		horiz_projectionsub[l]=xNewZeroInit<IssmDouble>(loads->nactivesubloads[l]);
	}

	//Convolution
	//av=0;
	for(av=0;av<n_activevertices;av++) { //vertices
		i=activevertices[av];
		//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
		b=i*nt;

		//GxL needs to be projected on the right axis before summation into the grd field
		//here we apply the projection scalar to the load prior to the actual convolution loop for more efficiency

		//get projection
		if (spatial_component==1){ //north
			for(ae=0;ae<loads->nactiveloads;ae++){
				e=loads->combined_loads_index[ae];
				horiz_projection[ae]=cos(2.0*M_PI*reCast<IssmDouble,int>(AzimIndex[av*nel+e])/65535.0); // 65535=2^16-1 = max value of 16 bits unsigned int
			}
			for(l=0;l<SLGEOM_NUMLOADS;l++){
				nbar=slgeom->nbar[l];
				for(ae=0;ae<loads->nactivesubloads[l];ae++){
					e=loads->combined_subloads_index[l][ae];
					horiz_projectionsub[l][ae]=cos(2.0*M_PI*reCast<IssmDouble,int>(AzimIndexsub[l][av*nbar+e])/65535.0);
				}
			}
		}
		else if (spatial_component==2){ //east
			for(ae=0;ae<loads->nactiveloads;ae++){
				e=loads->combined_loads_index[ae];
				horiz_projection[ae]=sin(2.0*M_PI*reCast<IssmDouble,int>(AzimIndex[av*nel+e])/65535.0);
			}
			for(l=0;l<SLGEOM_NUMLOADS;l++){
				nbar=slgeom->nbar[l];
				for(ae=0;ae<loads->nactivesubloads[l];ae++){
					e=loads->combined_subloads_index[l][ae];
					horiz_projectionsub[l][ae]=sin(2.0*M_PI*reCast<IssmDouble,int>(AzimIndexsub[l][av*nbar+e])/65535.0);
				}
			}
		}

		//project load in the right direction 
		for (ae=0;ae<loads->nactiveloads;ae++){
			projected_loads[ae]=loads->combined_loads[ae]*horiz_projection[ae];
		}
		for(l=0;l<SLGEOM_NUMLOADS;l++){
			nbar=slgeom->nbar[l];
			for(ae=0;ae<loads->nactivesubloads[l];ae++){
				projected_subloads[l][ae]=loads->combined_subloads[l][ae]*horiz_projectionsub[l][ae];
			}
		}

		//do the convolution
		c=av*nel;
		for(ae=0;ae<loads->nactiveloads;ae++){
			e=loads->combined_loads_index[ae];
			a=AlphaIndex[c+e]*viscousnumsteps;
			for(t=0;t<nt;t++){
				grdfield[b+t]+=G[a+t]*projected_loads[ae];
			}
		}
		for(l=0;l<SLGEOM_NUMLOADS;l++){
			nbar=slgeom->nbar[l];
			c=av*nbar;
			for(ae=0;ae<loads->nactivesubloads[l];ae++){
				e=loads->combined_subloads_index[l][ae];
				a=AlphaIndexsub[l][c+e]*viscousnumsteps;
				for(t=0;t<nt;t++){
					grdfield[b+t]+=G[a+t]*projected_subloads[l][ae];
				}
			}
		}
		//av+=1;
	} 

	//free resources
	xDelete<IssmDouble>(horiz_projection);
	xDelete<IssmDouble>(projected_loads);
	for(l=0;l<SLGEOM_NUMLOADS;l++) {
		xDelete<IssmDouble>(projected_subloads[l]);
		xDelete<IssmDouble>(horiz_projectionsub[l]);
	}
	return grdfield;

} /*}}}*/
void       Tria::SealevelchangeCollectGrdfield(IssmDouble* grdfieldout, IssmDouble* grdfield, SealevelGeometry* slgeom, int nel, bool percpu, int viscousenum, bool computefuture) { /*{{{*/

	//This function aligns grdfield with the requested output format: in a size 3 vector or in a size numberofvertices vector
	// if compute viscous is on, we also interpolate the field timewise given the current timestepping as well as collect viscous deformation and update the viscous deformation time series for future time steps
	int i,e,l,t,a, index, nbar, av, n_activevertices;
	int nt=1;

	//viscous
	bool computeviscous=false;
	int viscousindex=0; //important
	int viscousnumsteps=1; //important
	int* activevertices = NULL;
	IssmDouble* viscousfield=NULL;
	IssmDouble* grdfieldinterp=NULL;
	IssmDouble* viscoustimes=NULL;
	IssmDouble  final_time;
	IssmDouble  lincoeff;
	IssmDouble  timeacc;

	//parameters & initialization
	this->parameters->FindParam(&computeviscous,SolidearthSettingsViscousEnum);
	this->inputs->GetIntArrayPtr(SealevelchangeConvolutionVerticesEnum,this->lid,&activevertices,&n_activevertices);

	if(computeviscous){
		this->parameters->FindParam(&viscousnumsteps,SealevelchangeViscousNumStepsEnum);
		this->parameters->FindParam(&viscousindex,SealevelchangeViscousIndexEnum);
		this->parameters->FindParam(&viscoustimes,NULL,SealevelchangeViscousTimesEnum);
		this->parameters->FindParam(&final_time,TimesteppingFinalTimeEnum);
		this->parameters->FindParam(&timeacc,SolidearthSettingsTimeAccEnum);
		this->inputs->GetArrayPtr(viscousenum,this->lid,&viscousfield,NULL);
		if(computefuture) {
			nt=viscousnumsteps-viscousindex+2; //number of time steps remaining to reach final_time, +1 is sufficient with no adaptative time stepping, +2 necessary otherwise; we assume the safe choice here for the sake of simplicity
			if (nt>viscousnumsteps) nt=viscousnumsteps;
			grdfieldinterp=xNewZeroInit<IssmDouble>(3*viscousnumsteps); 
		}
		else nt=1;
	}

	if(!computeviscous){ /*{{{*/
		/*elastic or self attraction only case
		  store values computed on vertices, but don't repeat the computation if another element already computed it!:*/
		if(percpu){
			for(av=0;av<n_activevertices;av++) { //vertices
				i=activevertices[av];
				//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
				grdfieldout[this->vertices[i]->lid]=grdfield[i];
				//}
			}
		}
		else{
			for(i=0;i<NUMVERTICES;i++) grdfieldout[i]=grdfield[i];
		}
		//free resources
		xDelete<IssmDouble>(grdfield);
		return;
	}
	else { //viscous case
		// we need to do up to 3 things (* = only if computefuture)
		// 1: collect viscous grdfield from past loads due at present-day and add it to grdfield[current_time]
		// 2*: add new grdfield contribution to the viscous stack for future time steps
		// 3*: subtract from viscous stack the grdfield that has already been accounted for so we don't add it again at the next time step

		/*update grdfield at present time using viscous stack at present time: */
		for(av=0;av<n_activevertices;av++) { //vertices
			i=activevertices[av];
			//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
			grdfield[i*nt+0]+=viscousfield[i*viscousnumsteps+viscousindex]; 
		}

		/* Map new grdfield generated by present-day loads onto viscous time vector*/
		if(computefuture){
			//viscousindex time and first time step of grdfield coincide, so just copy that value
			for(av=0;av<n_activevertices;av++) { //vertices
				i=activevertices[av];
				//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
				grdfieldinterp[i*viscousnumsteps+viscousindex]=grdfield[i*nt+0];
			}
			if(viscoustimes[viscousindex]<final_time){
				//And interpolate the rest of the points in the future
				lincoeff=(viscoustimes[viscousindex+1]-viscoustimes[viscousindex])/timeacc;
				for(av=0;av<n_activevertices;av++) { //vertices
					i=activevertices[av];
					//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
					int i_time1= i*nt-viscousindex;
					int i_time2= i*viscousnumsteps;
					for(int t=viscousindex+1;t<viscousnumsteps;t++){
						grdfieldinterp[i_time2+t] = (1-lincoeff)*grdfield[i_time1+t-1]
									  +    lincoeff *grdfield[i_time1+t]
									  +          viscousfield[i_time2+t];
						/*update viscous stack with future deformation from present load: */
						viscousfield[i_time2+t]=grdfieldinterp[i_time2+t]
								       -grdfieldinterp[i_time2+viscousindex];
					}
				}
			}
			/*Save viscous stack now that we updated the values:*/
			this->inputs->SetArrayInput(viscousenum,this->lid,viscousfield,3*viscousnumsteps);
		}

		/*store values computed on vertices*/
		if(percpu){
			for(av=0;av<n_activevertices;av++) { //vertices
				i=activevertices[av];
				//if(slgeom->lids[this->vertices[i]->lid]!=this->lid)continue;
				grdfieldout[this->vertices[i]->lid]=grdfield[i*nt+0];
				//}
			}
		}
		else{
			for(i=0;i<NUMVERTICES;i++) grdfieldout[i]=grdfield[i*nt+0];
		}
		//free resources
		xDelete<IssmDouble>(grdfield);
		xDelete<IssmDouble>(viscoustimes);
		if (computefuture){
			xDelete<IssmDouble>(grdfieldinterp);
		}
		/*}}}*/
	}
} /*}}}*/

void       Tria::SealevelchangeShift(GrdLoads* loads,  IssmDouble offset, SealevelGeometry* slgeom){ /*{{{*/

	offset*=slgeom->LoadArea[SLGEOM_OCEAN][this->lid]; //kg.m^-2 to kg 
	if(slgeom->isoceanin[this->lid]){
		if(slgeom->issubelement[SLGEOM_OCEAN][this->lid]){
			int intj=slgeom->subelementmapping[SLGEOM_OCEAN][this->lid];
			loads->vsubsealevelloads->SetValue(intj,offset,ADD_VAL);
		}
		else loads->vsealevelloads->SetValue(this->sid,offset,ADD_VAL);
	}

} /*}}}*/
#endif

#ifdef _HAVE_DAKOTA_
void       Tria::InputScaleFromDakota(IssmDouble* distributed_values, IssmDouble* partition, int npart, int nt, int name){/*{{{*/

	int interp;
	int type;

	/*Branch according to whether we have a transient or not input: */
	type=this->inputs->GetInputObjectEnum(name);
	if(type==TriaInputEnum){
		/*Figure out if we are P0 or P1 interpolation: */
		TriaInput* triainput = this->inputs->GetTriaInput(name);
		TriaInput* triainput2 = this->inputs->GetTriaInput(DummyEnum);
		this->InputServe(triainput);
		interp=triainput->GetInterpolation();

		if (interp==P0Enum){
			/*Update the value if this element belongs to the partition: */
			if(partition[this->Sid()]!=-1){
				/*scale P0 value  for this element, corresponding to the partition:*/
				IssmDouble value = triainput->element_values[0];
				value*=distributed_values[(int)partition[this->Sid()]];
				triainput2->SetInput(P0Enum,this->lid,value);
			}
		}
		else if (interp==P1Enum){
			IssmDouble values[NUMVERTICES];
			int lidlist[NUMVERTICES];
			this->GetVerticesLidList(&lidlist[0]);
			for (int i=0;i<NUMVERTICES;i++){
				values[i]=triainput->element_values[i];
				if(partition[this->vertices[i]->Sid()]!=-1) values[i]*=distributed_values[(int)partition[this->vertices[i]->Sid()]];
			}
			triainput2->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
		}
		else _error_("Tria::InputScaleFromDakota error message: input interpolation " << EnumToStringx(interp) << " not supported yet!");
	}
	else if(type==TransientInputEnum){

		IssmDouble* steps=NULL;
		int nsteps;
		TransientInput* transientinput = NULL;
		TransientInput* transientinput2 = NULL;

		/*retrieve transient input:*/
		transientinput= this->inputs->GetTransientInput(name); _assert_(transientinput);
		transientinput2= this->inputs->GetTransientInput(DummyEnum); _assert_(transientinput2);

		/*retrieve time steps: */
		transientinput->GetAllTimes(&steps,&nsteps);

		/*double check:*/
		if (nsteps!=nt && nt!=1) _error_("Tria:InputScaleFromDakota error message: transient input " << EnumToStringx(name) <<  
				" should have the same number of time steps as the number of time values distributed by Dakota: " << nt << "\n");

		/*needed to update inputs:*/
		int lidlist[NUMVERTICES];
		this->GetVerticesLidList(&lidlist[0]);

		/*go through the transient inputs, and update:*/
		for (int i=0;i<nsteps;i++){
			TriaInput* triainput=transientinput->GetTriaInput(i);
			TriaInput* triainput2=transientinput2->GetTriaInput(i);
			this->InputServe(triainput);
			interp=triainput->GetInterpolation();

			if (interp==P0Enum){
				/*Update the value if this element belongs to the partition: */
				if(partition[this->Sid()]!=-1){
					/*scale P0 value  for this element, corresponding to the partition:*/
					IssmDouble value = triainput->element_values[0];
					if(nt==1) value*=distributed_values[(int)partition[this->Sid()]]; //we scale all the time steps  with the same distributed_value 
					else value*=distributed_values[(int)partition[this->Sid()]*nsteps+i]; //we scale all the time steps with distributed value for each step

					triainput2->SetInput(P0Enum,this->lid,value);
				}
			}
			else if (interp==P1Enum){
				IssmDouble values[NUMVERTICES];
				for (int j=0;j<NUMVERTICES;j++){
					values[j]=triainput->element_values[j];
					if(partition[this->vertices[j]->Sid()]!=-1){
						if(nt==1) values[j]*=distributed_values[(int)partition[this->vertices[j]->Sid()]];//we scale all the time steps  with the same distributed_value 
						else values[j]*=distributed_values[(int)partition[this->vertices[j]->Sid()]*nsteps+i];//we scale all the time steps with distributed value for each step
					}
				}
				triainput2->SetInput(P1Enum,NUMVERTICES,&lidlist[0],&values[0]);
			}
			else _error_("Tria::InputScaleFromDakota error message: input interpolation " << EnumToStringx(interp) << " not supported yet!");
		}
	}
	else _error_("Tria::InputScaleFromDakota error message: input type " << EnumToStringx(name) << " not supported yet!");
}
/*}}}*/
void       Tria::InputUpdateFromMatrixDakota(IssmDouble* matrix, int nrows, int ncols, int name, int type){/*{{{*/

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
				for(int i=0;i<3;i++){
					int row=this->vertices[i]->Sid();
					values[i]=matrix[ncols*row+t];
				}

				/*time:*/
				IssmDouble time=matrix[(nrows-1)*ncols+t];

				transientinput->AddTriaTimeInput(t,NUMVERTICES,&lidlist[0],&values[0],P1Enum);
			}
			break;

		case ElementEnum:
			/*Get value for the element: */
			for(int t=0;t<ncols;t++){ //ncols is the number of times
				IssmDouble value=matrix[ncols*(this->Sid())+t];
				IssmDouble time=matrix[(nrows-1)*ncols+t];
				transientinput->AddTriaTimeInput(t,1,&(this->lid),&value,P0Enum);
			}
			break;

		default:
			_error_("type " << type << " (" << EnumToStringx(type) << ") not implemented yet");
	}
}
/*}}}*/
void       Tria::InputUpdateFromVectorDakota(IssmDouble* vector, int name, int type){/*{{{*/

	int i,j;

	/*Check that name is an element input*/
	if(!IsInputEnum(name)) _error_("Enum "<<EnumToStringx(name)<<" is not in IsInput");

	switch(type){

		case VertexEnum:

			/*New TriaInput*/
			IssmDouble values[3];

			/*Get values on the 3 vertices*/
			for (i=0;i<3;i++){
				values[i]=vector[this->vertices[i]->Sid()]; //careful, vector of values here is not parallel distributed, but serial distributed (from a serial Dakota core!)
			}

			/*Branch on the specified type of update: */
			switch(name){
				case ThicknessEnum:
					IssmDouble  thickness[3];
					IssmDouble  thickness_init[3];
					IssmDouble  hydrostatic_ratio[3];
					IssmDouble  surface[3];
					IssmDouble  bed[3];

					/*retrieve inputs: */
					Element::GetInputListOnVertices(&thickness_init[0],ThicknessEnum);
					Element::GetInputListOnVertices(&hydrostatic_ratio[0],GeometryHydrostaticRatioEnum);
					Element::GetInputListOnVertices(&bed[0],BaseEnum);
					Element::GetInputListOnVertices(&surface[0],SurfaceEnum);

					/*build new bed and surface: */
					if (this->IsAllFloating()){
						/*hydrostatic equilibrium: */
						IssmDouble rho_ice,rho_water,di;
						rho_ice   = this->FindParam(MaterialsRhoIceEnum);
						rho_water = this->FindParam(MaterialsRhoSeawaterEnum);
						di        = rho_ice/rho_water;

						/*build new thickness: */
						for (j=0; j<3; j++) {
							/*  for observed/interpolated/hydrostatic thickness, remove scaling from any hydrostatic thickness  */
							if (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j]-(values[j]/thickness_init[j]-1.)*hydrostatic_ratio[j]*surface[j]/(1.-di);
							/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];

							/*  check the computed thickness and update bed*/
							if (thickness[j] < 0.) thickness[j]=1./(1.-di);
							bed[j]=surface[j]-thickness[j];
						}
					}
					else{
						/*build new thickness: */
						for (j=0; j<3; j++) {
							/*  for observed thickness, use scaled value  */
							if (hydrostatic_ratio[j] >= 0.)
								thickness[j]=values[j];
							/*  for minimum thickness, don't scale  */
							else
								thickness[j]=thickness_init[j];
						}

						/*update bed on grounded ice: */
						for(j=0;j<3;j++)bed[j]=surface[j]-thickness[j];
					}

					/*Add new inputs: */
					this->AddInput(ThicknessEnum,thickness,P1Enum);
					this->AddInput(BaseEnum,bed,P1Enum);
					this->AddInput(SurfaceEnum,surface,P1Enum);

					break;
				case MaterialsRheologyBEnum:
					this->AddInput(MaterialsRheologyBbarEnum,values,P1Enum);
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
