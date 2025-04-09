/*!\file Channel.c
 * \brief: implementation of the Channel object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <float.h> /*defines DBL_EPSILON*/
#include "shared/shared.h"
#include "../classes.h"
/*}}}*/	

/*Macros*/
#define NUMNODES    2
#define NUMVERTICES 2
#define C_W         4.22e3   /*specific heat capacity of water (J/kg/K)*/

/*Channel constructors and destructor*/
Channel::Channel(){/*{{{*/
	this->id         = -1;
	this->sid        = -1;
	this->parameters = NULL;
	this->helement   = NULL;
	this->element    = NULL;
	this->hnodes     = NULL;
	this->hvertices  = NULL;
	this->nodes      = NULL;
}
/*}}}*/
Channel::Channel(int channel_id,IssmDouble channelarea,int index,IoModel* iomodel){/*{{{*/
//Channel::Channel(int channel_id,int i,int index,IoModel* iomodel)

	this->id=channel_id;
	this->sid=channel_id-1;
	this->parameters = NULL;
	this->element    = NULL;
	this->nodes      = NULL;

	/*Set channel cross section to 0*/
	//this->S    = 0.;
	//this->Sold = 0.;
	this->S    = channelarea;
	this->Sold = channelarea;
	this->discharge = 0.;/*for output only*/

	/*Get edge info*/
	int i1 = iomodel->faces[4*index+0];
	int i2 = iomodel->faces[4*index+1];
	int e1 = iomodel->faces[4*index+2];
	int e2 = iomodel->faces[4*index+3];

	if(e2==-1){
		this->boundary = true;
	}
	else{
		this->boundary = false;
	}

	/*Set Element hook (4th column may be -1 for boundary edges)*/
	this->helement  = new Hook(&e1,1);

	/*Set Vertices hooks (4th column may be -1 for boundary edges)*/
	int channel_vertex_ids[2];
	channel_vertex_ids[0]=i1;
	channel_vertex_ids[1]=i2;
	this->hvertices =new Hook(&channel_vertex_ids[0],2);

	/*Set Nodes hooks (!! Assumes P1 CG)*/
	int channel_node_ids[2];
	channel_node_ids[0]=i1;
	channel_node_ids[1]=i2;
	this->hnodes=new Hook(&channel_node_ids[0],2);

}/*}}}*/
Channel::~Channel(){/*{{{*/
	this->parameters=NULL;
	delete helement;
	delete hnodes;
	delete hvertices;
}/*}}}*/

/*Object virtual functions definitions:*/
Object* Channel::copy() {/*{{{*/

	Channel* channel=NULL;

	channel=new Channel();

	/*copy fields: */
	channel->id=this->id;
	channel->S=this->S;

	/*point parameters: */
	channel->parameters=this->parameters;

	/*now deal with hooks and objects: */
	channel->hnodes    = (Hook*)this->hnodes->copy();
	channel->hvertices = (Hook*)this->hvertices->copy();
	channel->helement  = (Hook*)this->helement->copy();

	/*corresponding fields*/
	channel->nodes    = (Node**)channel->hnodes->deliverp();
	channel->vertices = (Vertex**)channel->hvertices->deliverp();
	channel->element  = (Element*)channel->helement->delivers();

	return channel;
}
/*}}}*/
void    Channel::DeepEcho(void){/*{{{*/

	_printf_("Channel:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   S:  " << S << "\n");
	hnodes->DeepEcho();
	hvertices->DeepEcho();
	helement->DeepEcho();
	_printf_("   parameters\n");
	if(parameters)
	 parameters->DeepEcho();
	else
	 _printf_("      NULL\n");
}		
/*}}}*/
void    Channel::Echo(void){/*{{{*/
	_printf_("Channel:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   S:  " << S << "\n");
	hnodes->Echo();
	hvertices->Echo();
	helement->Echo();
	_printf_("   parameters: " << parameters << "\n");
}
/*}}}*/
int     Channel::Id(void){/*{{{*/
	return id;
}
/*}}}*/
void    Channel::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	int object_enum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->id);
	marshallhandle->call(this->S);
	marshallhandle->call(this->Sold);
	marshallhandle->call(this->boundary);
	marshallhandle->call(this->discharge);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->hnodes      = new Hook();
		this->hvertices   = new Hook();
		this->helement    = new Hook();
	}

	this->hnodes->Marshall(marshallhandle);
	this->helement->Marshall(marshallhandle);
	this->hvertices->Marshall(marshallhandle);

	/*corresponding fields*/
	nodes    =(Node**)this->hnodes->deliverp();
	vertices =(Vertex**)this->hvertices->deliverp();
	element  =(Element*)this->helement->delivers();
}
/*}}}*/
int     Channel::ObjectEnum(void){/*{{{*/
	return ChannelEnum;
}/*}}}*/

/*Load virtual functions definitions:*/
void  Channel::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

	/*Take care of hooking up all objects for this element, ie links the objects in the hooks to their respective 
	 * datasets, using internal ids and offsets hidden in hooks: */
	hnodes->configure((DataSet*)nodesin);
	hvertices->configure((DataSet*)verticesin);
	helement->configure((DataSet*)elementsin);

	/*Initialize hooked fields*/
	this->nodes    = (Node**)hnodes->deliverp();
	this->vertices = (Vertex**)hvertices->deliverp();
	this->element  = (Element*)helement->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
}
/*}}}*/
void  Channel::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*recover some parameters*/
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case HydrologyGlaDSAnalysisEnum:
			Ke = this->CreateKMatrixHydrologyGlaDS();
			break;
		default:
			_error_("Don't know why we should be here");
	}

	/*Add to global matrix*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}

}
/*}}}*/
void  Channel::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	/*recover some parameters*/
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case HydrologyGlaDSAnalysisEnum:
			pe = this->CreatePVectorHydrologyGlaDS();
			break;
		default:
			_error_("Don't know why we should be here");
	}

	/*Add to global matrix*/
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}

}
/*}}}*/
void  Channel::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	for(int i=0;i<NUMNODES;i++) lidlist[i]=nodes[i]->Lid();
}
/*}}}*/
void  Channel::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	for(int i=0;i<NUMNODES;i++) sidlist[i]=nodes[i]->Sid();
}
/*}}}*/
int   Channel::GetNumberOfNodes(void){/*{{{*/
	return NUMNODES;
}
/*}}}*/
bool  Channel::IsPenalty(void){/*{{{*/
	return false;
}
/*}}}*/
void  Channel::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*No stiffness loads applied, do nothing: */
	return;

}
/*}}}*/
void  Channel::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*No penalty loads applied, do nothing: */
	return;

}
/*}}}*/
void  Channel::ResetHooks(){/*{{{*/

	this->nodes=NULL;
	this->vertices=NULL;
	this->element=NULL;
	this->parameters=NULL;

	/*Get Element type*/
	this->hnodes->reset();
	this->hvertices->reset();
	this->helement->reset();

}
/*}}}*/
void  Channel::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

}
/*}}}*/
void  Channel::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum){/*{{{*/

	/*Output */
	int d_nz = 0;
	int o_nz = 0;

	/*Loop over all nodes*/
	for(int i=0;i<this->GetNumberOfNodes();i++){

		if(!flags[this->nodes[i]->Lid()]){

			/*flag current node so that no other element processes it*/
			flags[this->nodes[i]->Lid()]=true;

			flagsindices[flagsindices_counter[0]]=this->nodes[i]->Lid();
         flagsindices_counter[0]++;

			/*if node is clone, we have an off-diagonal non-zero, else it is a diagonal non-zero*/
			switch(set2_enum){
				case FsetEnum:
					if(nodes[i]->FSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case GsetEnum:
					if(nodes[i]->gsize){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				case SsetEnum:
					if(nodes[i]->SSize()){
						if(this->nodes[i]->IsClone())
						 o_nz += 1;
						else
						 d_nz += 1;
					}
					break;
				default: _error_("not supported");
			}
		}
	}

	/*Assign output pointers: */
	*pd_nz=d_nz;
	*po_nz=o_nz;
}
/*}}}*/

/*Channel specific functions*/
ElementMatrix* Channel::CreateKMatrixHydrologyGlaDS(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceOnlyInElement()) return NULL;
	_assert_(tria->FiniteElement()==P1Enum); 
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);

	/*Intermediaries */
	IssmDouble  Jdet,v1,qc,fFactor,Afactor,Bfactor,Xifactor;
	IssmDouble  A,B,n,phi_old,phi,phi_0,dPw,ks,kc,Ngrad;
	IssmDouble  h_r;
	IssmDouble  H,h,b,dphi[2],dphids,dphimds,db[2],dbds;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_tria[3][3];
	const int   numnodes = NUMNODES;

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke=new ElementMatrix(this->nodes,NUMNODES,this->parameters);
	IssmDouble     basis[NUMNODES];
	IssmDouble     dbasisdx[2*NUMNODES];
	IssmDouble     dbasisds[NUMNODES];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0]     ,this->vertices,NUMVERTICES);
	GetVerticesCoordinates(&xyz_list_tria[0][0],tria->vertices,3);

	bool istransition;
	element->FindParam(&istransition,HydrologyIsTransitionEnum);
	bool isincludesheetthickness;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	IssmDouble L         = element->FindParam(MaterialsLatentheatEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	IssmDouble lc        = element->FindParam(HydrologyChannelSheetWidthEnum);
	IssmDouble c_t       = element->FindParam(HydrologyPressureMeltCoefficientEnum);
	IssmDouble alpha_c   = element->FindParam(HydrologyChannelAlphaEnum);
	IssmDouble beta_c    = element->FindParam(HydrologyChannelBetaEnum);
	IssmDouble alpha_s   = element->FindParam(HydrologySheetAlphaEnum);
	IssmDouble beta_s    = element->FindParam(HydrologySheetBetaEnum);
	IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);

	Input* h_input      = element->GetInput(HydrologySheetThicknessEnum);      _assert_(h_input);
	Input* H_input      = element->GetInput(ThicknessEnum);                    _assert_(H_input);
	Input* b_input      = element->GetInput(BedEnum);                          _assert_(b_input);
	Input* B_input      = element->GetInput(HydrologyRheologyBBaseEnum);       _assert_(B_input);
	Input* n_input      = element->GetInput(MaterialsRheologyNEnum);           _assert_(n_input);
	Input* ks_input     = element->GetInput(HydrologySheetConductivityEnum);   _assert_(ks_input);
	Input* kc_input     = element->GetInput(HydrologyChannelConductivityEnum); _assert_(kc_input);
	Input* hr_input     = element->GetInput(HydrologyBumpHeightEnum);          _assert_(hr_input);
	Input* phi_input    = element->GetInput(HydraulicPotentialEnum);           _assert_(phi_input);

	/*Get tangent vector*/
	IssmDouble tx = xyz_list_tria[index2][0] - xyz_list_tria[index1][0];
	IssmDouble ty = xyz_list_tria[index2][1] - xyz_list_tria[index1][1];
	IssmDouble Lt = sqrt(tx*tx+ty*ty);
	tx = tx/Lt;
	ty = ty/Lt;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());
		tria->GetSegmentNodalFunctionsDerivatives(&dbasisdx[0],&xyz_list_tria[0][0],gauss,index1,index2,tria->FiniteElement());
		dbasisds[0] = dbasisdx[0*2+0]*tx + dbasisdx[0*2+1]*ty;
		dbasisds[1] = dbasisdx[1*2+0]*tx + dbasisdx[1*2+1]*ty;

		/*Get input values at gauss points*/
		phi_input->GetInputDerivativeValue(&dphi[0],&xyz_list_tria[0][0],gauss);
		b_input->GetInputDerivativeValue(&db[0],&xyz_list_tria[0][0],gauss);
		phi_input->GetInputValue(&phi,gauss);
		h_input->GetInputValue(&h,gauss);
		ks_input->GetInputValue(&ks,gauss);
		kc_input->GetInputValue(&kc,gauss);
		hr_input->GetInputValue(&h_r,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);

		/*Get values for a few potentials*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		if(isincludesheetthickness) phi_0 += rho_water*g*h;
		dphids  = dphi[0]*tx + dphi[1]*ty;
		dphimds = rho_water*g*(db[0]*tx + db[1]*ty);
		Ngrad   = fabs(dphids);
		if(Ngrad<DBL_EPSILON) Ngrad = DBL_EPSILON;

		/*Compute the effective conductivity Kc = k h^alpha |grad Phi|^{beta-2} (same for sheet) and use transition model if specified*/
		IssmDouble Kc;
		IssmDouble Ks;
		IssmDouble nu = mu_water/rho_water;
		if(istransition==1 && omega>=DBL_EPSILON){
			IssmDouble hratio = h/h_r;
			IssmDouble coarg = 1. + 4.*omega*pow(hratio,3-2*alpha_s)*ks*pow(h,3)*Ngrad/nu;
			Ks = nu/2./omega*pow(hratio,2*alpha_s-3) * (-1 + pow(coarg, 0.5))/Ngrad;
			Kc = kc * pow(this->S,alpha_c) * pow(Ngrad,beta_c-2.);
		}
		else {
			Ks = ks*pow(h,alpha_s)*pow(Ngrad,beta_s-2.);
			Kc = kc * pow(this->S,alpha_c) * pow(Ngrad,beta_c-2.);
		}

		/*Approx. discharge in the sheet flowing folwing in the direction of the channel ofver a width lc*/
		qc = - Ks * dphids;

		/*d(phi - phi_m)/ds*/
		dPw = dphids - dphimds;

		/*Compute f factor*/
		fFactor = 0.;
		if(this->S>0. || qc*dPw>0.){
			fFactor = lc * qc;
		}

		/*Compute Afactor and Bfactor*/
		Afactor = C_W*c_t*rho_water;
		Bfactor = 1./L * (1./rho_ice - 1./rho_water);
		if(dphids>0){
			Xifactor = + Bfactor * (fabs(-Kc*dphids) + fabs(lc*qc));
		}
		else{
			Xifactor = - Bfactor * (fabs(-Kc*dphids) + fabs(lc*qc));
		}

		/*Diffusive term*/
		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				/*GlaDSCoupledSolver.F90 line 1659*/
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(
							+Kc*dbasisds[i]*dbasisds[j]                               /*Diffusion term*/
							- Afactor * Bfactor* Kc * dPw * basis[i] * dbasisds[j]    /*First part of Pi*/
							+ Afactor * fFactor * Bfactor * basis[i] * dbasisds[j]    /*Second part of Pi*/
							+ Xifactor* basis[i] * dbasisds[j]                        /*Xi term*/
							);
			}
		}

		/*Closing rate term*/ 
		/*See Gagliardini and Werder 2018 eq. A2 (v = v1*phi_i + v2(phi_{i+1}))*/
		A = pow(B,-n);
		if(phi_0-phi<0){
			v1 = 0.;
		}
		else{
			v1 = 2./pow(n,n)*A*S*(pow(fabs(phi_0-phi),n-1.)*( - n));

		}

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j] += gauss->weight*Jdet*(-v1)*basis[i]*basis[j];
			}
		}
	}

	/*Clean up and return*/
	delete gauss;
	return Ke;
}
/*}}}*/
ElementVector* Channel::CreatePVectorHydrologyGlaDS(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria* tria=(Tria*)element;
	if(!tria->IsIceOnlyInElement()) return NULL;
	_assert_(tria->FiniteElement()==P1Enum); 
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);

	/*Intermediaries */
	IssmDouble  Jdet,v2,Afactor,Bfactor,fFactor;
	IssmDouble  A,B,n,phi_old,phi,phi_0,dphimds,dphi[2];
	IssmDouble  H,h,b,db[2],dphids,qc,dPw,ks,kc,Ngrad;
	IssmDouble  h_r;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_tria[3][3];
	const int   numnodes = NUMNODES;

	/*Initialize Element vector and other vectors*/
	ElementVector* pe = new ElementVector(this->nodes,NUMNODES,this->parameters);
	IssmDouble     basis[NUMNODES];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);
	GetVerticesCoordinates(&xyz_list_tria[0][0],tria->vertices,3);

	bool istransition;
	element->FindParam(&istransition,HydrologyIsTransitionEnum);
	bool isincludesheetthickness;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	IssmDouble L         = element->FindParam(MaterialsLatentheatEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	IssmDouble lc        = element->FindParam(HydrologyChannelSheetWidthEnum);
	IssmDouble c_t       = element->FindParam(HydrologyPressureMeltCoefficientEnum);
	IssmDouble alpha_s   = element->FindParam(HydrologySheetAlphaEnum);
	IssmDouble beta_s    = element->FindParam(HydrologySheetBetaEnum);
	IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);

	Input* h_input      = element->GetInput(HydrologySheetThicknessEnum);      _assert_(h_input);
	Input* H_input      = element->GetInput(ThicknessEnum);                    _assert_(H_input);
	Input* b_input      = element->GetInput(BedEnum);                          _assert_(b_input);
	Input* B_input      = element->GetInput(HydrologyRheologyBBaseEnum);       _assert_(B_input);
	Input* n_input      = element->GetInput(MaterialsRheologyNEnum);           _assert_(n_input);
	Input* ks_input     = element->GetInput(HydrologySheetConductivityEnum);   _assert_(ks_input);
	Input* kc_input     = element->GetInput(HydrologyChannelConductivityEnum); _assert_(kc_input);
	Input* phi_input    = element->GetInput(HydraulicPotentialEnum);           _assert_(phi_input);
	Input* hr_input     = element->GetInput(HydrologyBumpHeightEnum);          _assert_(hr_input);

	/*Get tangent vector*/
	IssmDouble tx = xyz_list_tria[index2][0] - xyz_list_tria[index1][0];
	IssmDouble ty = xyz_list_tria[index2][1] - xyz_list_tria[index1][1];
	IssmDouble Lt = sqrt(tx*tx+ty*ty);
	tx = tx/Lt;
	ty = ty/Lt;

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());

		/*Get input values at gauss points*/
		b_input->GetInputDerivativeValue(&db[0],&xyz_list_tria[0][0],gauss);
		phi_input->GetInputDerivativeValue(&dphi[0],&xyz_list_tria[0][0],gauss);
		h_input->GetInputValue(&h,gauss);
		ks_input->GetInputValue(&ks,gauss);
		kc_input->GetInputValue(&kc,gauss);
		B_input->GetInputValue(&B,gauss);
		n_input->GetInputValue(&n,gauss);
		phi_input->GetInputValue(&phi,gauss);
		b_input->GetInputValue(&b,gauss);
		H_input->GetInputValue(&H,gauss);
		hr_input->GetInputValue(&h_r,gauss);

		/*Get values for a few potentials*/
		phi_0   = rho_water*g*b + rho_ice*g*H;
		if(isincludesheetthickness) phi_0 += rho_water*g*h;
		dphids  = dphi[0]*tx + dphi[1]*ty;
		dphimds = rho_water*g*(db[0]*tx + db[1]*ty);
		Ngrad   = fabs(dphids);
		if(Ngrad<DBL_EPSILON) Ngrad = DBL_EPSILON;

		/*Approx. discharge in the sheet flowing folwing in the direction of the channel ofver a width lc, use transition model if specified*/
		IssmDouble Ks;
		if (istransition==1 && omega>=DBL_EPSILON){
		IssmDouble hratio = h/h_r;
			IssmDouble nu = mu_water/rho_water;
			IssmDouble coarg = 1. + 4.*omega*pow(hratio,3-2*alpha_s)*ks*pow(h,3)*Ngrad/nu;
			Ks = nu/2./omega*pow(hratio,2*alpha_s-3) * (-1 + pow(coarg, 0.5))/Ngrad;
		}
		else {
			Ks = ks * pow(h,alpha_s) * pow(Ngrad,beta_s-2.);
		}

		/*Approx. discharge in the sheet flowing folwing in the direction of the channel ofver a width lc*/
		qc = - Ks * dphids;

		/*d(phi - phi_m)/ds*/
		dPw = dphids - dphimds;

		/*Compute f factor*/
		fFactor = 0.;
		if(this->S>0. || qc*dPw>0.){
			fFactor = lc * qc;
		}

		/*Compute Afactor and Bfactor*/
		Afactor = C_W*c_t*rho_water;
		Bfactor = 1./L * (1./rho_ice - 1./rho_water);

		/*Compute closing rate*/
		/*See Gagliardini and Werder 2018 eq. A2 (v = v2(phi_i) + v1*phi_{i+1})*/
		A = pow(B,-n);
		if(phi_0-phi<0){
			v2 = 0.;
		}
		else{
			v2 = 2./pow(n,n)*A*this->S*(pow(fabs(phi_0 - phi),n-1.)*(phi_0 +(n-1.)*phi));
		}

		for(int i=0;i<numnodes;i++){
			pe->values[i]+= - Jdet*gauss->weight*(-v2)*basis[i];
			pe->values[i]+= + Jdet*gauss->weight*Afactor*Bfactor*fFactor*dphimds*basis[i];
		}
	}

	/*Clean up and return*/
	delete gauss;
	return pe;
}
/*}}}*/
void           Channel::SetChannelCrossSectionOld(void){/*{{{*/

	this->Sold = this->S;

} /*}}}*/
void           Channel::UpdateChannelCrossSection(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(this->boundary || !tria->IsIceOnlyInElement()){
		this->S = 0.;
		return;
	}
	_assert_(tria->FiniteElement()==P1Enum); 

	/*Evaluate all fields on center of edge*/
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);
	GaussTria* gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);

	/*Set to 0 if inactive*/
	IssmDouble active;
	Input* active_input = element->GetInput(HydrologyMaskNodeActivationEnum); _assert_(active_input);
	active_input->GetInputValue(&active,gauss);
	if(active!=1.){
		this->S = 0.;
		delete gauss;
		return;
	}

	/*Intermediaries */
	IssmDouble  A,B,n,phi,phi_0,ks,kc,Ngrad;
	IssmDouble  h_r;
	IssmDouble  H,h,b,dphi[2],dphids,dphimds,db[2],dbds;
	IssmDouble  xyz_list[NUMVERTICES][3];
	IssmDouble  xyz_list_tria[3][3];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0]     ,this->vertices,NUMVERTICES);
	GetVerticesCoordinates(&xyz_list_tria[0][0],tria->vertices,3);

	bool istransition;
	element->FindParam(&istransition,HydrologyIsTransitionEnum);
	bool isincludesheetthickness;
	element->FindParam(&isincludesheetthickness,HydrologyIsIncludeSheetThicknessEnum);
	IssmDouble L         = element->FindParam(MaterialsLatentheatEnum);
	IssmDouble rho_ice   = element->FindParam(MaterialsRhoIceEnum);
	IssmDouble rho_water = element->FindParam(MaterialsRhoFreshwaterEnum);
	IssmDouble mu_water  = element->FindParam(MaterialsMuWaterEnum);
	IssmDouble g         = element->FindParam(ConstantsGEnum);
	IssmDouble lc        = element->FindParam(HydrologyChannelSheetWidthEnum);
	IssmDouble c_t       = element->FindParam(HydrologyPressureMeltCoefficientEnum);
	IssmDouble dt        = element->FindParam(TimesteppingTimeStepEnum);
	IssmDouble alpha_c   = element->FindParam(HydrologyChannelAlphaEnum);
	IssmDouble beta_c    = element->FindParam(HydrologyChannelBetaEnum);
	IssmDouble alpha_s   = element->FindParam(HydrologySheetAlphaEnum);
	IssmDouble beta_s    = element->FindParam(HydrologySheetBetaEnum);
	IssmDouble omega     = element->FindParam(HydrologyOmegaEnum);

	Input* h_input      = element->GetInput(HydrologySheetThicknessEnum);      _assert_(h_input);
	Input* H_input      = element->GetInput(ThicknessEnum);                    _assert_(H_input);
	Input* b_input      = element->GetInput(BedEnum);                          _assert_(b_input);
	Input* B_input      = element->GetInput(HydrologyRheologyBBaseEnum);       _assert_(B_input);
	Input* n_input      = element->GetInput(MaterialsRheologyNEnum);           _assert_(n_input);
	Input* ks_input     = element->GetInput(HydrologySheetConductivityEnum);   _assert_(ks_input);
	Input* kc_input     = element->GetInput(HydrologyChannelConductivityEnum); _assert_(kc_input);
	Input* phi_input    = element->GetInput(HydraulicPotentialEnum);           _assert_(phi_input);
	Input* hr_input     = element->GetInput(HydrologyBumpHeightEnum);          _assert_(hr_input);

	/*Get tangent vector*/
	IssmDouble tx = xyz_list_tria[index2][0] - xyz_list_tria[index1][0];
	IssmDouble ty = xyz_list_tria[index2][1] - xyz_list_tria[index1][1];
	IssmDouble Lt = sqrt(tx*tx+ty*ty);
	tx = tx/Lt;
	ty = ty/Lt;

	/*Get input values at gauss points*/
	phi_input->GetInputValue(&phi,gauss);
	phi_input->GetInputDerivativeValue(&dphi[0],&xyz_list_tria[0][0],gauss);
	h_input->GetInputValue(&h,gauss);
	ks_input->GetInputValue(&ks,gauss);
	kc_input->GetInputValue(&kc,gauss);
	B_input->GetInputValue(&B,gauss);
	n_input->GetInputValue(&n,gauss);
	b_input->GetInputValue(&b,gauss);
	b_input->GetInputDerivativeValue(&db[0],&xyz_list_tria[0][0],gauss);
	H_input->GetInputValue(&H,gauss);
	hr_input->GetInputValue(&h_r,gauss);

	/*Get values for a few potentials*/
	phi_0   = rho_water*g*b + rho_ice*g*H;
	if(isincludesheetthickness) phi_0 += rho_water*g*h;
	dphids  = dphi[0]*tx + dphi[1]*ty;
	dphimds = rho_water*g*(db[0]*tx + db[1]*ty);
	Ngrad   = fabs(dphids);
	if(Ngrad<DBL_EPSILON) Ngrad = DBL_EPSILON;

	/*d(phi - phi_m)/ds*/
	IssmDouble dPw = dphids - dphimds;

	/*Approx. discharge in the sheet flowing folwing in the direction of the channel ofver a width lc, use transition model if necessary*/
	IssmDouble qc;
	if (istransition==1 && omega>=DBL_EPSILON){
	IssmDouble hratio = h/h_r;
		IssmDouble nu = mu_water/rho_water;
		IssmDouble coarg = 1. + 4.*omega*pow(hratio,3-2*alpha_s)*ks*pow(h,3)*fabs(Ngrad)/nu;
		qc = -nu/2./omega*pow(hratio,2*alpha_s-3) * (-1 + pow(coarg, 0.5))*dphids/Ngrad;
	}
	else {
		qc = - ks * pow(h,alpha_s) * pow(Ngrad,beta_s-2.) * dphids;
	}

	/*Ice rate factor*/
	A = pow(B,-n);

	IssmDouble C = C_W*c_t*rho_water;
	IssmDouble Qprime = -kc * pow(Ngrad,beta_c-2.)*dphids;
	IssmDouble N = phi_0 - phi;

	bool converged  = false;
	int  count      = 0;

	while(!converged){

		IssmDouble Snew = this->S;

		/*Compute f factor*/
		IssmDouble fFactor = 0.;
		if(this->S>0. || qc*dPw>0.){
			fFactor = lc * qc;
		}

		IssmDouble alpha = 1./(rho_ice*L)*(
					fabs(Qprime*pow(Snew,alpha_c-1.)*dphids)
					+ C*Qprime*pow(Snew,alpha_c-1.)*dPw
					) - 2./pow(n,n)*A*pow(fabs(N),n-1.)*N;
		if(N<0){
			alpha = 1./(rho_ice*L)*(
               fabs(Qprime*pow(Snew,alpha_c-1.)*dphids)
               + C*Qprime*pow(Snew,alpha_c-1.)*dPw
               );
		}

		IssmDouble beta = 1./(rho_ice*L)*( fabs(lc*qc*dphids) + C*fFactor*dPw );

		/*Solve ODE*/
		this->S = ODE1(alpha,beta,this->Sold,dt,2);
		_assert_(!xIsNan<IssmDouble>(this->S)); 

		/*Constrain the cross section to be between 0 and 500 m^2*/
		if(this->S<0.)   this->S = 0.;
		if(this->S>500.) this->S = 500.;

		count++;

		if(fabs((this->S - Snew)/(Snew+DBL_EPSILON))<1e-8  || count>=10) converged = true;
	}

	/*Compute new channel discharge for output only*/
	IssmDouble Kc = kc * pow(this->S,alpha_c) * pow(Ngrad,beta_c-2.);
	this->discharge = -Kc*dphids;

	/*Clean up and return*/
	delete gauss;
}
/*}}}*/
void           Channel::WriteChannelCrossSection(IssmPDouble* values){/*{{{*/
	_assert_(values);
	values[this->sid] = reCast<IssmPDouble>(this->S);
}
/*}}}*/
void           Channel::WriteChannelDischarge(IssmPDouble* values){/*{{{*/
	_assert_(values);
	values[this->sid] = reCast<IssmPDouble>(this->discharge);
}
/*}}}*/
