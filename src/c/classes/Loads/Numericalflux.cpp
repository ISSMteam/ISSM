/*!\file Numericalflux.c
 * \brief: implementation of the Numericalflux object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "shared/shared.h"
#include "../classes.h"
/*}}}*/	

/*Load macros*/
#define NUMVERTICES 2

/*Numericalflux constructors and destructor*/
Numericalflux::Numericalflux(){/*{{{*/
	this->parameters = NULL;
	this->helement   = NULL;
	this->element    = NULL;
	this->hnodes     = NULL;
	this->hvertices  = NULL;
	this->nodes      = NULL;
}
/*}}}*/
Numericalflux::Numericalflux(int numericalflux_id,int i,int index,IoModel* iomodel){/*{{{*/

	/* Intermediary */
	int pos1,pos2,pos3,pos4;
	int numnodes;

	/*numericalflux constructor data: */
	int numericalflux_elem_ids[2];
	int numericalflux_vertex_ids[2];
	int numericalflux_node_ids[4];
	int numericalflux_type;
   int numericalflux_degree;

	/*Get edge*/
	int i1 = iomodel->faces[4*index+0];
	int i2 = iomodel->faces[4*index+1];
	int e1 = iomodel->faces[4*index+2];
	int e2 = iomodel->faces[4*index+3];

	/*First, see wether this is an internal or boundary edge (if e2=-1)*/
	if(e2==-1){
		/* Boundary edge, only one element */
		numericalflux_type=BoundaryEnum;
		numericalflux_elem_ids[0]=e1;
	}
	else{
		/* internal edge: connected to 2 elements */
		numericalflux_type=InternalEnum;
		numericalflux_elem_ids[0]=e1;
		numericalflux_elem_ids[1]=e2;
	}

   /*FIXME: hardcode element degree for now*/
   this->flux_degree= P1DGEnum;
   //this->flux_degree= P0DGEnum;

	/*1: Get vertices ids*/
	numericalflux_vertex_ids[0]=i1;
	numericalflux_vertex_ids[1]=i2;

	/*2: Get node ids*/
	if(numericalflux_type==InternalEnum){
		/*Get the column where these ids are located in the index*/
		pos1=pos2=pos3=pos4=UNDEF;
		for(int j=0;j<3;j++){
			if(iomodel->elements[3*(e1-1)+j]==i1) pos1=j+1;
			if(iomodel->elements[3*(e1-1)+j]==i2) pos2=j+1;
			if(iomodel->elements[3*(e2-1)+j]==i1) pos3=j+1;
			if(iomodel->elements[3*(e2-1)+j]==i2) pos4=j+1;
		}
		_assert_(pos1!=UNDEF && pos2!=UNDEF && pos3!=UNDEF && pos4!=UNDEF);

		/* We have the id of the elements and the position of the vertices in the index
		 * we can compute their dofs!*/
		numericalflux_node_ids[0]=3*(e1-1)+pos1;
		numericalflux_node_ids[1]=3*(e1-1)+pos2;
		numericalflux_node_ids[2]=3*(e2-1)+pos3;
		numericalflux_node_ids[3]=3*(e2-1)+pos4;
	}
	else{
		/*Get the column where these ids are located in the index*/
		pos1=pos2=UNDEF;
		for(int j=0;j<3;j++){
			if(iomodel->elements[3*(e1-1)+j]==i1) pos1=j+1;
			if(iomodel->elements[3*(e1-1)+j]==i2) pos2=j+1;
		}
		_assert_(pos1!=UNDEF && pos2!=UNDEF);

		/* We have the id of the elements and the position of the vertices in the index
		 * we can compute their dofs!*/
		numericalflux_node_ids[0]=3*(e1-1)+pos1;
		numericalflux_node_ids[1]=3*(e1-1)+pos2;
	}

   switch(this->flux_degree){
      case P0DGEnum:
         if(numericalflux_type==InternalEnum) numnodes = 2;
         else numnodes = 1;
			for(int i=0;i<numnodes;i++) numericalflux_node_ids[i] = numericalflux_elem_ids[i]; 
         numericalflux_node_ids[1] = numericalflux_elem_ids[1];
         break;
      case P1DGEnum:
         if(numericalflux_type==InternalEnum) numnodes = 4;
         else numnodes = 2;
         for(int i=0;i<numnodes;i++) numericalflux_node_ids[i] = numericalflux_node_ids[i]; //FIXME: to be improved...
         break;
      default:
         _error_("not supported yet");

   }

	/*Assign object fields: */
	this->id          = numericalflux_id;
	this->flux_type   = numericalflux_type;

	/*Hooks: */
	this->hnodes    = new Hook(numericalflux_node_ids,numnodes);
	this->hvertices = new Hook(&numericalflux_vertex_ids[0],2);
	this->helement  = new Hook(numericalflux_elem_ids,1); // take only the first element for now

	/*other fields*/
	this->parameters = NULL;
	this->element    = NULL;
	this->nodes      = NULL;
}
/*}}}*/
Numericalflux::~Numericalflux(){/*{{{*/
	this->parameters=NULL;
	delete helement;
	delete hnodes;
	delete hvertices;
}
/*}}}*/

/*Object virtual functions definitions:*/
Object* Numericalflux::copy() {/*{{{*/

	Numericalflux* numericalflux=NULL;

	numericalflux=new Numericalflux();

	/*copy fields: */
	numericalflux->id=this->id;
	numericalflux->flux_type=this->flux_type;
	numericalflux->flux_degree=this->flux_degree;

	/*point parameters: */
	numericalflux->parameters=this->parameters;

	/*now deal with hooks and objects: */
	numericalflux->hnodes    = (Hook*)this->hnodes->copy();
	numericalflux->hvertices = (Hook*)this->hvertices->copy();
	numericalflux->helement  = (Hook*)this->helement->copy();

	/*corresponding fields*/
	numericalflux->nodes    = (Node**)numericalflux->hnodes->deliverp();
	numericalflux->vertices = (Vertex**)numericalflux->hvertices->deliverp();
	numericalflux->element  = (Element*)numericalflux->helement->delivers();

	return numericalflux;
}
/*}}}*/
void    Numericalflux::DeepEcho(void){/*{{{*/

	_printf_("Numericalflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   flux_type: " << this->flux_type<< "\n");
	_printf_("   flux_degree: " << this->flux_degree<< "\n");
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
void    Numericalflux::Echo(void){/*{{{*/
	_printf_("Numericalflux:\n");
	_printf_("   id: " << id << "\n");
	_printf_("   flux_type: " << this->flux_type<< "\n");
	_printf_("   flux_degree: " << this->flux_degree<< "\n");
	hnodes->Echo();
	hvertices->Echo();
	helement->Echo();
	_printf_("   parameters: " << parameters << "\n");
}
/*}}}*/
int     Numericalflux::Id(void){/*{{{*/
	return id;
}
/*}}}*/
void    Numericalflux::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	_assert_(this);

	/*ok, marshall operations: */
	int object_enum = NumericalfluxEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->id);
	marshallhandle->call(this->flux_type);
	marshallhandle->call(this->flux_degree);

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
int     Numericalflux::ObjectEnum(void){/*{{{*/

	return NumericalfluxEnum;

}
/*}}}*/

/*Load virtual functions definitions:*/
void  Numericalflux::Configure(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/

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
void  Numericalflux::CreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs){/*{{{*/

	/*recover some parameters*/
	ElementMatrix* Ke=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*Just branch to the correct element stiffness matrix generator, according to the type of analysis we are carrying out: */
	switch(analysis_type){
		case MasstransportAnalysisEnum:
			Ke=CreateKMatrixMasstransport();
			break;
		case BalancethicknessAnalysisEnum:
			Ke=CreateKMatrixBalancethickness();
			break;
		case AdjointBalancethicknessAnalysisEnum:
			Ke=CreateKMatrixAdjointBalancethickness();
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global matrix*/
	if(Ke){
		Ke->AddToGlobal(Kff,Kfs);
		delete Ke;
	}

}
/*}}}*/
void  Numericalflux::CreatePVector(Vector<IssmDouble>* pf){/*{{{*/

	/*recover some parameters*/
	ElementVector* pe=NULL;
	int analysis_type;
	this->parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	switch(analysis_type){
		case MasstransportAnalysisEnum:
			pe=CreatePVectorMasstransport();
			break;
		case BalancethicknessAnalysisEnum:
			pe=CreatePVectorBalancethickness();
			break;
		case AdjointBalancethicknessAnalysisEnum:
			pe=CreatePVectorAdjointBalancethickness();
			break;
		default:
			_error_("analysis " << analysis_type << " (" << EnumToStringx(analysis_type) << ") not supported yet");
	}

	/*Add to global matrix*/
	if(pe){
		pe->AddToGlobal(pf);
		delete pe;
	}

}
/*}}}*/
void  Numericalflux::GetNodesLidList(int* lidlist){/*{{{*/

	_assert_(lidlist);
	_assert_(nodes);

	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++) lidlist[i]=nodes[i]->Lid();
}
/*}}}*/
void  Numericalflux::GetNodesSidList(int* sidlist){/*{{{*/

	_assert_(sidlist);
	_assert_(nodes);

	int numnodes = this->GetNumberOfNodes();
	for(int i=0;i<numnodes;i++) sidlist[i]=nodes[i]->Sid();
}
/*}}}*/
int   Numericalflux::GetNumberOfNodes(void){/*{{{*/

	if(this->flux_degree==P0DGEnum){
		switch(this->flux_type){
			case InternalEnum:
				return 2;
			case BoundaryEnum:
				return 1;
			default:
				_error_("Numericalflux type " << EnumToStringx(this->flux_type) << " not supported yet");
		}
	}
	else if(this->flux_degree==P1DGEnum){
		switch(this->flux_type){
			case InternalEnum:
				return 4;
			case BoundaryEnum:
				return 2;
			default:
				_error_("Numericalflux type " << EnumToStringx(this->flux_type) << " not supported yet");
		}
	}
	else{
		_error_("Numericalflux " << EnumToStringx(this->flux_degree) << " not supported yet");
	}

}
/*}}}*/
int   Numericalflux::GetNumberOfNodesOneSide(void){/*{{{*/

	if(this->flux_degree==P0DGEnum){
		return 1;
	}
	else if(this->flux_degree==P1DGEnum){
		return 2;
	}
	else{
		_error_("Numericalflux " << EnumToStringx(this->flux_degree) << " not supported yet");
	}

}
/*}}}*/
bool  Numericalflux::IsPenalty(void){/*{{{*/
	return false;
}
/*}}}*/
void  Numericalflux::PenaltyCreateKMatrix(Matrix<IssmDouble>* Kff, Matrix<IssmDouble>* Kfs,IssmDouble kmax){/*{{{*/

	/*No stiffness loads applied, do nothing: */
	return;

}
/*}}}*/
void  Numericalflux::PenaltyCreatePVector(Vector<IssmDouble>* pf,IssmDouble kmax){/*{{{*/

	/*No penalty loads applied, do nothing: */
	return;

}
/*}}}*/
void  Numericalflux::ResetHooks(){/*{{{*/

	this->nodes      = NULL;
	this->vertices   = NULL;
	this->element    = NULL;
	this->parameters = NULL;

	/*Get Element type*/
	this->hnodes->reset();
	this->hvertices->reset();
	this->helement->reset();

}
/*}}}*/
void  Numericalflux::SetCurrentConfiguration(Elements* elementsin,Loads* loadsin,Nodes* nodesin,Vertices* verticesin,Materials* materialsin,Parameters* parametersin){/*{{{*/
   /*Nothing to do :)*/

}
/*}}}*/
void  Numericalflux::SetwiseNodeConnectivity(int* pd_nz,int* po_nz,Node* node,bool* flags,int* flagsindices,int* flagsindices_counter,int set1_enum,int set2_enum){/*{{{*/

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

/*Numericalflux management*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixAdjointBalancethicknessInternal();
		case BoundaryEnum:
			return CreateKMatrixAdjointBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethicknessBoundary(void){/*{{{*/

	ElementMatrix* Ke=CreateKMatrixBalancethicknessBoundary();
	if(Ke) Ke->Transpose();
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixAdjointBalancethicknessInternal(void){/*{{{*/

	ElementMatrix* Ke=CreateKMatrixBalancethicknessInternal();
	if (Ke) Ke->Transpose();
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixBalancethicknessInternal();
		case BoundaryEnum:
			return CreateKMatrixBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethicknessBoundary(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria* tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble DL,Jdet,vx,vy,mean_vx,mean_vy;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input=tria->GetInput(VxEnum); _assert_(vxaverage_input); 
	Input* vyaverage_input=tria->GetInput(VyEnum); _assert_(vyaverage_input); 
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);

	GaussTria* gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;

	IssmDouble UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if(UdotN<=0){
		return NULL; /*(u,n)<0 -> inflow, PenaltyCreatePVector will take care of it*/
	}

	/*Initialize Element vector and other vectors*/
   int            numnodes = this->GetNumberOfNodes();
   ElementMatrix *Ke       = new ElementMatrix(nodes,numnodes,this->parameters);
   IssmDouble    *basis    = xNew<IssmDouble>(numnodes);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL=gauss->weight*Jdet*UdotN;

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j]+=DL*basis[i]*basis[j];
			}
		}
	} 

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixBalancethicknessInternal(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble A1,A2,Jdet,vx,vy,UdotN;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Fetch number of nodes for this flux*/
	int numnodes       = this->GetNumberOfNodes();
	int numnodes_plus  = this->GetNumberOfNodesOneSide();
	int numnodes_minus = numnodes_plus; /*For now we are not doing p-adaptive DG*/
	_assert_(numnodes==numnodes_plus+numnodes_minus);

	/*Initialize variables*/
	ElementMatrix *Ke = new ElementMatrix(nodes,numnodes,this->parameters);
	IssmDouble    *basis_plus  = xNew<IssmDouble>(numnodes_plus);
	IssmDouble    *basis_minus = xNew<IssmDouble>(numnodes_minus);

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input=tria->GetInput(VxEnum); _assert_(vxaverage_input); 
	Input* vyaverage_input=tria->GetInput(VyEnum); _assert_(vyaverage_input); 
	GetNormal(&normal[0],xyz_list);

	/* Start  looping on the number of gaussian points: */
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);
	GaussTria* gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis_plus[0] ,gauss,index1,index2,tria->FiniteElement());
		tria->GetSegmentNodalFunctions(&basis_minus[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		A1=gauss->weight*Jdet*UdotN/2;
		A2=gauss->weight*Jdet*fabs(UdotN)/2;

		/*Term 1 (numerical flux): {Hv}.[[phi]] = 0.5(H+v+ + H-v-)(phi+n+ + phi-n-)
		 *                                      = v.n/2 (H+phi+ + H-phi+ -H+phi- -H-phi-)
		 *                                      = v.n/2 (H+phi+ + H-phi+ -H+phi- -H-phi-)
		 *
		 *Term 2 (stabilization)  |v.n|/2 [[H]].[[phi]] = |v.n|/2 (H+n+ + H-n-)(phi+n+ + phi-n-)
		 *                                      = |v.n|/2 (H+phi+ -H-phi+ -H+phi- +H-phi-)
		 *     | A++ | A+- |
		 * K = |-----------|
		 *     | A-+ | A-- |
		 *
		 *These 4 terms for each expressions are added independently*/

		/*First term A++*/
		for(int i=0;i<numnodes_plus;i++){
			for(int j=0;j<numnodes_plus;j++){
				Ke->values[i*numnodes+j] += A1*(basis_plus[j]*basis_plus[i]);
				Ke->values[i*numnodes+j] += A2*(basis_plus[j]*basis_plus[i]);
			}
		}
		/*Second term A+-*/
		for(int i=0;i<numnodes_plus;i++){
			for(int j=0;j<numnodes_minus;j++){
				Ke->values[i*numnodes+numnodes_plus+j] +=  A1*(basis_minus[j]*basis_plus[i]);
				Ke->values[i*numnodes+numnodes_plus+j] += -A2*(basis_minus[j]*basis_plus[i]);
			}
		}
		/*Third term A-+*/
		for(int i=0;i<numnodes_minus;i++){
			for(int j=0;j<numnodes_plus;j++){
				Ke->values[(numnodes_plus+i)*numnodes+j] += -A1*(basis_plus[j]*basis_minus[i]);
				Ke->values[(numnodes_plus+i)*numnodes+j] += -A2*(basis_plus[j]*basis_minus[i]);
			}
		}
		/*Fourth term A-+*/
		for(int i=0;i<numnodes_minus;i++){
			for(int j=0;j<numnodes_minus;j++){
				Ke->values[(numnodes_plus+i)*numnodes+numnodes_plus+j] += -A1*(basis_minus[j]*basis_minus[i]);
				Ke->values[(numnodes_plus+i)*numnodes+numnodes_plus+j] +=  A2*(basis_minus[j]*basis_minus[i]);
			}
		}
	}

	/*Clean up and return*/
   xDelete<IssmDouble>(basis_plus);
   xDelete<IssmDouble>(basis_minus);
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransport(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreateKMatrixMasstransportInternal();
		case BoundaryEnum:
			return CreateKMatrixMasstransportBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransportBoundary(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble DL,Jdet,vx,vy,mean_vx,mean_vy;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	IssmDouble dt = parameters->FindParam(TimesteppingTimeStepEnum);
	Input* vxaverage_input=tria->GetInput(VxEnum); _assert_(vxaverage_input);
	Input* vyaverage_input=tria->GetInput(VyEnum); _assert_(vyaverage_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);

	GaussTria* gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;

	IssmDouble UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if(UdotN<=0){
		return NULL; /*(u,n)<0 -> inflow, PenaltyCreatePVector will take care of it*/
	}

	/*Initialize Element vector and other vectors*/
   int            numnodes = this->GetNumberOfNodes();
   ElementMatrix *Ke       = new ElementMatrix(nodes,numnodes,this->parameters);
   IssmDouble    *basis    = xNew<IssmDouble>(numnodes);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL=gauss->weight*Jdet*dt*UdotN;

		for(int i=0;i<numnodes;i++){
			for(int j=0;j<numnodes;j++){
				Ke->values[i*numnodes+j]+=DL*basis[i]*basis[j];
			}
		}
	} 

	/*Clean up and return*/
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Ke;
}
/*}}}*/
ElementMatrix* Numericalflux::CreateKMatrixMasstransportInternal(void){/*{{{*/

	/*Initialize Element matrix and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble A1,A2,Jdet,vx,vy,UdotN;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Fetch number of nodes for this flux*/
	int numnodes       = this->GetNumberOfNodes();
	int numnodes_plus  = this->GetNumberOfNodesOneSide();
	int numnodes_minus = numnodes_plus; /*For now we are not doing p-adaptive DG*/
	_assert_(numnodes==numnodes_plus+numnodes_minus);

	/*Initialize variables*/
	ElementMatrix *Ke = new ElementMatrix(nodes,numnodes,this->parameters);
	IssmDouble    *basis_plus  = xNew<IssmDouble>(numnodes_plus);
	IssmDouble    *basis_minus = xNew<IssmDouble>(numnodes_minus);

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	IssmDouble dt = parameters->FindParam(TimesteppingTimeStepEnum);
	Input* vxaverage_input=tria->GetInput(VxEnum); _assert_(vxaverage_input); 
	Input* vyaverage_input=tria->GetInput(VyEnum); _assert_(vyaverage_input); 
	GetNormal(&normal[0],xyz_list);

	/* Start  looping on the number of gaussian points: */
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);
	GaussTria* gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis_plus[0] ,gauss,index1,index2,tria->FiniteElement());
		tria->GetSegmentNodalFunctions(&basis_minus[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		A1=gauss->weight*Jdet*dt*UdotN/2;
		A2=gauss->weight*Jdet*dt*fabs(UdotN)/2;

		/*Term 1 (numerical flux): {Hv}.[[phi]] = 0.5(H+v+ + H-v-)(phi+n+ + phi-n-)
		 *                                      = v.n/2 (H+phi+ + H-phi+ -H+phi- -H-phi-)
		 *                                      = v.n/2 (H+phi+ + H-phi+ -H+phi- -H-phi-)
		 *
		 *Term 2 (stabilization)  |v.n|/2 [[H]].[[phi]] = |v.n|/2 (H+n+ + H-n-)(phi+n+ + phi-n-)
		 *                                      = |v.n|/2 (H+phi+ -H-phi+ -H+phi- +H-phi-)
		 *     | A++ | A+- |
		 * K = |-----------|
		 *     | A-+ | A-- |
		 *
		 *These 4 terms for each expressions are added independently*/

		/*First term A++*/
		for(int i=0;i<numnodes_plus;i++){
			for(int j=0;j<numnodes_plus;j++){
				Ke->values[i*numnodes+j] += A1*(basis_plus[j]*basis_plus[i]);
				Ke->values[i*numnodes+j] += A2*(basis_plus[j]*basis_plus[i]);
			}
		}
		/*Second term A+-*/
		for(int i=0;i<numnodes_plus;i++){
			for(int j=0;j<numnodes_minus;j++){
				Ke->values[i*numnodes+numnodes_plus+j] +=  A1*(basis_minus[j]*basis_plus[i]);
				Ke->values[i*numnodes+numnodes_plus+j] += -A2*(basis_minus[j]*basis_plus[i]);
			}
		}
		/*Third term A-+*/
		for(int i=0;i<numnodes_minus;i++){
			for(int j=0;j<numnodes_plus;j++){
				Ke->values[(numnodes_plus+i)*numnodes+j] += -A1*(basis_plus[j]*basis_minus[i]);
				Ke->values[(numnodes_plus+i)*numnodes+j] += -A2*(basis_plus[j]*basis_minus[i]);
			}
		}
		/*Fourth term A-+*/
		for(int i=0;i<numnodes_minus;i++){
			for(int j=0;j<numnodes_minus;j++){
				Ke->values[(numnodes_plus+i)*numnodes+numnodes_plus+j] += -A1*(basis_minus[j]*basis_minus[i]);
				Ke->values[(numnodes_plus+i)*numnodes+numnodes_plus+j] +=  A2*(basis_minus[j]*basis_minus[i]);
			}
		}
	}

	/*Clean up and return*/
   xDelete<IssmDouble>(basis_plus);
   xDelete<IssmDouble>(basis_minus);
	delete gauss;
	return Ke;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorAdjointBalancethickness(void){/*{{{*/

	/*No PVector for the Adjoint*/
	return NULL;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethickness(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreatePVectorBalancethicknessInternal();
		case BoundaryEnum:
			return CreatePVectorBalancethicknessBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethicknessBoundary(void){/*{{{*/

	/*Initialize Load Vector and return if necessary*/
	Tria*  tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble DL,Jdet,vx,vy,mean_vx,mean_vy,thickness;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	Input* vxaverage_input = tria->GetInput(VxEnum);        _assert_(vxaverage_input);
	Input* vyaverage_input = tria->GetInput(VyEnum);        _assert_(vyaverage_input);
	Input* thickness_input = tria->GetInput(ThicknessEnum); _assert_(thickness_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);
	GaussTria* gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;
	IssmDouble UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if(UdotN>0){
		return NULL; /*(u,n)>0 -> outflow, PenaltyCreateKMatrix will take care of it*/
	}

	/*Initialize Load Vector */
	int            numnodes = this->GetNumberOfNodes();
	ElementVector *pe       = new ElementVector(nodes,numnodes,this->parameters);
	IssmDouble    *basis    = xNew<IssmDouble>(numnodes);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		thickness_input->GetInputValue(&thickness,gauss);

		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL= - gauss->weight*Jdet*UdotN*thickness;

		for(int i=0;i<numnodes;i++) pe->values[i] += DL*basis[i];
	}

	/*Clean up and return*/
   xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorBalancethicknessInternal(void){/*{{{*/

	/*Nothing added to PVector*/
	return NULL;

}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransport(void){/*{{{*/

	switch(this->flux_type){
		case InternalEnum:
			return CreatePVectorMasstransportInternal();
		case BoundaryEnum:
			return CreatePVectorMasstransportBoundary();
		default:
			_error_("type not supported yet");
	}
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransportBoundary(void){/*{{{*/

	/*Initialize Load Vector and return if necessary*/
	Tria* tria=(Tria*)element;
	if(!tria->IsIceInElement()) return NULL;

	/* Intermediaries*/
	IssmDouble DL,Jdet,vx,vy,mean_vx,mean_vy,thickness;
	IssmDouble xyz_list[NUMVERTICES][3];
	IssmDouble normal[2];

	/*Retrieve all inputs and parameters*/
	GetVerticesCoordinates(&xyz_list[0][0],vertices,NUMVERTICES);
	IssmDouble dt = parameters->FindParam(TimesteppingTimeStepEnum);
	Input* vxaverage_input    = tria->GetInput(VxEnum);                        _assert_(vxaverage_input);
	Input* vyaverage_input    = tria->GetInput(VyEnum);                        _assert_(vyaverage_input);
	Input* spcthickness_input = tria->GetInput(MasstransportSpcthicknessEnum); _assert_(spcthickness_input);
	GetNormal(&normal[0],xyz_list);

	/*Check wether it is an inflow or outflow BC (0 is the middle of the segment)*/
	int index1=tria->GetVertexIndex(vertices[0]);
	int index2=tria->GetVertexIndex(vertices[1]);
	GaussTria* gauss=new GaussTria();
	gauss->GaussEdgeCenter(index1,index2);
	vxaverage_input->GetInputValue(&mean_vx,gauss);
	vyaverage_input->GetInputValue(&mean_vy,gauss);
	delete gauss;
	IssmDouble UdotN=mean_vx*normal[0]+mean_vy*normal[1];
	if(UdotN>0){
		return NULL; /*(u,n)>0 -> outflow, PenaltyCreateKMatrix will take care of it*/
	}

	/*Initialize Load Vector */
	int            numnodes = this->GetNumberOfNodes();
	ElementVector *pe       = new ElementVector(nodes,numnodes,this->parameters);
	IssmDouble    *basis    = xNew<IssmDouble>(numnodes);

	/* Start  looping on the number of gaussian points: */
	gauss=new GaussTria(index1,index2,2);
	while(gauss->next()){

		tria->GetSegmentNodalFunctions(&basis[0],gauss,index1,index2,tria->FiniteElement());

		vxaverage_input->GetInputValue(&vx,gauss);
		vyaverage_input->GetInputValue(&vy,gauss);
		spcthickness_input->GetInputValue(&thickness,gauss);
		if(xIsNan<IssmDouble>(thickness)) _error_("Cannot weakly apply constraint because NaN was provided");

		UdotN=vx*normal[0]+vy*normal[1];
		tria->GetSegmentJacobianDeterminant(&Jdet,&xyz_list[0][0],gauss);
		DL= - gauss->weight*Jdet*dt*UdotN*thickness;

		for(int i=0;i<numnodes;i++) pe->values[i] += DL*basis[i];
	}

	/*Clean up and return*/
   xDelete<IssmDouble>(basis);
	delete gauss;
	return pe;
}
/*}}}*/
ElementVector* Numericalflux::CreatePVectorMasstransportInternal(void){/*{{{*/

	/*Nothing added to PVector*/
	return NULL;

}
/*}}}*/
void           Numericalflux::GetNormal(IssmDouble* normal,IssmDouble xyz_list[4][3]){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble vector[2];

	vector[0]=xyz_list[1][0] - xyz_list[0][0];
	vector[1]=xyz_list[1][1] - xyz_list[0][1];

	IssmDouble norm=sqrt(pow(vector[0],2.0)+pow(vector[1],2.0));

	normal[0]= + vector[1]/norm;
	normal[1]= - vector[0]/norm;
}
/*}}}*/
