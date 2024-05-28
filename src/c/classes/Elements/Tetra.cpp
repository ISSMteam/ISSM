/*!\file Tetra.cpp
 * \brief: implementation of the Tetrament object
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
#include "../classes.h"
#include "../Inputs/ElementInput.h"
#include "../../shared/shared.h"
/*}}}*/

/*Element macros*/
#define NUMVERTICES 4

/*Constructors/destructor/copy*/
Tetra::Tetra(int tet_id, int tet_sid,int tet_lid,IoModel* iomodel,int nummodels)/*{{{*/
		:ElementHook(nummodels,tet_id,NUMVERTICES,iomodel){

			/*id: */
			this->id  = tet_id;
			this->sid = tet_sid;
			this->lid = tet_lid;

			/*surface and base*/
			this->isonsurface = false;
			this->isonbase    = false;

			//this->parameters: we still can't point to it, it may not even exist. Configure will handle this.
			this->parameters = NULL;

			/*initialize pointers:*/
			this->nodes    = NULL;
			this->vertices = NULL;
			this->material = NULL;

			/*Only allocate pointer*/
			this->element_type_list=xNew<int>(nummodels);

			/*surface and base*/
			_assert_(iomodel->Data("md.mesh.vertexonsurface"));
			_assert_(iomodel->Data("md.mesh.vertexonbase"));
			this->isonsurface = false;
			this->isonbase    = false;
			IssmDouble sum = 0.;
			for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonsurface")[reCast<int>(iomodel->elements[(tet_id-1)*NUMVERTICES+i])-1];
			_assert_(sum>=0 && sum<4);
			if(sum>2.5) this->isonsurface = true;
			sum = 0.;
			for(int i=0;i<NUMVERTICES;i++) sum += iomodel->Data("md.mesh.vertexonbase")[reCast<int>(iomodel->elements[(tet_id-1)*NUMVERTICES+i])-1];
			_assert_(sum>=0 && sum<4);
			if(sum>2.5) this->isonbase = true;
		}
/*}}}*/
Tetra::~Tetra(){/*{{{*/
	this->parameters=NULL;
}
/*}}}*/
Object* Tetra::copy() {/*{{{*/

	int i;
	Tetra* tetra=NULL;

	tetra=new Tetra();

	//deal with TetraRef mother class
	int nanalyses = this->numanalyses;
	if(nanalyses > 0){
		tetra->element_type_list=xNew<int>(nanalyses);
		for(i=0;i<nanalyses;i++){
			if (this->element_type_list[i]) tetra->element_type_list[i]=this->element_type_list[i];
			else tetra->element_type_list[i] = 0;
		}
	}
	else tetra->element_type_list = NULL;
	tetra->element_type=this->element_type;
	tetra->numanalyses=nanalyses;

	//deal with ElementHook mother class
	if (this->hnodes){
		tetra->hnodes=xNew<Hook*>(tetra->numanalyses);
		for(i=0;i<tetra->numanalyses;i++){
			if (this->hnodes[i]) tetra->hnodes[i] = (Hook*)(this->hnodes[i]->copy());
			else tetra->hnodes[i] = NULL;
		}
	}
	else tetra->hnodes = NULL;

	tetra->hvertices = (Hook*)this->hvertices->copy();
	tetra->hmaterial = (Hook*)this->hmaterial->copy();
	tetra->hneighbors = NULL;

	/*deal with Tria fields: */
	tetra->id  = this->id;
	tetra->sid = this->sid;
	tetra->lid = this->lid;
	tetra->isonbase  = this->isonbase;
	tetra->isonsurface  = this->isonsurface;

	/*point parameters: */
	tetra->parameters=this->parameters;

	/*recover objects: */
	unsigned int num_nodes = 4;
	tetra->nodes = xNew<Node*>(num_nodes); //we cannot rely on an analysis_counter to tell us which analysis_type we are running, so we just copy the nodes.
	for(i=0;i<num_nodes;i++) if(this->nodes[i]) tetra->nodes[i]=this->nodes[i]; else tetra->nodes[i] = NULL;

	tetra->vertices = (Vertex**)this->hvertices->deliverp();
	tetra->material = (Material*)this->hmaterial->delivers();

	return tetra;
}
/*}}}*/
void Tetra::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = TetraEnum;
   marshallhandle->call(object_enum);

	marshallhandle->call(this->isonsurface);
	marshallhandle->call(this->isonbase);

	/*Call parent classes: */
	ElementHook::Marshall(marshallhandle);
	Element::MarshallElement2(marshallhandle,this->numanalyses);
	TetraRef::Marshall(marshallhandle);

	vertices = (Vertex**)this->hvertices->deliverp();
	material = (Material*)this->hmaterial->delivers();

}
/*}}}*/

void     Tetra::Configure(Elements* elementsin, Loads* loadsin, Nodes* nodesin,Vertices* verticesin, Materials* materialsin, Parameters* parametersin,Inputs* inputsin){/*{{{*/

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

	/*Now, go pick up the objects inside the hooks: */
	if (this->hnodes[analysis_counter]) this->nodes=(Node**)this->hnodes[analysis_counter]->deliverp();
	else this->nodes=NULL;
	this->vertices          = (Vertex**)this->hvertices->deliverp();
	this->material          = (Material*)this->hmaterial->delivers();

	/*point parameters to real dataset: */
	this->parameters=parametersin;
	this->inputs=inputsin;
}
/*}}}*/
void     Tetra::ElementSizes(IssmDouble* hx,IssmDouble* hy,IssmDouble* hz){/*{{{*/

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
void     Tetra::FaceOnBaseIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1. && values[indices[i][2]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
void     Tetra::FaceOnFrontIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MaskIceLevelsetEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 0. && values[indices[i][1]] == 0. && values[indices[i][2]] == 0.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
void     Tetra::FaceOnSurfaceIndices(int* pindex1,int* pindex2,int* pindex3){/*{{{*/

	IssmDouble values[NUMVERTICES];
	int        indices[4][3] = {{0,1,2},{0,3,1},{1,3,2},{0,2,3}};

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);

	for(int i=0;i<4;i++){
		if(values[indices[i][0]] == 1. && values[indices[i][1]] == 1. && values[indices[i][2]] == 1.){
			*pindex1 = indices[i][0];
			*pindex2 = indices[i][1];
			*pindex3 = indices[i][2];
			return;
		}
	}

	_error_("Could not find 3 vertices on bed");
}
/*}}}*/
int      Tetra::FiniteElement(void){/*{{{*/
	return this->element_type;
} /*}}}*/
int      Tetra::GetElementType(){/*{{{*/

	/*return TetraRef field*/
	return this->element_type;
}
/*}}}*/
Input*    Tetra::GetInput(int inputenum){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
Input*    Tetra::GetInput(int inputenum,IssmDouble time){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
void       Tetra::GetInputListOnVertices(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/* Start looping on the number of vertices: */
	if(input){
		GaussTetra gauss;
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
void       Tetra::GetInputListOnNodes(IssmDouble* pvalue,Input* input,IssmDouble default_value){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(pvalue);

	/*What type of finite element are we dealing with?*/
	int fe       = this->FiniteElement();
	int numnodes = this->GetNumberOfNodes();

	/* Start looping on the number of vertices: */
	if(input){
		GaussTetra gauss;
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
void     Tetra::GetInputValue(IssmDouble* pvalue,Node* node,int enumtype){/*{{{*/

	Input* input=this->GetInput(enumtype);
	if(!input) _error_("No input of type " << EnumToStringx(enumtype) << " found in tria");

	GaussTetra* gauss=new GaussTetra();
	gauss->GaussVertex(this->GetNodeIndex(node));

	input->GetInputValue(pvalue,gauss);
	delete gauss;
}
/*}}}*/
int      Tetra::GetNumberOfNodes(void){/*{{{*/
	return this->NumberofNodes(this->element_type);
}
/*}}}*/
int      Tetra::GetNumberOfVertices(void){/*{{{*/
	return NUMVERTICES;
}
/*}}}*/
void     Tetra::GetVerticesCoordinatesBase(IssmDouble** pxyz_list){/*{{{*/

	int        indices[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(3*3);
	this->FaceOnBaseIndices(&indices[0],&indices[1],&indices[2]);
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
void     Tetra::GetVerticesCoordinatesTop(IssmDouble** pxyz_list){/*{{{*/

	int        indices[3];
	IssmDouble xyz_list[NUMVERTICES][3];

	/*Element XYZ list*/
	::GetVerticesCoordinates(&xyz_list[0][0],this->vertices,NUMVERTICES);

	/*Allocate Output*/
	IssmDouble* xyz_list_edge = xNew<IssmDouble>(3*3);
	this->FaceOnSurfaceIndices(&indices[0],&indices[1],&indices[2]);
	for(int i=0;i<3;i++) for(int j=0;j<3;j++) xyz_list_edge[i*3+j]=xyz_list[indices[i]][j];

	/*Assign output pointer*/
	*pxyz_list = xyz_list_edge;

}/*}}}*/
bool     Tetra::HasFaceOnBase(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonbaseEnum);
	sum = values[0]+values[1]+values[2]+values[3];

	_assert_(sum==0. || sum==1. || sum==2. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
bool     Tetra::HasFaceOnSurface(){/*{{{*/

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	Element::GetInputListOnVertices(&values[0],MeshVertexonsurfaceEnum);
	sum = values[0]+values[1]+values[2]+values[3];

	_assert_(sum==0. || sum==1. || sum==2. || sum==3.);

	if(sum==3){
		return true;
	}
	else{
		return false;
	}
}
/*}}}*/
void     Tetra::InputUpdateFromIoModel(int index,IoModel* iomodel){ /*{{{*/

	/*Intermediaries*/
	int         i,j;
	int         tetra_vertex_ids[NUMVERTICES];
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
		tetra_vertex_ids[i]=iomodel->elements[NUMVERTICES*index+i]; //ids for vertices are in the elements array from Matlab
	}
}
/*}}}*/
void     Tetra::InputUpdateFromSolutionOneDof(IssmDouble* solution,int enum_type){/*{{{*/

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
		if(xIsInf<IssmDouble>(values[i])) _error_("Inf found in solution vector");
	}

	/*Add input to the element: */
	this->AddInput(enum_type,values,this->element_type);

	/*Free resources:*/
	xDelete<IssmDouble>(values);
	xDelete<int>(doflist);
}
/*}}}*/
bool     Tetra::IsIcefront(void){/*{{{*/

	/*Retrieve all inputs and parameters*/
	IssmDouble ls[NUMVERTICES];
	Element::GetInputListOnVertices(&ls[0],MaskIceLevelsetEnum);

	/* If only one vertex has ice, there is an ice front here */
	if(IsIceInElement()){
		int nrice=0;       
		for(int i=0;i<NUMVERTICES;i++) if(ls[i]<0.) nrice++;
		if(nrice==1) return true;
	}
	return false;
}/*}}}*/
void     Tetra::JacobianDeterminant(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminant(pJdet,xyz_list,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantBase(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list_base,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantSurface(IssmDouble* pJdet,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list,(GaussTetra*)gauss);

}
/*}}}*/
void     Tetra::JacobianDeterminantTop(IssmDouble* pJdet,IssmDouble* xyz_list_base,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetJacobianDeterminantFace(pJdet,xyz_list_base,(GaussTetra*)gauss);

}
/*}}}*/
Gauss*   Tetra::NewGauss(void){/*{{{*/
	return new GaussTetra();
}
/*}}}*/
Gauss*   Tetra::NewGauss(int order){/*{{{*/
	return new GaussTetra(order);
}
/*}}}*/
Gauss*   Tetra::NewGauss(IssmDouble* xyz_list, IssmDouble* xyz_list_front,int order_horiz,int order_vert){/*{{{*/
	/*FIXME: this is messed up, should provide indices and not xyz_list!*/
	int indices[3];
	this->FaceOnFrontIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],max(order_horiz,order_vert));
}
/*}}}*/
Gauss*   Tetra::NewGaussBase(int order){/*{{{*/

	int indices[3];
	this->FaceOnBaseIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],order);
}
/*}}}*/
Gauss*   Tetra::NewGaussTop(int order){/*{{{*/

	int indices[3];
	this->FaceOnSurfaceIndices(&indices[0],&indices[1],&indices[2]);
	return new GaussTetra(indices[0],indices[1],indices[2],order);
}
/*}}}*/
void     Tetra::NodalFunctions(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->element_type);

}
/*}}}*/
void     Tetra::NodalFunctionsDerivatives(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTetra*)gauss,this->element_type);

}
/*}}}*/
void     Tetra::NodalFunctionsDerivativesVelocity(IssmDouble* dbasis,IssmDouble* xyz_list,Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctionsDerivatives(dbasis,xyz_list,(GaussTetra*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsPressure(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->PressureInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsTensor(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->TensorInterpolation());

}
/*}}}*/
void     Tetra::NodalFunctionsVelocity(IssmDouble* basis, Gauss* gauss){/*{{{*/

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->GetNodalFunctions(basis,(GaussTetra*)gauss,this->VelocityInterpolation());

}
/*}}}*/
void     Tetra::NormalBase(IssmDouble* bed_normal,IssmDouble* xyz_list){/*{{{*/

	IssmDouble v13[3],v23[3];
	IssmDouble normal[3];
	IssmDouble normal_norm;

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0*3+i]-xyz_list[2*3+i];
		v23[i]=xyz_list[1*3+i]-xyz_list[2*3+i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];
	normal_norm=sqrt(normal[0]*normal[0]+ normal[1]*normal[1]+ normal[2]*normal[2]);

	/*Bed normal is opposite to surface normal*/
	bed_normal[0]=-normal[0]/normal_norm;
	bed_normal[1]=-normal[1]/normal_norm;
	bed_normal[2]=-normal[2]/normal_norm;

	_assert_(bed_normal[2]<0.);
}
/*}}}*/
void     Tetra::NormalSection(IssmDouble* normal,IssmDouble* xyz_list){/*{{{*/

	/*Build unit outward pointing vector*/
	IssmDouble AB[3];
	IssmDouble AC[3];
	IssmDouble norm;

	AB[0]=xyz_list[1*3+0] - xyz_list[0*3+0];
	AB[1]=xyz_list[1*3+1] - xyz_list[0*3+1];
	AB[2]=xyz_list[1*3+2] - xyz_list[0*3+2];
	AC[0]=xyz_list[2*3+0] - xyz_list[0*3+0];
	AC[1]=xyz_list[2*3+1] - xyz_list[0*3+1];
	AC[2]=xyz_list[2*3+2] - xyz_list[0*3+2];

	cross(normal,AB,AC);
	norm=sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

	for(int i=0;i<3;i++) normal[i]=normal[i]/(norm+1e-10);
}
/*}}}*/
void     Tetra::NormalTop(IssmDouble* top_normal,IssmDouble* xyz_list){/*{{{*/

	IssmDouble v13[3],v23[3];
	IssmDouble normal[3];
	IssmDouble normal_norm;

	for(int i=0;i<3;i++){
		v13[i]=xyz_list[0*3+i]-xyz_list[2*3+i];
		v23[i]=xyz_list[1*3+i]-xyz_list[2*3+i];
	}

	normal[0]=v13[1]*v23[2]-v13[2]*v23[1];
	normal[1]=v13[2]*v23[0]-v13[0]*v23[2];
	normal[2]=v13[0]*v23[1]-v13[1]*v23[0];
	normal_norm=sqrt(normal[0]*normal[0]+ normal[1]*normal[1]+ normal[2]*normal[2]);

	top_normal[0]=normal[0]/normal_norm;
	top_normal[1]=normal[1]/normal_norm;
	top_normal[2]=normal[2]/normal_norm;
	_assert_(top_normal[2]>0.);
}
/*}}}*/
int      Tetra::NumberofNodesPressure(void){/*{{{*/
	return TetraRef::NumberofNodes(this->PressureInterpolation());
}
/*}}}*/
int      Tetra::NumberofNodesVelocity(void){/*{{{*/
	return TetraRef::NumberofNodes(this->VelocityInterpolation());
}
/*}}}*/
int      Tetra::ObjectEnum(void){/*{{{*/

	return TetraEnum;

}/*}}}*/
int      Tetra::PressureInterpolation(void){/*{{{*/
	return TetraRef::PressureInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::ReduceMatrices(ElementMatrix* Ke,ElementVector* pe){/*{{{*/

	if(pe){
		if(this->element_type==MINIcondensedEnum){
			int indices[3]={12,13,14};
			pe->StaticCondensation(Ke,3,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[4]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<4;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			pe->StaticCondensation(Ke,size,indices);
			xDelete<int>(indices);
		}
	}

	if(Ke){
		if(this->element_type==MINIcondensedEnum){
			int indices[3]={12,13,14};
			Ke->StaticCondensation(3,&indices[0]);
		}
		else if(this->element_type==P1bubblecondensedEnum){
			int size   = nodes[4]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int offset = 0;
			for(int i=0;i<4;i++) offset+=nodes[i]->GetNumberOfDofs(NoneApproximationEnum,GsetEnum);
			int* indices=xNew<int>(size);
			for(int i=0;i<size;i++) indices[i] = offset+i;
			Ke->StaticCondensation(size,indices);
			xDelete<int>(indices);
		}
	}
}
/*}}}*/
void     Tetra::ResetFSBasalBoundaryCondition(void){/*{{{*/

	int numnodes = this->GetNumberOfNodes();

	int          approximation;
	IssmDouble*  vertexonbase= NULL;
	IssmDouble   slopex,slopey,groundedice;
	IssmDouble   xz_plane[6];

	/*For FS only: we want the CS to be tangential to the bedrock*/
	this->Element::GetInputValue(&approximation,ApproximationEnum);
	if(!HasNodeOnBase() ||  approximation!=FSApproximationEnum) return;

	//printf("element number %i \n",this->id);
	/*Get inputs*/
	Input* slopex_input=this->GetInput(BedSlopeXEnum); _assert_(slopex_input);
	Input* slopey_input=this->GetInput(BedSlopeYEnum); _assert_(slopey_input);
	Input* groundedicelevelset_input=this->GetInput(MaskOceanLevelsetEnum); _assert_(groundedicelevelset_input);
	vertexonbase = xNew<IssmDouble>(numnodes);
	this->GetInputListOnNodesVelocity(&vertexonbase[0],MeshVertexonbaseEnum);

	/*Loop over basal nodes and update their CS*/
	GaussTetra* gauss = new GaussTetra();
	for(int i=0;i<this->NumberofNodesVelocity();i++){

		if(vertexonbase[i]==1){
			gauss->GaussNode(this->VelocityInterpolation(),i);

			slopex_input->GetInputValue(&slopex,gauss);
			slopey_input->GetInputValue(&slopey,gauss);
			groundedicelevelset_input->GetInputValue(&groundedice,gauss);

			/*New X axis          New Z axis*/
			xz_plane[0]=1.;       xz_plane[3]=-slopex;  
			xz_plane[1]=0.;       xz_plane[4]=-slopey;  
			xz_plane[2]=slopex;   xz_plane[5]=1.;          

			if(groundedice>0){
				if(this->nodes[i]->GetApproximation()==FSvelocityEnum){
					this->nodes[i]->DofInSSet(2); //vz 
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[i]->GetApproximation())<<" not supported yet");
			}
			else{
				if(this->nodes[i]->GetApproximation()==FSvelocityEnum){
					this->nodes[i]->DofInFSet(2); //vz
				}
				else _error_("Flow equation approximation"<<EnumToStringx(this->nodes[i]->GetApproximation())<<" not supported yet");
			}

			XZvectorsToCoordinateSystem(&this->nodes[i]->coord_system[0][0],&xz_plane[0]);
		}
	}

	/*cleanup*/
	xDelete<IssmDouble>(vertexonbase);
	delete gauss;
}
/*}}}*/
void     Tetra::ResetHooks(){/*{{{*/

	if(this->nodes) xDelete<Node*>(this->nodes);
	this->nodes=NULL;
	this->vertices=NULL;
	this->material=NULL;
	this->parameters=NULL;

	//deal with ElementHook mother class
	for(int i=0;i<this->numanalyses;i++) if(this->hnodes[i]) this->hnodes[i]->reset();
	this->hvertices->reset();
	this->hmaterial->reset();
	if(this->hneighbors) this->hneighbors->reset();
}
/*}}}*/
void     Tetra::SetCurrentConfiguration(Elements* elementsin, Loads* loadsin, Nodes* nodesin, Materials* materialsin, Parameters* parametersin){/*{{{*/

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
Element* Tetra::SpawnBasalElement(bool depthaverage_materials){/*{{{*/

	_assert_(HasFaceOnBase());

	int index1,index2,index3;
	this->FaceOnBaseIndices(&index1,&index2,&index3);
	return SpawnTria(index1,index2,index3);
}/*}}}*/
Element* Tetra::SpawnTopElement(void){/*{{{*/

	_assert_(HasFaceOnSurface());

	int index1,index2,index3;
	this->FaceOnSurfaceIndices(&index1,&index2,&index3);
	return SpawnTria(index1,index2,index3);
}/*}}}*/
bool       Tetra::IsSpawnedElement(void){/*{{{*/

   /*Tetras cannot be collapsed elements*/
   return false;

}/*}}}*/
Tria*    Tetra::SpawnTria(int index1,int index2,int index3){/*{{{*/

	int analysis_counter;

	/*go into parameters and get the analysis_counter: */
	this->parameters->FindParam(&analysis_counter,AnalysisCounterEnum);

	/*Create Tria*/
	Tria* tria=new Tria();
	tria->id=this->id;
	tria->parameters=this->parameters;
	tria->element_type=P1Enum; //Only P1 CG for now (TO BE CHANGED)
	this->SpawnTriaHook(xDynamicCast<ElementHook*>(tria),index1,index2,index3);

	/*Spawn material*/
	tria->material=(Material*)this->material->copy2(tria);

	/*recover nodes, material*/
	tria->nodes    = (Node**)tria->hnodes[analysis_counter]->deliverp();
	tria->vertices = (Vertex**)tria->hvertices->deliverp();

	/*Return new Tria*/
	return tria;
}
/*}}}*/
int      Tetra::TensorInterpolation(void){/*{{{*/
	return TetraRef::TensorInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::Update(Inputs* inputs,int index,IoModel* iomodel,int analysis_counter,int analysis_type,int finiteelement_type){ /*{{{*/

	/*Intermediaries*/
	int        i;
	int        tetra_vertex_ids[6];
	IssmDouble nodeinputs[6];
	IssmDouble yts;
	bool       dakota_analysis;
	bool       isFS;
	int        numnodes;
	int*       tetra_node_ids = NULL;

	/*Fetch parameters: */
	iomodel->FindConstant(&yts,"md.constants.yts");
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");
	iomodel->FindConstant(&isFS,"md.flowequation.isFS");

	/*Checks if debuging*/
	_assert_(iomodel->elements);
	_assert_(index==this->sid); 

	/*Recover element type*/
	this->element_type_list[analysis_counter]=finiteelement_type;

	/*Recover vertices ids needed to initialize inputs*/
	for(i=0;i<4;i++) tetra_vertex_ids[i]=iomodel->elements[4*index+i]; //ids for vertices are in the elements array from Matlab

	/*Recover nodes ids needed to initialize the node hook.*/
	switch(finiteelement_type){
		case P1Enum:
			numnodes         = 4;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			break;
		case P1bubbleEnum: case P1bubblecondensedEnum:
			numnodes         = 5;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->numberofvertices+index+1;
			break;
		case P2Enum:
			numnodes        = 10;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;
			break;
		case MINIEnum: case MINIcondensedEnum:
			numnodes         = 9;
			tetra_node_ids   = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->numberofvertices+index+1;

			tetra_node_ids[5]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+0];
			tetra_node_ids[6]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+1];
			tetra_node_ids[7]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+2];
			tetra_node_ids[8]=iomodel->numberofvertices+iomodel->numberofelements+iomodel->elements[4*index+3];
			break;
		case TaylorHoodEnum:
		case XTaylorHoodEnum:
			numnodes        = 14;
			tetra_node_ids  = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;

			tetra_node_ids[10]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+0];
			tetra_node_ids[11]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+1];
			tetra_node_ids[12]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+2];
			tetra_node_ids[13]=iomodel->numberofvertices+iomodel->numberofedges+iomodel->elements[4*index+3];
			break;
		case LATaylorHoodEnum:
			numnodes        = 10;
			tetra_node_ids  = xNew<int>(numnodes);
			tetra_node_ids[0]=iomodel->elements[4*index+0];
			tetra_node_ids[1]=iomodel->elements[4*index+1];
			tetra_node_ids[2]=iomodel->elements[4*index+2];
			tetra_node_ids[3]=iomodel->elements[4*index+3];
			tetra_node_ids[4]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+0]+1;
			tetra_node_ids[5]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+1]+1;
			tetra_node_ids[6]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+2]+1;
			tetra_node_ids[7]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+3]+1;
			tetra_node_ids[8]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+4]+1;
			tetra_node_ids[9]=iomodel->numberofvertices+iomodel->elementtoedgeconnectivity[6*index+5]+1;
			break;
		default:
			_error_("Finite element "<<EnumToStringx(finiteelement_type)<<" not supported yet");
	}

	/*hooks: */
	this->SetHookNodes(tetra_node_ids,numnodes,analysis_counter); this->nodes=NULL;
	xDelete<int>(tetra_node_ids);

	/*Fill with IoModel*/
	this->InputUpdateFromIoModel(index,iomodel);
}
/*}}}*/
void     Tetra::ValueP1OnGauss(IssmDouble* pvalue,IssmDouble* values,Gauss* gauss){/*{{{*/
	TetraRef::GetInputValue(pvalue,values,gauss,P1Enum);
}
/*}}}*/
int      Tetra::VelocityInterpolation(void){/*{{{*/
	return TetraRef::VelocityInterpolation(this->element_type);
}
/*}}}*/
void     Tetra::ViscousHeating(IssmDouble* pphi,IssmDouble* xyz_list,Gauss* gauss,Input* vx_input,Input* vy_input,Input* vz_input){/*{{{*/

	/*Intermediaries*/
	IssmDouble phi;
	IssmDouble viscosity;
	IssmDouble epsilon[6];

	_assert_(gauss->Enum()==GaussTetraEnum);
	this->StrainRateFS(&epsilon[0],xyz_list,(GaussTetra*)gauss,vx_input,vy_input,vz_input);
	this->material->ViscosityFS(&viscosity,3,xyz_list,(GaussTetra*)gauss,vx_input,vy_input,vz_input);
	GetPhi(&phi,&epsilon[0],viscosity);

	/*Assign output pointer*/
	*pphi = phi;
}
/*}}}*/
