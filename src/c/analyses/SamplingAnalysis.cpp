#include "./SamplingAnalysis.h" //
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include <random>

#define FINITEELEMENT P1Enum
#define 	NUMVERTICES   3

/* Reimplementation of IsFaceOnBoundary with two edges on boundary supported */
bool IsFaceOnBoundary(Element*  element)
{

	IssmDouble values[NUMVERTICES];
	IssmDouble sum;

	/*Retrieve all inputs and parameters*/
	element->GetInputListOnVertices(&values[0],MeshVertexonboundaryEnum);
	sum = values[0]+values[1]+values[2];

	_assert_(sum==0. || sum==1. || sum==2. || sum==3.);

	if(sum>1.){
	 	return true;
	}
	else{
	 return false;
	}
}

/*Model processing*/
void SamplingAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
void SamplingAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
void SamplingAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/

	int finiteelement;
	finiteelement = FINITEELEMENT;

	/*Check in 2d*/
	if(iomodel->domaintype!=Domain2DhorizontalEnum) _error_("Only 2D horizontal domain is implemented.");

	::CreateNodes(nodes,iomodel,SamplingAnalysisEnum,finiteelement);

}/*}}}*/
int  SamplingAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void SamplingAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	int    finiteelement;

	/*Finite element type*/
	finiteelement = FINITEELEMENT;

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,finiteelement);
			counter++;
		}
	}

	/*Create inputs: */
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.mesh.vertexonboundary",MeshVertexonboundaryEnum);

  iomodel->FetchDataToInput(inputs,elements,"md.sampling.kappa",SamplingKappaEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.sampling.beta",SamplingBetaEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.sampling.tau",SamplingTauEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sample",SampleEnum,0.);
	iomodel->FetchDataToInput(inputs,elements,"md.initialization.sample",SampleNoiseEnum,0.);
	if(iomodel->solution_enum==TransientSolutionEnum) iomodel->FetchDataToInput(inputs,elements,"md.sampling.phi",SamplingPhiEnum,0.);

}/*}}}*/
void SamplingAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int     numoutputs;
	char**  requestedoutputs = NULL;

	parameters->AddObject(iomodel->CopyConstantObject("md.sampling.alpha",SamplingAlphaEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.sampling.robin",SamplingRobinEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.sampling.seed",SamplingSeedEnum));

	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.sampling.requested_outputs");
	parameters->AddObject(new IntParam(SamplingNumRequestedOutputsEnum,numoutputs));
	if(numoutputs)parameters->AddObject(new StringArrayParam(SamplingRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.sampling.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           SamplingAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           SamplingAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* SamplingAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* SamplingAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("not implemented");
}/*}}}*/
ElementMatrix* SamplingAnalysis::CreateKMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	//if(!element->IsIceInElement()) return NULL; // Ice in element not required for sampling capability

	/*Intermediaries*/
	bool robin;
	int      domaintype;

	/*compute all stiffness matrices for this element*/

	ElementMatrix* Ke = NULL;
	element->FindParam(&robin,SamplingRobinEnum);

  if(!robin) Ke=CreateKMatrixModifiedHelmholtz(element);
  else{
		ElementMatrix* Ke1=CreateKMatrixModifiedHelmholtz(element);
		ElementMatrix* Ke2=CreateKMatrixRobinBC(element);
		Ke =new ElementMatrix(Ke1,Ke2);

		delete Ke1;
		delete Ke2;
  }

	/*clean-up and return*/
	return Ke;

}/*}}}*/
ElementMatrix* SamplingAnalysis::CreateKMatrixModifiedHelmholtz(Element* element){/*{{{*/

	/* Check if ice in element */
	//if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	int dim;
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble kappa;

	/*Get problem dimension*/
	dim=2;	// Only 2D horizontal domain is implemented so far

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);
	IssmDouble*    dbasis = xNew<IssmDouble>(dim*numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	Input*	kappa_input=element->GetInput(SamplingKappaEnum); _assert_(kappa_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);
		element->NodalFunctionsDerivatives(dbasis,xyz_list,gauss);

		kappa_input->GetInputValue(&kappa,gauss); _assert_(kappa>0);

		D=gauss->weight*Jdet;

		/* Laplacian operator */
		for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += D*(dbasis[0*numnodes+j]*dbasis[0*numnodes+i] + dbasis[1*numnodes+j]*dbasis[1*numnodes+i]);
			}
 		}

		/* Identity identity */
		IssmDouble factor = D*kappa*kappa;
		for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
					Ke->values[i*numnodes+j] += factor*(basis[j]*basis[i]);
			}
 		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	xDelete<IssmDouble>(dbasis);
	delete gauss;
	return Ke;

}/*}}}*/
ElementMatrix* SamplingAnalysis::CreateKMatrixRobinBC(Element* element){/*{{{*/

	/* Check if ice in element */
	//if(!element->IsIceInElement()) return NULL;

	/*If no boundary, return NULL*/
	if(!IsFaceOnBoundary(element)) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble* xyz_list_boundary = NULL;
	IssmDouble beta;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Ke     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);
	element->GetLevelCoordinates(&xyz_list_boundary,xyz_list,MeshVertexonboundaryEnum,1.);
	Input*	beta_input=element->GetInput(SamplingBetaEnum); _assert_(beta_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(xyz_list,xyz_list_boundary,3);
	while(gauss->next()){

		element->JacobianDeterminantSurface(&Jdet,xyz_list_boundary,gauss);
		element->NodalFunctions(basis,gauss);

		beta_input->GetInputValue(&beta,gauss); _assert_(beta>=0);

		D=gauss->weight*Jdet*beta;

		for(int i=0;i<numnodes;i++){
				for(int j=0;j<numnodes;j++){
						Ke->values[i*numnodes+j] += D*(basis[j]*basis[i]);
			}
 		}

	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(xyz_list_boundary);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Ke;

}/*}}}*/
ElementVector* SamplingAnalysis::CreatePVector(Element* element){/*{{{*/
	//_error_("not supported");
	return NULL;
}/*}}}*/
void           SamplingAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	element->GetSolutionFromInputsOneDof(solution,SampleEnum);
}/*}}}*/
void           SamplingAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("not supported");
}/*}}}*/
void           SamplingAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/

    /*Fetch number of nodes and dof for this finite element*/
    int numnodes = element->GetNumberOfNodes();

		/*Fetch dof list and allocate solution vector*/
		int* doflist=NULL;
    element->GetDofListLocal(&doflist,NoneApproximationEnum,GsetEnum);
    IssmDouble* newsample = xNew<IssmDouble>(numnodes);
		IssmDouble* tau = xNew<IssmDouble>(numnodes);
		element->GetInputListOnNodes(&tau[0],SamplingTauEnum);

    /*Use the dof list to index into the solution vector: */
    for(int i=0;i<numnodes;i++){
			newsample[i]=solution[doflist[i]] / tau[i]; // new

			/*Check solution*/
			if(xIsNan<IssmDouble>(newsample[i])) _error_("NaN found in solution vector");
			if(xIsInf<IssmDouble>(newsample[i])) _error_("Inf found in solution vector");
		}

    /*Add sample inputs to the tria element: */
    element->AddInput(SampleEnum,newsample,element->GetElementType());

    /*Free resources:*/
    xDelete<IssmDouble>(newsample);
		xDelete<IssmDouble>(tau);
    xDelete<int>(doflist);
 }/*}}}*/
void           SamplingAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
  return;
}/*}}}*/

ElementMatrix* SamplingAnalysis::CreateMassMatrix(Element* element){/*{{{*/

	/* Check if ice in element */
	//if(!element->IsIceInElement()) return NULL;

	/*Intermediaries*/
	IssmDouble  D,Jdet;
	IssmDouble* xyz_list = NULL;

	/*Fetch number of nodes and dof for this finite element*/
	int numnodes = element->GetNumberOfNodes();

	/*Initialize Element vector and other vectors*/
	ElementMatrix* Me     = element->NewElementMatrix();
	IssmDouble*    basis  = xNew<IssmDouble>(numnodes);

	/*Retrieve all inputs and parameters*/
	element->GetVerticesCoordinates(&xyz_list);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss(2);
	while(gauss->next()){

		element->JacobianDeterminant(&Jdet,xyz_list,gauss);
		element->NodalFunctions(basis,gauss);

		D=gauss->weight*Jdet;
		TripleMultiply(basis,1,numnodes,1,
					&D,1,1,0,
					basis,1,numnodes,0,
					&Me->values[0],1);
	}

	/*Clean up and return*/
	xDelete<IssmDouble>(xyz_list);
	xDelete<IssmDouble>(basis);
	delete gauss;
	return Me;
}/*}}}*/
void           SamplingAnalysis::LumpedKMatrix(Vector<IssmDouble>** pKlff,FemModel* femmodel){/*{{{*/

	/*Initialize Lumped mass matrix (actually we just save its diagonal)*/
	int fsize      = femmodel->nodes->NumberOfDofs(FsetEnum);
	int flocalsize = femmodel->nodes->NumberOfDofsLocal(FsetEnum);
	Vector<IssmDouble>* Klff = new Vector<IssmDouble>(flocalsize,fsize);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* KLe     = this->CreateKMatrix(element);
		if(KLe){
			KLe->Lump();
			KLe->AddDiagonalToGlobal(Klff);
		}
		delete KLe;
	}
	Klff->Assemble();

	/*Assign output pointer*/
	*pKlff=Klff;
}/*}}}*/
void           SamplingAnalysis::LumpedMassMatrix(Vector<IssmDouble>** pMlff,FemModel* femmodel){/*{{{*/

	/*Initialize Lumped mass matrix (actually we just save its diagonal)*/
	int fsize      = femmodel->nodes->NumberOfDofs(FsetEnum);
	int flocalsize = femmodel->nodes->NumberOfDofsLocal(FsetEnum);
	Vector<IssmDouble>* Mlff = new Vector<IssmDouble>(flocalsize,fsize);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->Lump();
			MLe->AddDiagonalToGlobal(Mlff);
		}
		delete MLe;
	}
	Mlff->Assemble();

	/*Assign output pointer*/
	*pMlff=Mlff;
}/*}}}*/
void           SamplingAnalysis::MassMatrix(Matrix<IssmDouble>** pMff,FemModel* femmodel){/*{{{*/

	/*Initialize Mass matrix*/
	Matrix<IssmDouble> *Mff = NULL;
	AllocateSystemMatricesx(&Mff,NULL,NULL,NULL,femmodel);

	/*Create and assemble matrix*/
	for(int i=0;i<femmodel->elements->Size();i++){
		Element*       element = xDynamicCast<Element*>(femmodel->elements->GetObjectByOffset(i));
		ElementMatrix* MLe     = this->CreateMassMatrix(element);
		if(MLe){
			MLe->AddToGlobal(Mff);
		}
		delete MLe;
	}
	Mff->Assemble();

	/*Assign output pointer*/
	*pMff=Mff;
}/*}}}*/
void 						SamplingAnalysis::UpdateTransientSample(FemModel *	femmodel){

	for(int j=0;j<femmodel->elements->Size();j++){
			Element* element=(Element*)femmodel->elements->GetObjectByOffset(j);
			UpdateTransientSample(element);
	 }

}
void 						SamplingAnalysis::UpdateTransientSample(Element *  	element){

	/*Intermediaries */
	IssmDouble phi, sample, noise;

	/*Fetch number vertices for this element*/
	int numvertices = element->GetNumberOfVertices();

	/*Initialize new sample*/
	IssmDouble* sample_new = xNew<IssmDouble>(numvertices);

	/*Retrieve all inputs and parameters*/
	Input*	sample_input=element->GetInput(SampleOldEnum); _assert_(sample_input);
	Input*	noise_input=element->GetInput(SampleNoiseEnum); _assert_(noise_input);
	Input*	phi_input=element->GetInput(SamplingPhiEnum); _assert_(phi_input);

	/* Start  looping on the number of gaussian points: */
	Gauss* gauss=element->NewGauss();
  for(int iv=0;iv<numvertices;iv++){
  	gauss->GaussVertex(iv);

		/*Get input values at gauss points*/
		sample_input->GetInputValue(&sample,gauss);
		noise_input->GetInputValue(&noise,gauss);
		phi_input->GetInputValue(&phi,gauss);

		/*Get new sample*/
		sample_new[iv] = phi*sample + noise;

  }

  element->AddInput(SampleEnum,sample_new,element->GetElementType());

  /*Clean up and return*/
  xDelete<IssmDouble>(sample_new);
  delete gauss;

}
