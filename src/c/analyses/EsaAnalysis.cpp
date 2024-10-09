#include "./EsaAnalysis.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"

/*Model processing*/
void EsaAnalysis::CreateConstraints(Constraints* constraints,IoModel* iomodel){/*{{{*/
	/*No constraints*/
}/*}}}*/
void EsaAnalysis::CreateLoads(Loads* loads, IoModel* iomodel){/*{{{*/
	/*No loads*/
}/*}}}*/
void EsaAnalysis::CreateNodes(Nodes* nodes,IoModel* iomodel,bool isamr){/*{{{*/
	::CreateNodes(nodes,iomodel,EsaAnalysisEnum,P1Enum);
}/*}}}*/
int  EsaAnalysis::DofsPerNode(int** doflist,int domaintype,int approximation){/*{{{*/
	return 1;
}/*}}}*/
void EsaAnalysis::UpdateElements(Elements* elements,Inputs* inputs,IoModel* iomodel,int analysis_counter,int analysis_type){/*{{{*/

	/*Update elements: */
	int counter=0;
	for(int i=0;i<iomodel->numberofelements;i++){
		if(iomodel->my_elements[i]){
			Element* element=(Element*)elements->GetObjectByOffset(counter);
			element->Update(inputs,i,iomodel,analysis_counter,analysis_type,P1Enum);
			counter++;
		}
	}

	/*Create inputs: */
	iomodel->FetchDataToInput(inputs,elements,"md.mask.ice_levelset",MaskIceLevelsetEnum);
	iomodel->FetchDataToInput(inputs,elements,"md.esa.deltathickness",DeltaIceThicknessEnum);

}/*}}}*/
void EsaAnalysis::UpdateParameters(Parameters* parameters,IoModel* iomodel,int solution_enum,int analysis_enum){/*{{{*/

	int         nl;
	IssmDouble* love_h=NULL;
	IssmDouble* love_l=NULL;
	IssmDouble* love_k=NULL;

	IssmDouble* U_elastic = NULL;
	IssmDouble* U_elastic_local = NULL;
	IssmDouble* H_elastic = NULL;
	IssmDouble* H_elastic_local = NULL;
	IssmDouble* G_elastic = NULL;
	IssmDouble* G_elastic_local = NULL;
	int         M,m,lower_row,upper_row;
	IssmDouble  degacc=.01;
	IssmDouble  planetradius=0;
	IssmDouble  planetarea=0;

	int     numoutputs;
	char**  requestedoutputs = NULL;

	/*transition vectors: */
	IssmDouble **transitions    = NULL;
	int         *transitions_M    = NULL;
	int         *transitions_N    = NULL;
	int          ntransitions;

	/*some constant parameters: */
	parameters->AddObject(iomodel->CopyConstantObject("md.esa.hemisphere",EsaHemisphereEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.selfattraction",SolidearthSettingsSelfAttractionEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.horiz",SolidearthSettingsHorizEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.elastic",SolidearthSettingsElasticEnum));
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.settings.rotation",SolidearthSettingsRotationEnum));

	/*deal with planet radius and area: */
	parameters->AddObject(iomodel->CopyConstantObject("md.solidearth.planetradius",SolidearthPlanetRadiusEnum));
	iomodel->FetchData(&planetradius,"md.solidearth.planetradius");
	planetarea=4*PI*planetradius*planetradius;
	parameters->AddObject(new DoubleParam(SolidearthPlanetAreaEnum,planetarea));

	/*love numbers: */
	iomodel->FetchData(&love_h,&nl,NULL,"md.solidearth.lovenumbers.h");
	iomodel->FetchData(&love_l,&nl,NULL,"md.solidearth.lovenumbers.l");
	iomodel->FetchData(&love_k,&nl,NULL,"md.solidearth.lovenumbers.k");

	/*compute elastic green function for a range of angles*/
	iomodel->FetchData(&degacc,"md.esa.degacc");
	M=reCast<int,IssmDouble>(180./degacc+1.);
	U_elastic=xNew<IssmDouble>(M);
	H_elastic=xNew<IssmDouble>(M);
	G_elastic=xNew<IssmDouble>(M);

	/*compute combined legendre + love number (elastic green function:*/
	m=DetermineLocalSize(M,IssmComm::GetComm());
	GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,IssmComm::GetComm());
	U_elastic_local=xNew<IssmDouble>(m);
	H_elastic_local=xNew<IssmDouble>(m);
	G_elastic_local=xNew<IssmDouble>(m);

	/*compute U_elastic_local, H_elastic_local, and G_elastic_local {{{ */
	for(int i=lower_row;i<upper_row;i++){
		IssmDouble alpha,x;
		alpha= reCast<IssmDouble>(i)*degacc * PI / 180.0;

		U_elastic_local[i-lower_row]= 0.5*love_h[nl-1]/sin(alpha/2.0);
		H_elastic_local[i-lower_row]= -love_l[nl-1]*(nl-1) * cos(alpha/2)*(1 + 2*sin(alpha/2)) / (2*sin(alpha/2)*(1 + sin(alpha/2)));  
		G_elastic_local[i-lower_row]= -(-0.25 + love_h[nl-1] - 0.5*love_k[nl-1]*(nl-1)) / sin(alpha/2.0); // negative sign is imposed to mean g_deformed_earth minus g_initial_undeformed_earth. Ferrell defined it as the difference in g between the undeformed initial Earth and deformed Earth. 
		//IssmDouble Pn,Pn1,Pn2;
		//IssmDouble Pn_p,Pn_p1,Pn_p2;
		IssmDouble Pn = 0.; 
		IssmDouble Pn1 = 0.; 
		IssmDouble Pn2 = 0.; 
		IssmDouble Pn_p = 0.; 
		IssmDouble Pn_p1 = 0.; 
		IssmDouble Pn_p2 = 0.; 

		for (int n=0;n<nl;n++) {
			IssmDouble deltalove_U;
			IssmDouble deltalove_H;
			IssmDouble deltalove_G;
			
			deltalove_U = love_h[n]-love_h[nl-1];
			deltalove_H = love_l[n] - (love_l[nl-1]*(nl-1)/(n+1e-12));
			deltalove_G = 2*(love_h[n]-love_h[nl-1]) - (n*love_k[n]-love_k[nl-1]*(nl-1)) - love_k[n];  

			/*compute legendre polynomials: P_n(cos\theta) & d P_n(cos\theta)/ d\theta: */
			if(n==0){
				Pn=1; 
				Pn_p=0; 
			}
			else if(n==1){ 
				Pn = cos(alpha); 
				Pn_p = 1; 
			}
			else{
				Pn = ( (2*n-1)*cos(alpha)*Pn1 - (n-1)*Pn2 ) /n;
				Pn_p = ( (2*n-1)*(Pn1+cos(alpha)*Pn_p1) - (n-1)*Pn_p2 ) /n;
			}
			Pn2=Pn1; Pn1=Pn;
			Pn_p2=Pn_p1; Pn_p1=Pn_p;

			U_elastic_local[i-lower_row] += deltalove_U*Pn;		// vertical (up) displacement 
			H_elastic_local[i-lower_row] -= sin(alpha)*deltalove_H*Pn_p;		// horizontal displacements 
			G_elastic_local[i-lower_row] -= deltalove_G*Pn;		// change in gravitational acceleration => negative sign is imposed to mean g_deformed_earth minus g_initial_undeformed_earth. Ferrell defined it as the difference in g between the undeformed initial Earth and deformed Earth. 
		//IssmDouble Pn,Pn1,Pn2;
		}
	} 
	/* }}} */

	/*merge U_elastic_local into U_elastic; H_elastic_local to H_elastic:{{{*/
	int* recvcounts=xNew<int>(IssmComm::GetSize());
	int* displs=xNew<int>(IssmComm::GetSize());

	//recvcounts:
	ISSM_MPI_Allgather(&m,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,IssmComm::GetComm());

	/*displs: */
	ISSM_MPI_Allgather(&lower_row,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,IssmComm::GetComm());

	/*All gather:*/
	ISSM_MPI_Allgatherv(U_elastic_local, m, ISSM_MPI_DOUBLE, U_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
	ISSM_MPI_Allgatherv(H_elastic_local, m, ISSM_MPI_DOUBLE, H_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
	ISSM_MPI_Allgatherv(G_elastic_local, m, ISSM_MPI_DOUBLE, G_elastic, recvcounts, displs, ISSM_MPI_DOUBLE,IssmComm::GetComm());
	/*Free resources: */
	xDelete<int>(recvcounts);
	xDelete<int>(displs);

	/*}}}*/

	/*Avoid singularity at 0: */
	U_elastic[0]=U_elastic[1];
	parameters->AddObject(new DoubleVecParam(EsaUElasticEnum,U_elastic,M));
	H_elastic[0]=H_elastic[1];
	parameters->AddObject(new DoubleVecParam(EsaHElasticEnum,H_elastic,M));
	G_elastic[0]=G_elastic[1];
	parameters->AddObject(new DoubleVecParam(EsaGElasticEnum,G_elastic,M));

	/*Free resources: */
	xDelete<IssmDouble>(love_h);
	xDelete<IssmDouble>(love_l);
	xDelete<IssmDouble>(love_k);
	xDelete<IssmDouble>(U_elastic);
	xDelete<IssmDouble>(U_elastic_local);
	xDelete<IssmDouble>(H_elastic);
	xDelete<IssmDouble>(H_elastic_local);
	xDelete<IssmDouble>(G_elastic);
	xDelete<IssmDouble>(G_elastic_local);

	/*Transitions: */
	iomodel->FetchData(&transitions,&transitions_M,&transitions_N,&ntransitions,"md.esa.transitions");
	if(transitions){
		parameters->AddObject(new DoubleMatArrayParam(EsaTransitionsEnum,transitions,ntransitions,transitions_M,transitions_N));

		for(int i=0;i<ntransitions;i++){
			IssmDouble* transition=transitions[i];
			xDelete<IssmDouble>(transition);
		}
		xDelete<IssmDouble*>(transitions);
		xDelete<int>(transitions_M);
		xDelete<int>(transitions_N);
	}

	/*Requested outputs*/
	iomodel->FindConstant(&requestedoutputs,&numoutputs,"md.esa.requested_outputs");
	if(numoutputs)parameters->AddObject(new StringArrayParam(EsaRequestedOutputsEnum,requestedoutputs,numoutputs));
	iomodel->DeleteData(&requestedoutputs,numoutputs,"md.esa.requested_outputs");

}/*}}}*/

/*Finite Element Analysis*/
void           EsaAnalysis::Core(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
void           EsaAnalysis::PreCore(FemModel* femmodel){/*{{{*/
	_error_("not implemented");
}/*}}}*/
ElementVector* EsaAnalysis::CreateDVector(Element* element){/*{{{*/
	/*Default, return NULL*/
	return NULL;
}/*}}}*/
ElementMatrix* EsaAnalysis::CreateJacobianMatrix(Element* element){/*{{{*/
_error_("Not implemented");
}/*}}}*/
ElementMatrix* EsaAnalysis::CreateKMatrix(Element* element){/*{{{*/
	_error_("not implemented yet");
}/*}}}*/
ElementVector* EsaAnalysis::CreatePVector(Element* element){/*{{{*/
_error_("not implemented yet");
}/*}}}*/
void           EsaAnalysis::GetSolutionFromInputs(Vector<IssmDouble>* solution,Element* element){/*{{{*/
	   _error_("not implemented yet");
}/*}}}*/
void           EsaAnalysis::GradientJ(Vector<IssmDouble>* gradient,Element*  element,int control_type,int control_interp,int control_index){/*{{{*/
	_error_("Not implemented yet");
}/*}}}*/
void           EsaAnalysis::InputUpdateFromSolution(IssmDouble* solution,Element* element){/*{{{*/
	/*Default, do nothing*/
	return;

}/*}}}*/
void           EsaAnalysis::UpdateConstraints(FemModel* femmodel){/*{{{*/
	/*Default, do nothing*/
	return;
}/*}}}*/
