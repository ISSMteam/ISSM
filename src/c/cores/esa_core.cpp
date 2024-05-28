/*!\file: esa_core.cpp
 * \brief: core of the ESA solution 
 */ 

#include "./cores.h"
#include "../toolkits/toolkits.h"
#include "../classes/classes.h"
#include "../shared/shared.h"
#include "../modules/modules.h"
#include "../solutionsequences/solutionsequences.h"

void esa_core(FemModel* femmodel){ /*{{{*/

	/*Start profiler*/
	femmodel->profiler->Start(ESACORE);

	Vector<IssmDouble> *U_radial  = NULL; 
	Vector<IssmDouble> *U_north   = NULL; 
	Vector<IssmDouble> *U_east    = NULL; 
	Vector<IssmDouble> *U_x   = NULL; 
	Vector<IssmDouble> *U_y    = NULL; 
	bool save_results,isesa;
	int iscoupler;
	int domaintype;
	int solution_type;
	int        numoutputs        = 0;
	char     **requested_outputs = NULL;

	/*additional parameters: */
	int  gsize;
	bool spherical=true;
	IssmDouble          *latitude   = NULL;
	IssmDouble          *longitude  = NULL;
	IssmDouble          *radius     = NULL;
	IssmDouble          *xx     = NULL;
	IssmDouble          *yy     = NULL;
	IssmDouble          *zz     = NULL;

	/*Recover some parameters: */
	femmodel->parameters->FindParam(&domaintype,DomainTypeEnum);
	femmodel->parameters->FindParam(&solution_type,SolutionTypeEnum);
	femmodel->parameters->FindParam(&save_results,SaveResultsEnum);
	femmodel->parameters->FindParam(&isesa,TransientIsesaEnum);
	femmodel->parameters->FindParam(&iscoupler,IsSlcCouplingEnum);

	/* recover coordinates of vertices: */
	VertexCoordinatesx(&latitude,&longitude,&radius,femmodel->vertices,spherical); 
	VertexCoordinatesx(&xx,&yy,&zz,femmodel->vertices); 

	/*Figure out size of g-set deflection vector and allocate solution vector: */
	gsize      = femmodel->nodes->NumberOfDofs(GsetEnum);

	/*several cases here, depending on value of iscoupler and isesa: 
	solution_type == EsaSolutionEnum)       we are running elastic adjustment core (no coupler)
	( !iscoupler & !isesa)       we are not interested in being here :) 
	( !iscoupler & isesa)        we are running in uncoupled mode
	( iscoupler & isesa)         we are running in coupled mode, and better be earth
	( iscoupler & !isesa)        we are running in coupled mode, and better be an ice cap
	*/

	if(solution_type==EsaSolutionEnum){
		isesa=1;
		iscoupler=0;
	}

	/*early return: */
	if(!iscoupler & !isesa) return;  //we are not interested in being here :) 

	/*In what follows we assume we are all running esa, either in coupled, or uncoupled mode:*/
	if(VerboseSolution()) _printf0_("   computing elastic adjustment\n");

	/*set configuration: */
	if(isesa)femmodel->SetCurrentConfiguration(EsaAnalysisEnum);

	if(VerboseSolution()) _printf0_("   computing elastic geodetic core\n");
	if(isesa){

		/*compute components of 3-D crustal motion: */
		/*Initialize:*/
		U_radial = new Vector<IssmDouble>(gsize);
		U_north = new Vector<IssmDouble>(gsize);
		U_east = new Vector<IssmDouble>(gsize);
		U_x = new Vector<IssmDouble>(gsize);
		U_y = new Vector<IssmDouble>(gsize);

		/*call the geodetic main modlule:*/ 
		if(domaintype==Domain3DsurfaceEnum){
			femmodel->EsaGeodetic3D(U_radial,U_north,U_east,latitude,longitude,radius,xx,yy,zz); 
		}
		if(domaintype==Domain2DhorizontalEnum){
			femmodel->EsaGeodetic2D(U_radial,U_north,U_east,U_x,U_y,xx,yy); 
			InputUpdateFromVectorx(femmodel,U_x,EsaXmotionEnum,VertexSIdEnum);
			InputUpdateFromVectorx(femmodel,U_y,EsaYmotionEnum,VertexSIdEnum);
		}

		/*get results into elements:*/
		InputUpdateFromVectorx(femmodel,U_radial,EsaUmotionEnum,VertexSIdEnum);	// radial displacement 
		InputUpdateFromVectorx(femmodel,U_north,EsaNmotionEnum,VertexSIdEnum);	// north motion 
		InputUpdateFromVectorx(femmodel,U_east,EsaEmotionEnum,VertexSIdEnum);		// east motion 

		if(save_results){
			femmodel->parameters->FindParam(&requested_outputs,&numoutputs,EsaRequestedOutputsEnum);
			femmodel->RequestedOutputsx(&femmodel->results,requested_outputs,numoutputs);
		}

		if(solution_type==EsaSolutionEnum)femmodel->RequestedDependentsx();

		/*Free resources:*/	
		delete U_radial;
		delete U_north;
		delete U_east;
		delete U_x; 
		delete U_y; 
		if(numoutputs){for(int i=0;i<numoutputs;i++){xDelete<char>(requested_outputs[i]);} xDelete<char*>(requested_outputs);}
	}

	/*End profiler*/
	femmodel->profiler->Stop(ESACORE);

} 
/*}}}*/
