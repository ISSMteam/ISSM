/*!\file GroundinglineMigrationx
 * \brief: migration grounding line position.
 */

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "./GroundinglineMigrationx.h"

void GroundinglineMigrationx(Elements* elements,Nodes* nodes, Vertices* vertices,Loads* loads,Materials* materials, Parameters* parameters){

   int         migration_style;
	IssmDouble *vertices_potentially_ungrounding = NULL;
	IssmDouble *phi_ungrounding                  = NULL;

	/*retrieve parameters: */
	parameters->FindParam(&migration_style,GroundinglineMigrationEnum);
	if(migration_style==NoneEnum) return;

	if(VerboseModule()) _printf0_("   Migrating grounding line based on "<<EnumToStringx(migration_style)<<"\n");


	switch(migration_style){
		case SoftMigrationEnum:
			ToolkitsOptionsFromAnalysis(parameters,DefaultAnalysisEnum);
			/*Create flag for grounded vertices above the hydrostatic equilibrium: */
			vertices_potentially_ungrounding=PotentialUngrounding(elements,vertices,parameters);
			/*propagate ice shelf into connex areas of the ice sheet that potentially want to unground: */
			phi_ungrounding=PropagateFloatingiceToGroundedNeighbors(elements,nodes,vertices,parameters,vertices_potentially_ungrounding);
			break;
		case ContactEnum:
			ToolkitsOptionsFromAnalysis(parameters,DefaultAnalysisEnum);
			phi_ungrounding=ContactFSLevelset(elements,vertices);
			break;
		case SubelementMigrationEnum:
		case AggressiveMigrationEnum:
		case GroundingOnlyEnum:
			/*Nothing additional to do here, MigrateGroundingLine takes care of everything*/
			break;
		default:
			_error_("Grounding line migration "<<EnumToStringx(migration_style) << " not supported yet!");
	}

	/*Migrate grounding line : */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->MigrateGroundingLine(phi_ungrounding);
	}

	/*free ressouces: */
	xDelete<IssmDouble>(vertices_potentially_ungrounding);
	xDelete<IssmDouble>(phi_ungrounding);
}

IssmDouble*    ContactFSLevelset(Elements* elements,Vertices* vertices){ /*{{{*/

	Vector<IssmDouble>* vertex_sigmann = NULL;
	Vector<IssmDouble>* vertex_waterpressure = NULL;
	IssmDouble*  serial_vertex_sigmann = NULL;
	IssmDouble*  serial_vertex_waterpressure = NULL;
	IssmDouble*  phi                   = NULL;

	/*Initialize vector with number of vertices*/
	int numberofvertices = vertices->NumberOfVertices();
	vertex_sigmann = new Vector<IssmDouble>(numberofvertices);
	vertex_waterpressure = new Vector<IssmDouble>(numberofvertices);
	phi            = xNew<IssmDouble>(numberofvertices);

	/*Fill vector vertices_potentially_floating: */
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->FSContactMigration(vertex_sigmann,vertex_waterpressure);
	}
	/*Assemble vector and serialize */
	vertex_sigmann->Assemble();
	vertex_waterpressure->Assemble();
	serial_vertex_sigmann=vertex_sigmann->ToMPISerial();
	serial_vertex_waterpressure=vertex_waterpressure->ToMPISerial();

	for(int i=0;i<numberofvertices;i++){
		if (serial_vertex_waterpressure[i] > serial_vertex_sigmann[i]) phi[i]=-1;
		else phi[i]=1;
	}

	/*free ressouces and return: */
	delete vertex_sigmann;
	delete vertex_waterpressure;
	xDelete<IssmDouble>(serial_vertex_sigmann);
	xDelete<IssmDouble>(serial_vertex_waterpressure);

	return phi;
}
/*}}}*/
IssmDouble*    PotentialUngrounding(Elements* elements,Vertices* vertices,Parameters* parameters){ /*{{{*/

	int                 i,numberofvertices;
	IssmDouble*         vertices_potentially_ungrounding      = NULL;
	Vector<IssmDouble>* vec_vertices_potentially_ungrounding  = NULL;
	Element*            element                               = NULL;

	/*Initialize vector with number of vertices*/
	numberofvertices=vertices->NumberOfVertices();
	vec_vertices_potentially_ungrounding=new Vector<IssmDouble>(numberofvertices); //grounded vertex that could start floating

	/*Fill vector vertices_potentially_floating: */
	for(Object* & object : elements->objects){
		element=xDynamicCast<Element*>(object);
		element->PotentialUngrounding(vec_vertices_potentially_ungrounding);
	}

	/*Assemble vector and serialize */
	vec_vertices_potentially_ungrounding->Assemble();
	vertices_potentially_ungrounding=vec_vertices_potentially_ungrounding->ToMPISerial();

	/*free ressouces and return: */
	delete vec_vertices_potentially_ungrounding;
	return vertices_potentially_ungrounding;
}
/*}}}*/
IssmDouble*    PropagateFloatingiceToGroundedNeighbors(Elements* elements,Nodes* nodes,Vertices* vertices,Parameters* parameters,IssmDouble* vertices_potentially_ungrounding){ /*{{{*/
	int                 i,analysis_type;
	int                 nflipped,local_nflipped;
	IssmDouble*         phi                                  = NULL;
	IssmDouble*         elements_neighboring_floatingce      = NULL;
	Vector<IssmDouble>* vec_elements_neighboring_floatingice = NULL;
	Vector<IssmDouble>* vec_phi                              = NULL;
	Element*            element                               = NULL;

	/*recover parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);

	/*recover vec_phi*/
	vec_phi=new Vector<IssmDouble>(vertices->NumberOfVertices());
	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		element->GetVectorFromInputs(vec_phi,MaskOceanLevelsetEnum,VertexPIdEnum);
	}
	vec_phi->Assemble();
	phi=vec_phi->ToMPISerial();

	nflipped=1; //bootstrap
	while(nflipped){

		/*Vector of size number of elements*/
		vec_elements_neighboring_floatingice=new Vector<IssmDouble>(elements->NumberOfElements(),true);

		/*Figure out if any of the nodes of the element will be floating -> elements neighbouting the floating ice*/
		for(Object* & object : elements->objects){
			element=xDynamicCast<Element*>(object);
			vec_elements_neighboring_floatingice->SetValue(element->Sid(),element->IsNodeOnShelfFromFlags(phi)?1.0:0.0,INS_VAL);
		}

		/*Assemble vector and serialize: */
		vec_elements_neighboring_floatingice->Assemble();
		elements_neighboring_floatingce=vec_elements_neighboring_floatingice->ToMPISerial();

		/*Go through elements_neighboring_floatingce, and update vector of the nodes that will start floating*/
		local_nflipped=0;
		for(Object* & object : elements->objects){
			element=xDynamicCast<Element*>(object);
			if(reCast<int,IssmDouble>(elements_neighboring_floatingce[element->Sid()])){
				local_nflipped+=element->UpdatePotentialUngrounding(vertices_potentially_ungrounding,vec_phi,phi);
			}
		}
		vec_phi->Assemble();

		ISSM_MPI_Allreduce(&local_nflipped,&nflipped,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
		if(VerboseConvergence()) _printf0_("   Additional number of vertices allowed to unground: " << nflipped << "\n");

		/*Avoid leaks: */
		xDelete<IssmDouble>(elements_neighboring_floatingce);
		xDelete<IssmDouble>(phi);

		/*Assemble and serialize:*/
		delete vec_elements_neighboring_floatingice;
		phi=vec_phi->ToMPISerial();
	}

	/*Free resources:*/
	delete vec_phi;
	xDelete<IssmDouble>(elements_neighboring_floatingce);

	return phi;
}
/*}}}*/
