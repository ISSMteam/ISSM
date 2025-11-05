/*!\file ModelProcessorx
 * \brief: create datasets using input binary file and a set of requested analyses
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void ModelProcessorx(Elements** pelements, Nodes*** pnodes, Vertices** pvertices, Materials** pmaterials, Constraints*** pconstraints, Loads*** ploads, Parameters** pparameters,Inputs** pinputs,IoModel* iomodel,FILE* toolkitfile, char* rootpath,const int solution_enum,const int nummodels,const int* analysis_enum_list){
	_assert_(nummodels>0);

	/*Set Verbosity once for all*/
	int verbose;
	iomodel->FindConstant(&verbose,"md.verbose");
	SetVerbosityLevel(verbose);

	if(VerboseMProcessor()) _printf0_("   starting model processor \n");

	/*Initialize datasets*/
	Elements    *elements     = new Elements();
	Vertices    *vertices     = new Vertices();
	Materials   *materials    = new Materials();
	Parameters  *parameters   = new Parameters();
	Constraints **constraints = xNew<Constraints*>(nummodels);
	Loads       **loads       = xNew<Loads*>(nummodels);
	Nodes       **nodes       = xNew<Nodes*>(nummodels);
	for(int i = 0;i<nummodels;i++) constraints[i] = new Constraints();
	for(int i = 0;i<nummodels;i++) loads[i]       = new Loads();
	for(int i = 0;i<nummodels;i++) nodes[i]       = new Nodes();

	/*Partition Elements and Nodes*/
	if (iomodel->numberofelements > 0) ElementsAndVerticesPartitioning(iomodel);

	/*Create elements, vertices and materials, independent of analysis_enum: */
	CreateElements(elements,iomodel,nummodels);
	CreateVertices(elements,vertices,iomodel,solution_enum);
	CreateParameters(parameters,iomodel,rootpath,toolkitfile,solution_enum);

	/*Should move to CreateInputs*/
	Inputs *inputs = new Inputs(elements->Size(),vertices->Size());
	if (iomodel->domaintype != Domain3DsurfaceEnum) iomodel->FetchDataToInput(inputs,elements,"md.mesh.scale_factor",MeshScaleFactorEnum,1.);

	/*Can now do Materials since we have created Inputs*/
	CreateMaterials(elements,inputs,materials,iomodel,nummodels);

	/*Update datasets based on each analysis (and add nodes, constrains and loads)*/
	for(int i=0;i<nummodels;i++){

		int analysis_enum=analysis_enum_list[i];
		parameters->AddObject(new IntParam(AnalysisCounterEnum,i));

		if(VerboseMProcessor()) _printf0_("   creating datasets for analysis " << EnumToStringx(analysis_enum) << "\n");
		Analysis* analysis = EnumToAnalysis(analysis_enum);
		analysis->UpdateParameters(parameters,iomodel,solution_enum,analysis_enum);
		analysis->CreateNodes(nodes[i],iomodel);
		analysis->UpdateElements(elements,inputs,iomodel,i,analysis_enum);
		analysis->CreateConstraints(constraints[i],iomodel);
		analysis->CreateLoads(loads[i],iomodel);
		delete analysis;

		/*Tell datasets that Ids are already sorted*/
		constraints[i]->Presort();
		loads[i]->Presort();
		nodes[i]->Presort();

		/*Finalize loads (count pengrids,penpairs,rifts,etc)*/
		loads[i]->Finalize();
	}

	/*Transient specific updates*/
	if(solution_enum==TransientSolutionEnum){
		UpdateParametersTransient(parameters,iomodel);
		UpdateElementsTransient(elements,parameters,inputs,iomodel);
	}

	/*Control specific updates*/
	if(VerboseMProcessor()) _printf0_("   updating elements and materials for control parameters" << "\n");
	UpdateElementsAndMaterialsControl(elements,parameters,inputs,materials,iomodel);

	/*Dakota specific updates*/
	#ifdef _HAVE_DAKOTA_
	if(VerboseMProcessor()) _printf0_("   updating elements and materials for uncertainty quantification" << "\n");
	UpdateElementsAndMaterialsDakota(elements,inputs,materials,iomodel);
	#endif

	/*Output definitions dataset: */
	if(VerboseMProcessor()) _printf0_("   creating output definitions" << "\n");
	CreateOutputDefinitions(elements,parameters,inputs,iomodel);

	/* Sort datasets:
	 * All our datasets are already ordered by ids. Set presort flag so that
	 * later on, when sorting is requested on these datasets, it will not be
	 * redone: */
	elements->Presort();
	vertices->Presort();
	materials->Presort();

	/*Assign output pointers:*/
	*pelements    = elements;
	*pnodes       = nodes;
	*pvertices    = vertices;
	*pmaterials   = materials;
	*pconstraints = constraints;
	*ploads       = loads;
	*pparameters  = parameters;
	*pinputs     = inputs;

	if(VerboseMProcessor()) _printf0_("   done with model processor \n");
}
