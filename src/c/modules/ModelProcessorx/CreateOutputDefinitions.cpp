/*!\file: CreateParametersOutputDefinitions.cpp
 * \brief driver for creating output definitions dataset, and including it into the parameters dataset
 */ 

#include "../../toolkits/toolkits.h"
#include "../../classes/classes.h"
#include "../../shared/shared.h"
#include "./ModelProcessorx.h"

void CreateOutputDefinitions(Elements* elements,Parameters* parameters,Inputs* inputs,IoModel* iomodel){

	int i,j;

	DataSet*     output_definitions      = NULL;
	int*         output_definition_enums = NULL;
	int          num_output_definitions;

	/*Create output_definitions dataset: */
	output_definitions=new DataSet();
	char** out_strings = NULL;
	iomodel->FetchData(&out_strings,&num_output_definitions,"md.outputdefinition.list");
	if(num_output_definitions>0){
		output_definition_enums=xNew<int>(num_output_definitions);
		for(int i=0;i<num_output_definitions;i++){
			output_definition_enums[i]=StringToEnumx(out_strings[i]);
		}
	}
	// free data:
	for(int i=0;i<num_output_definitions;i++) xDelete<char>(out_strings[i]);
	xDelete<char*>(out_strings);

	if(num_output_definitions){
		for (i=0;i<num_output_definitions;i++){
			if (output_definition_enums[i]==MassfluxatgateEnum){
				/*Deal with mass flux gates:{{{ */

				/*massfluxatgate variables: */
				int          temp,numgates;
				char       **gatenames           = NULL;
				char		  **gatedefinitionstrings = NULL;
				IssmDouble **gatesegments        = NULL;
				int         *gatesegments_M      = NULL;

				/*Fetch segments and names: */
				iomodel->FetchMultipleData(&gatenames,&numgates,                     "md.massfluxatgate.name");
				iomodel->FetchMultipleData(&gatedefinitionstrings,&temp,             "md.massfluxatgate.definitionstring"); _assert_(temp==numgates);
				iomodel->FetchMultipleData(&gatesegments,&gatesegments_M,NULL,&temp, "md.massfluxatgate.segments");         _assert_(temp==numgates);

				for(j=0;j<numgates;j++){
					output_definitions->AddObject(new Massfluxatgate<IssmDouble>(gatenames[j],StringToEnumx(gatedefinitionstrings[j]),gatesegments_M[j],gatesegments[j]));
				}
				/*Free resources:*/
				for(j=0;j<numgates;j++){
					char*       string  = gatenames[j];             xDelete<char>(string);
					char*       string2 = gatedefinitionstrings[j]; xDelete<char>(string2);
					IssmDouble* gate    = gatesegments[j];          xDelete<IssmDouble>(gate);
				}
				xDelete<char*>(gatenames);
				xDelete<IssmDouble*>(gatesegments);
				xDelete<int>(gatesegments_M);
				xDelete<char*>(gatedefinitionstrings);
				/*}}}*/
			}
			else if (output_definition_enums[i]==MisfitEnum){
				/*Deal with misfits: {{{*/

				/*misfit variables: */
				int          nummisfits;
				char**       misfit_name_s						= NULL;    
				char**		 misfit_definitionstring_s		= NULL;    
				char**       misfit_model_string_s			= NULL;
				IssmDouble** misfit_observation_s			= NULL;
				char**		 misfit_observation_string_s	= NULL;
				int*         misfit_observation_M_s			= NULL;
				int*         misfit_observation_N_s			= NULL;
				int*         misfit_local_s					= NULL;
				char**       misfit_timeinterpolation_s	= NULL;
				IssmDouble** misfit_weights_s					= NULL;
				int*         misfit_weights_M_s				= NULL;
				int*         misfit_weights_N_s				= NULL;
				char**       misfit_weights_string_s		= NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/misfit.m): */
				iomodel->FetchMultipleData(&misfit_name_s,&nummisfits,                                                        "md.misfit.name");
				iomodel->FetchMultipleData(&misfit_definitionstring_s,&nummisfits,                                            "md.misfit.definitionstring");
				iomodel->FetchMultipleData(&misfit_model_string_s,&nummisfits,                                                "md.misfit.model_string");
				iomodel->FetchMultipleData(&misfit_observation_s,&misfit_observation_M_s,&misfit_observation_N_s,&nummisfits, "md.misfit.observation");
				iomodel->FetchMultipleData(&misfit_observation_string_s,&nummisfits,                                          "md.misfit.observation_string");
				iomodel->FetchMultipleData(&misfit_timeinterpolation_s,&nummisfits,                                           "md.misfit.timeinterpolation");
				iomodel->FetchMultipleData(&misfit_local_s,&nummisfits,                                                       "md.misfit.local");
				iomodel->FetchMultipleData(&misfit_weights_s,&misfit_weights_M_s,&misfit_weights_N_s,&nummisfits,             "md.misfit.weights");
				iomodel->FetchMultipleData(&misfit_weights_string_s,&nummisfits,                                              "md.misfit.weights_string");

				for(j=0;j<nummisfits;j++){

					int obs_vector_type=0;
					if ((misfit_observation_M_s[j]==iomodel->numberofvertices) || (misfit_observation_M_s[j]==iomodel->numberofvertices+1)){
						obs_vector_type=1;
					}
					else if ((misfit_observation_M_s[j]==iomodel->numberofelements) || (misfit_observation_M_s[j]==iomodel->numberofelements+1)){
						obs_vector_type=2;
					}
					else
					 _error_("misfit observation size not supported yet");

					int weight_vector_type=0;
					if ((misfit_weights_M_s[j]==iomodel->numberofvertices) || (misfit_weights_M_s[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((misfit_weights_M_s[j]==iomodel->numberofelements) || (misfit_weights_M_s[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("misfit weight size not supported yet");

					/*First create a misfit object for that specific string (misfit_model_string_s[j]):*/
					output_definitions->AddObject(new Misfit(misfit_name_s[j],StringToEnumx(misfit_definitionstring_s[j]),StringToEnumx(misfit_model_string_s[j]),StringToEnumx(misfit_observation_string_s[j]),misfit_timeinterpolation_s[j],misfit_local_s[j],StringToEnumx(misfit_weights_string_s[j])));

					/*Now, for this particular misfit object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->InputCreate(misfit_observation_s[j],inputs,iomodel,misfit_observation_M_s[j],misfit_observation_N_s[j],obs_vector_type,StringToEnumx(misfit_observation_string_s[j]),7);
						element->InputCreate(misfit_weights_s[j],inputs,iomodel,misfit_weights_M_s[j],misfit_weights_N_s[j],weight_vector_type,StringToEnumx(misfit_weights_string_s[j]),7);
					}

				}

				/*Free resources:*/
				for(j=0;j<nummisfits;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;
					string = misfit_definitionstring_s[j];		xDelete<char>(string);
					string = misfit_observation_string_s[j];	xDelete<char>(string);
					string = misfit_model_string_s[j];			xDelete<char>(string);
					string = misfit_weights_string_s[j];		xDelete<char>(string);
					string = misfit_name_s[j];    xDelete<char>(string);
					string = misfit_timeinterpolation_s[j];    xDelete<char>(string);
					matrix = misfit_observation_s[j]; xDelete<IssmDouble>(matrix);
					matrix = misfit_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(misfit_name_s);
				xDelete<char*>(misfit_model_string_s);
				xDelete<char*>(misfit_definitionstring_s);
				xDelete<IssmDouble*>(misfit_observation_s);
				xDelete<char*>(misfit_observation_string_s);
				xDelete<int>(misfit_observation_M_s);
				xDelete<int>(misfit_observation_N_s);
				xDelete<int>(misfit_local_s);
				xDelete<char*>(misfit_timeinterpolation_s);
				xDelete<IssmDouble*>(misfit_weights_s);
				xDelete<int>(misfit_weights_M_s);
				xDelete<int>(misfit_weights_N_s);
				xDelete<char*>(misfit_weights_string_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CfsurfacesquareEnum){
				/*Deal with cfsurfacesquare: {{{*/

				/*cfsurfacesquare variables: */
				int          num_cfsurfacesquares;
				char       **cfsurfacesquare_name_s               = NULL;
				char       **cfsurfacesquare_definitionstring_s   = NULL;
				int         *cfsurfacesquare_surfaceid_s          = NULL;
				char       **cfsurfacesquare_model_string_s       = NULL;
				IssmDouble **cfsurfacesquare_observation_s        = NULL;
				char       **cfsurfacesquare_observation_string_s = NULL;
				int         *cfsurfacesquare_observation_M_s      = NULL;
				int         *cfsurfacesquare_observation_N_s      = NULL;
				IssmDouble **cfsurfacesquare_weights_s            = NULL;
				int         *cfsurfacesquare_weights_M_s          = NULL;
				int         *cfsurfacesquare_weights_N_s          = NULL;
				char       **cfsurfacesquare_weights_string_s     = NULL;
				IssmDouble  *cfsurfacesquare_datatime_s           = NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfsurfacesquare.m): */
				iomodel->FetchMultipleData(&cfsurfacesquare_name_s,&num_cfsurfacesquares,             "md.cfsurfacesquare.name");
				iomodel->FetchMultipleData(&cfsurfacesquare_definitionstring_s,&num_cfsurfacesquares, "md.cfsurfacesquare.definitionstring");
				iomodel->FetchMultipleData(&cfsurfacesquare_surfaceid_s,&num_cfsurfacesquares,        "md.cfsurfacesquare.surfaceid");
				iomodel->FetchMultipleData(&cfsurfacesquare_model_string_s,&num_cfsurfacesquares,     "md.cfsurfacesquare.model_string");
				iomodel->FetchMultipleData(&cfsurfacesquare_observation_s,&cfsurfacesquare_observation_M_s,&cfsurfacesquare_observation_N_s,&num_cfsurfacesquares, "md.cfsurfacesquare.observation");
				iomodel->FetchMultipleData(&cfsurfacesquare_observation_string_s,&num_cfsurfacesquares,"md.cfsurfacesquare.observation_string");
				iomodel->FetchMultipleData(&cfsurfacesquare_weights_s,&cfsurfacesquare_weights_M_s,&cfsurfacesquare_weights_N_s,&num_cfsurfacesquares,"md.cfsurfacesquare.weights");
				iomodel->FetchMultipleData(&cfsurfacesquare_weights_string_s,&num_cfsurfacesquares,    "md.cfsurfacesquare.weights_string");
				iomodel->FetchMultipleData(&cfsurfacesquare_datatime_s,&num_cfsurfacesquares,				"md.cfsurfacesquare.datatime");

				for(j=0;j<num_cfsurfacesquares;j++){

					int obs_vector_type=0;
					if ((cfsurfacesquare_observation_M_s[j]==iomodel->numberofvertices) || (cfsurfacesquare_observation_M_s[j]==iomodel->numberofvertices+1)){
						obs_vector_type=1;
					}
					else if ((cfsurfacesquare_observation_M_s[j]==iomodel->numberofelements) || (cfsurfacesquare_observation_M_s[j]==iomodel->numberofelements+1)){
						obs_vector_type=2;
					}
					else
					 _error_("cfsurfacesquare observation size not supported yet");

					int weight_vector_type=0;
					if ((cfsurfacesquare_weights_M_s[j]==iomodel->numberofvertices) || (cfsurfacesquare_weights_M_s[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((cfsurfacesquare_weights_M_s[j]==iomodel->numberofelements) || (cfsurfacesquare_weights_M_s[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("cfsurfacesquare weight size not supported yet");

					/*First create a cfsurfacesquare object for that specific string (cfsurfacesquare_model_string_s[j]):*/
					output_definitions->AddObject(new Cfsurfacesquare(cfsurfacesquare_name_s[j],StringToEnumx(cfsurfacesquare_definitionstring_s[j]),StringToEnumx(cfsurfacesquare_model_string_s[j]),cfsurfacesquare_datatime_s[j], cfsurfacesquare_surfaceid_s[j]));

					/*Now, for this particular cfsurfacesquare object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputAdd(StringToEnumx(cfsurfacesquare_definitionstring_s[j]),cfsurfacesquare_observation_s[j],inputs,iomodel,cfsurfacesquare_observation_M_s[j],cfsurfacesquare_observation_N_s[j],obs_vector_type,StringToEnumx(cfsurfacesquare_observation_string_s[j]),SurfaceObservationEnum);
						element->DatasetInputAdd(StringToEnumx(cfsurfacesquare_definitionstring_s[j]),cfsurfacesquare_weights_s[j],inputs,iomodel,cfsurfacesquare_weights_M_s[j],cfsurfacesquare_weights_N_s[j],weight_vector_type,StringToEnumx(cfsurfacesquare_weights_string_s[j]),WeightsSurfaceObservationEnum);

					}

				}

				  /*Free resources:*/
				for(j=0;j<num_cfsurfacesquares;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = cfsurfacesquare_definitionstring_s[j];		xDelete<char>(string);
					string = cfsurfacesquare_observation_string_s[j];	xDelete<char>(string);
					string = cfsurfacesquare_model_string_s[j];			xDelete<char>(string);
					string = cfsurfacesquare_weights_string_s[j];		xDelete<char>(string);
					string = cfsurfacesquare_name_s[j];    xDelete<char>(string);
					matrix = cfsurfacesquare_observation_s[j]; xDelete<IssmDouble>(matrix);
					matrix = cfsurfacesquare_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cfsurfacesquare_name_s);
				xDelete<int>(cfsurfacesquare_surfaceid_s);
				xDelete<char*>(cfsurfacesquare_model_string_s);
				xDelete<char*>(cfsurfacesquare_definitionstring_s);
				xDelete<IssmDouble*>(cfsurfacesquare_observation_s);
				xDelete<char*>(cfsurfacesquare_observation_string_s);
				xDelete<int>(cfsurfacesquare_observation_M_s);
				xDelete<int>(cfsurfacesquare_observation_N_s);
				xDelete<IssmDouble*>(cfsurfacesquare_weights_s);
				xDelete<int>(cfsurfacesquare_weights_M_s);
				xDelete<int>(cfsurfacesquare_weights_N_s);
				xDelete<char*>(cfsurfacesquare_weights_string_s);
				xDelete<IssmDouble>(cfsurfacesquare_datatime_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CfsurfacesquaretransientEnum){
				/*Deal with cfsurfacesquaretransient: {{{*/

				/*cfsurfacesquaretransient variables: */
				int          num_cfsurfacesquaretransients,test;
				char       **cfssqt_name_s                = NULL;
				char       **cfssqt_definitionstring_s    = NULL;
				char       **cfssqt_model_string_s        = NULL;
				IssmDouble **cfssqt_observations_s        = NULL;
				int         *cfssqt_observations_M_s      = NULL;
				int         *cfssqt_observations_N_s      = NULL;
				IssmDouble **cfssqt_weights_s             = NULL;
				int         *cfssqt_weights_M_s           = NULL;
				int         *cfssqt_weights_N_s           = NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfsurfacesquaretransient.m): */
				iomodel->FetchMultipleData(&cfssqt_name_s,&num_cfsurfacesquaretransients,"md.cfsurfacesquaretransient.name");
				iomodel->FetchMultipleData(&cfssqt_definitionstring_s,&test,"md.cfsurfacesquaretransient.definitionstring"); _assert_(test==num_cfsurfacesquaretransients);
				iomodel->FetchMultipleData(&cfssqt_model_string_s,&test,"md.cfsurfacesquaretransient.model_string"); _assert_(test==num_cfsurfacesquaretransients);
				iomodel->FetchMultipleData(&cfssqt_observations_s,&cfssqt_observations_M_s,&cfssqt_observations_N_s,&test, "md.cfsurfacesquaretransient.observations"); _assert_(test==num_cfsurfacesquaretransients);
				iomodel->FetchMultipleData(&cfssqt_weights_s,&cfssqt_weights_M_s,&cfssqt_weights_N_s, &test,"md.cfsurfacesquaretransient.weights"); _assert_(test==num_cfsurfacesquaretransients);

				for(j=0;j<num_cfsurfacesquaretransients;j++){

               /*Check that we can use P1 inputs*/
					if (cfssqt_observations_M_s[j]!=(iomodel->numberofvertices+1)) _error_("observations should be a P1 time series");
               if (cfssqt_weights_M_s[j]!=iomodel->numberofvertices+1)        _error_("weights should be a P1 time series");
					_assert_(cfssqt_observations_N_s[j]>0);

					/*extract data times from last row of observations*/
					IssmDouble *datatimes = xNew<IssmDouble>(cfssqt_observations_N_s[j]);
					for(int k=0;k<cfssqt_observations_N_s[j];k++) datatimes[k] = (cfssqt_observations_s[j])[cfssqt_observations_N_s[j]*(cfssqt_weights_M_s[j]-1)+k];

					/*First create a cfsurfacesquaretransient object for that specific string (cfssqt_model_string_s[j]):*/
					output_definitions->AddObject(new Cfsurfacesquaretransient(cfssqt_name_s[j], StringToEnumx(cfssqt_definitionstring_s[j]), StringToEnumx(cfssqt_model_string_s[j]), cfssqt_observations_N_s[j],datatimes ));
					xDelete<IssmDouble>(datatimes);

					/*Now, for this particular cfsurfacesquaretransient object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputAdd(StringToEnumx(cfssqt_definitionstring_s[j]),cfssqt_observations_s[j],inputs,iomodel,cfssqt_observations_M_s[j],cfssqt_observations_N_s[j],1,SurfaceObservationEnum,SurfaceObservationEnum);
						element->DatasetInputAdd(StringToEnumx(cfssqt_definitionstring_s[j]),cfssqt_weights_s[j],inputs,iomodel,cfssqt_weights_M_s[j],cfssqt_weights_N_s[j],1,WeightsSurfaceObservationEnum,WeightsSurfaceObservationEnum);

					}
				}

				/*Free resources:*/
				for(j=0;j<num_cfsurfacesquaretransients;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;
					string = cfssqt_definitionstring_s[j];		xDelete<char>(string);
					string = cfssqt_model_string_s[j];			xDelete<char>(string);
					string = cfssqt_name_s[j];    xDelete<char>(string);
					matrix = cfssqt_observations_s[j]; xDelete<IssmDouble>(matrix);
					matrix = cfssqt_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cfssqt_name_s);
				xDelete<char*>(cfssqt_model_string_s);
				xDelete<char*>(cfssqt_definitionstring_s);
				xDelete<IssmDouble*>(cfssqt_observations_s);
				xDelete<int>(cfssqt_observations_M_s);
				xDelete<int>(cfssqt_observations_N_s);
				xDelete<IssmDouble*>(cfssqt_weights_s);
				xDelete<int>(cfssqt_weights_M_s);
				xDelete<int>(cfssqt_weights_N_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CfdragcoeffabsgradEnum){
				/*Deal with cfdragcoeffabsgrad: {{{*/

				/*cfdragcoeffabsgrad variables: */
				int          num_cfdragcoeffabsgrads;
				char**       cfdragcoeffabsgrad_name_s						= NULL;    
				char**		 cfdragcoeffabsgrad_definitionstring_s		= NULL;    
				IssmDouble** cfdragcoeffabsgrad_weights_s					= NULL;
				int*         cfdragcoeffabsgrad_weights_M_s				= NULL;
				int*         cfdragcoeffabsgrad_weights_N_s				= NULL;
				char**       cfdragcoeffabsgrad_weights_string_s		= NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfdragcoeffabsgrad.m): */
				iomodel->FetchMultipleData(&cfdragcoeffabsgrad_name_s,&num_cfdragcoeffabsgrads,                                                        "md.cfdragcoeffabsgrad.name");
				iomodel->FetchMultipleData(&cfdragcoeffabsgrad_definitionstring_s,&num_cfdragcoeffabsgrads,                                            "md.cfdragcoeffabsgrad.definitionstring");
				iomodel->FetchMultipleData(&cfdragcoeffabsgrad_weights_s,&cfdragcoeffabsgrad_weights_M_s,&cfdragcoeffabsgrad_weights_N_s,&num_cfdragcoeffabsgrads,             "md.cfdragcoeffabsgrad.weights");
				iomodel->FetchMultipleData(&cfdragcoeffabsgrad_weights_string_s,&num_cfdragcoeffabsgrads,                                              "md.cfdragcoeffabsgrad.weights_string");

				for(j=0;j<num_cfdragcoeffabsgrads;j++){

					int weight_vector_type=0;
					if ((cfdragcoeffabsgrad_weights_M_s[j]==iomodel->numberofvertices) || (cfdragcoeffabsgrad_weights_M_s[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((cfdragcoeffabsgrad_weights_M_s[j]==iomodel->numberofelements) || (cfdragcoeffabsgrad_weights_M_s[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("cfdragcoeffabsgrad weight size not supported yet");

					/*First create a cfdragcoeffabsgrad object for that specific string (cfdragcoeffabsgrad_model_string_s[j]):*/
					output_definitions->AddObject(new Cfdragcoeffabsgrad(cfdragcoeffabsgrad_name_s[j],StringToEnumx(cfdragcoeffabsgrad_definitionstring_s[j])));

					/*Now, for this particular cfdragcoeffabsgrad object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){

						Element* element=xDynamicCast<Element*>(object);

						element->DatasetInputAdd(StringToEnumx(cfdragcoeffabsgrad_definitionstring_s[j]),cfdragcoeffabsgrad_weights_s[j],inputs,iomodel,cfdragcoeffabsgrad_weights_M_s[j],cfdragcoeffabsgrad_weights_N_s[j],weight_vector_type,StringToEnumx(cfdragcoeffabsgrad_weights_string_s[j]),WeightsSurfaceObservationEnum);

					}

				}

				/*Free resources:*/
				for(j=0;j<num_cfdragcoeffabsgrads;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = cfdragcoeffabsgrad_definitionstring_s[j];		xDelete<char>(string);
					string = cfdragcoeffabsgrad_weights_string_s[j];		xDelete<char>(string);
					string = cfdragcoeffabsgrad_name_s[j];    xDelete<char>(string);
					matrix = cfdragcoeffabsgrad_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cfdragcoeffabsgrad_name_s);
				xDelete<char*>(cfdragcoeffabsgrad_definitionstring_s);
				xDelete<IssmDouble*>(cfdragcoeffabsgrad_weights_s);
				xDelete<int>(cfdragcoeffabsgrad_weights_M_s);
				xDelete<int>(cfdragcoeffabsgrad_weights_N_s);
				xDelete<char*>(cfdragcoeffabsgrad_weights_string_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CfdragcoeffabsgradtransientEnum){
				/*Deal with cfdragcoeffabsgradtransient: {{{*/

				/*cfdragcoeffabsgrad variables: */
				int          num_cfdragcoeffabsgradtransients, test;
				char**       cfdraggradt_name_s						= NULL;    
				char**		 cfdraggradt_definitionstring_s		= NULL;    
				IssmDouble** cfdraggradt_weights_s					= NULL;
				int*         cfdraggradt_weights_M_s				= NULL;
				int*         cfdraggradt_weights_N_s				= NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfdragcoeffabsgradtransient.m): */
				iomodel->FetchMultipleData(&cfdraggradt_name_s,&num_cfdragcoeffabsgradtransients,                                                        "md.cfdragcoeffabsgradtransient.name");
				iomodel->FetchMultipleData(&cfdraggradt_definitionstring_s,&num_cfdragcoeffabsgradtransients,                                            "md.cfdragcoeffabsgradtransient.definitionstring");
				iomodel->FetchMultipleData(&cfdraggradt_weights_s,&cfdraggradt_weights_M_s,&cfdraggradt_weights_N_s,&test,             "md.cfdragcoeffabsgradtransient.weights");

				for(j=0;j<num_cfdragcoeffabsgradtransients;j++){

					/*Check that we can use P1 inputs*/
					if (cfdraggradt_weights_M_s[j]!=iomodel->numberofvertices+1)  _error_("weights should be a P1 time series");

					/*extract data times from last row of observations*/
					IssmDouble *datatimes = xNew<IssmDouble>(cfdraggradt_weights_N_s[j]);
					for(int k=0;k<cfdraggradt_weights_N_s[j];k++) datatimes[k] = (cfdraggradt_weights_s[j])[cfdraggradt_weights_N_s[j]*(cfdraggradt_weights_M_s[j]-1)+k];

					 /*First create a cfdragcoeffabsgradtransient object for that specific string:*/
					output_definitions->AddObject(new Cfdragcoeffabsgradtransient(cfdraggradt_name_s[j],StringToEnumx(cfdraggradt_definitionstring_s[j]), cfdraggradt_weights_N_s[j], datatimes));

					/*Now, for this particular cfdragcoeffabsgrad object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){

						Element* element=xDynamicCast<Element*>(object);

						element->DatasetInputAdd(StringToEnumx(cfdraggradt_definitionstring_s[j]),cfdraggradt_weights_s[j],inputs,iomodel,cfdraggradt_weights_M_s[j],cfdraggradt_weights_N_s[j],1,WeightsSurfaceObservationEnum,WeightsSurfaceObservationEnum);

					}
				}

				/*Free resources:*/
				for(j=0;j<num_cfdragcoeffabsgradtransients;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = cfdraggradt_definitionstring_s[j];		xDelete<char>(string);
					string = cfdraggradt_name_s[j];    xDelete<char>(string);
					matrix = cfdraggradt_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cfdraggradt_name_s);
				xDelete<char*>(cfdraggradt_definitionstring_s);
				xDelete<IssmDouble*>(cfdraggradt_weights_s);
				xDelete<int>(cfdraggradt_weights_M_s);
				xDelete<int>(cfdraggradt_weights_N_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CfrheologybbarabsgradEnum){
				/*Deal with cfrheologybbarabsgrad: {{{*/

				/*cfrheologybbarabsgrad variables: */
				int          num_cfrheologybbarabsgrads;
				char**       cfrheologybbarabsgrad_name_s                = NULL;
				char**       cfrheologybbarabsgrad_definitionstring_s    = NULL;
				IssmDouble** cfrheologybbarabsgrad_weights_s             = NULL;
				int*         cfrheologybbarabsgrad_weights_M_s           = NULL;
				int*         cfrheologybbarabsgrad_weights_N_s           = NULL;
				char**       cfrheologybbarabsgrad_weights_string_s      = NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfrheologybbarabsgrad.m): */
				iomodel->FetchMultipleData(&cfrheologybbarabsgrad_name_s,&num_cfrheologybbarabsgrads,                                                        "md.cfrheologybbarabsgrad.name");
				iomodel->FetchMultipleData(&cfrheologybbarabsgrad_definitionstring_s,&num_cfrheologybbarabsgrads,                                            "md.cfrheologybbarabsgrad.definitionstring");
				iomodel->FetchMultipleData(&cfrheologybbarabsgrad_weights_s,&cfrheologybbarabsgrad_weights_M_s,&cfrheologybbarabsgrad_weights_N_s,&num_cfrheologybbarabsgrads,             "md.cfrheologybbarabsgrad.weights");
				iomodel->FetchMultipleData(&cfrheologybbarabsgrad_weights_string_s,&num_cfrheologybbarabsgrads,                                              "md.cfrheologybbarabsgrad.weights_string");

				for(j=0;j<num_cfrheologybbarabsgrads;j++){

					int weight_vector_type=0;
					if ((cfrheologybbarabsgrad_weights_M_s[j]==iomodel->numberofvertices) || (cfrheologybbarabsgrad_weights_M_s[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((cfrheologybbarabsgrad_weights_M_s[j]==iomodel->numberofelements) || (cfrheologybbarabsgrad_weights_M_s[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("cfrheologybbarabsgrad weight size not supported yet");

					/*First create a cfrheologybbarabsgrad object for that specific string (cfrheologybbarabsgrad_model_string_s[j]):*/
					output_definitions->AddObject(new Cfrheologybbarabsgrad(cfrheologybbarabsgrad_name_s[j],StringToEnumx(cfrheologybbarabsgrad_definitionstring_s[j]),StringToEnumx(cfrheologybbarabsgrad_weights_string_s[j])));

					/*Now, for this particular cfrheologybbarabsgrad object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){

						Element* element=xDynamicCast<Element*>(object);

						element->DatasetInputAdd(StringToEnumx(cfrheologybbarabsgrad_definitionstring_s[j]),cfrheologybbarabsgrad_weights_s[j],inputs,iomodel,cfrheologybbarabsgrad_weights_M_s[j],cfrheologybbarabsgrad_weights_N_s[j],weight_vector_type,StringToEnumx(cfrheologybbarabsgrad_weights_string_s[j]),WeightsSurfaceObservationEnum);

					}

				}
				    /*Free resources:*/
            for(j=0;j<num_cfrheologybbarabsgrads;j++){
               char* string=NULL;
               IssmDouble* matrix = NULL;

               string = cfrheologybbarabsgrad_definitionstring_s[j];    xDelete<char>(string);
               string = cfrheologybbarabsgrad_weights_string_s[j];      xDelete<char>(string);
               string = cfrheologybbarabsgrad_name_s[j];    xDelete<char>(string);
               matrix = cfrheologybbarabsgrad_weights_s[j]; xDelete<IssmDouble>(matrix);
            }
            xDelete<char*>(cfrheologybbarabsgrad_name_s);
            xDelete<char*>(cfrheologybbarabsgrad_definitionstring_s);
            xDelete<IssmDouble*>(cfrheologybbarabsgrad_weights_s);
            xDelete<int>(cfrheologybbarabsgrad_weights_M_s);
            xDelete<int>(cfrheologybbarabsgrad_weights_N_s);
            xDelete<char*>(cfrheologybbarabsgrad_weights_string_s);
            /*}}}*/
         }
			else if (output_definition_enums[i]==CfrheologybbarabsgradtransientEnum){
				/*Deal with cfrheologybbarabsgradtransient: {{{*/

				/*cfrheologybbarabsgrad variables: */
				int          num_cfrheologybbarabsgradtransients, test;
				char**       cfrheogradt_name_s                = NULL;
				char**       cfrheogradt_definitionstring_s    = NULL;
				IssmDouble** cfrheogradt_weights_s             = NULL;
				int*         cfrheogradt_weights_M_s           = NULL;
				int*         cfrheogradt_weights_N_s           = NULL;
				char**       cfrheogradt_weights_string_s      = NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cfrheologybbarabsgradtransient.m): */
				iomodel->FetchMultipleData(&cfrheogradt_name_s,&num_cfrheologybbarabsgradtransients,                                                        "md.cfrheologybbarabsgradtransient.name");
				iomodel->FetchMultipleData(&cfrheogradt_definitionstring_s,&num_cfrheologybbarabsgradtransients,                                            "md.cfrheologybbarabsgradtransient.definitionstring");
				iomodel->FetchMultipleData(&cfrheogradt_weights_s,&cfrheogradt_weights_M_s,&cfrheogradt_weights_N_s,&test,             "md.cfrheologybbarabsgradtransient.weights");

				for(j=0;j<num_cfrheologybbarabsgradtransients;j++){

					if (cfrheogradt_weights_M_s[j]!=iomodel->numberofvertices+1) _error_("weights should be a P1 time series");

					/*extract data times from last row of observations*/
					IssmDouble *datatimes = xNew<IssmDouble>(cfrheogradt_weights_N_s[j]);
					for(int k=0;k<cfrheogradt_weights_N_s[j];k++) datatimes[k] = (cfrheogradt_weights_s[j])[cfrheogradt_weights_N_s[j]*(cfrheogradt_weights_M_s[j]-1)+k];

					/*First create a cfrheologybbarabsgradtransient object for that specific string:*/
					output_definitions->AddObject(new Cfrheologybbarabsgradtransient(cfrheogradt_name_s[j],StringToEnumx(cfrheogradt_definitionstring_s[j]), cfrheogradt_weights_N_s[j], datatimes));

					/*Now, for this particular cfrheologybbarabsgrad object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){

						Element* element=xDynamicCast<Element*>(object);

						element->DatasetInputAdd(StringToEnumx(cfrheogradt_definitionstring_s[j]),cfrheogradt_weights_s[j],inputs,iomodel,cfrheogradt_weights_M_s[j],cfrheogradt_weights_N_s[j],1,WeightsSurfaceObservationEnum,WeightsSurfaceObservationEnum);

					}
				}

				/*Free resources:*/
            for(j=0;j<num_cfrheologybbarabsgradtransients;j++){
               char* string=NULL;
               IssmDouble* matrix = NULL;

               string = cfrheogradt_definitionstring_s[j];    xDelete<char>(string);
               string = cfrheogradt_name_s[j];    xDelete<char>(string);
               matrix = cfrheogradt_weights_s[j]; xDelete<IssmDouble>(matrix);
            }
            xDelete<char*>(cfrheogradt_name_s);
            xDelete<char*>(cfrheogradt_definitionstring_s);
            xDelete<IssmDouble*>(cfrheogradt_weights_s);
            xDelete<int>(cfrheogradt_weights_M_s);
            xDelete<int>(cfrheogradt_weights_N_s);
            /*}}}*/
         }
			else if (output_definition_enums[i]==CfsurfacelogvelEnum){
				/*Deal with cfsurfacelogvel: {{{*/

				/*cfsurfacelogvel variables: */
				int          num_cfsurfacelogvels;
				char       **cfsurfacelogvel_name             = NULL;
				char       **cfsurfacelogvel_definitionstring = NULL;
				IssmDouble **cfsurfacelogvel_vxobs            = NULL;
				IssmDouble **cfsurfacelogvel_vyobs            = NULL;
				char       **cfsurfacelogvel_vxobs_string     = NULL;
				char       **cfsurfacelogvel_vyobs_string     = NULL;
				int         *cfsurfacelogvel_observation_M    = NULL;
				int         *cfsurfacelogvel_observation_N    = NULL;
				IssmDouble **cfsurfacelogvel_weights          = NULL;
				int         *cfsurfacelogvel_weights_M        = NULL;
				int         *cfsurfacelogvel_weights_N        = NULL;
				char       **cfsurfacelogvel_weightstring     = NULL;
				IssmDouble  *cfsurfacelogvel_datatime         = NULL;

            /*Fetch name, modeltring, observation, observationtring, etc ... (see src/m/classes/cfsurfacelogvel.m): */
            iomodel->FetchMultipleData(&cfsurfacelogvel_name,&num_cfsurfacelogvels,"md.cfsurfacelogvel.name");
            iomodel->FetchMultipleData(&cfsurfacelogvel_definitionstring,&num_cfsurfacelogvels,"md.cfsurfacelogvel.definitionstring");
            iomodel->FetchMultipleData(&cfsurfacelogvel_vxobs,&cfsurfacelogvel_observation_M,&cfsurfacelogvel_observation_N,&num_cfsurfacelogvels,"md.cfsurfacelogvel.vxobs");
            iomodel->FetchMultipleData(&cfsurfacelogvel_vxobs_string,&num_cfsurfacelogvels,"md.cfsurfacelogvel.vxobs_string");
            iomodel->FetchMultipleData(&cfsurfacelogvel_vyobs,NULL,NULL,&num_cfsurfacelogvels,"md.cfsurfacelogvel.vyobs");
            iomodel->FetchMultipleData(&cfsurfacelogvel_vyobs_string,&num_cfsurfacelogvels,"md.cfsurfacelogvel.vyobs_string");
            iomodel->FetchMultipleData(&cfsurfacelogvel_weights,&cfsurfacelogvel_weights_M,&cfsurfacelogvel_weights_N,&num_cfsurfacelogvels,"md.cfsurfacelogvel.weights");
            iomodel->FetchMultipleData(&cfsurfacelogvel_weightstring,&num_cfsurfacelogvels,"md.cfsurfacelogvel.weights_string");
            iomodel->FetchMultipleData(&cfsurfacelogvel_datatime,&num_cfsurfacelogvels,"md.cfsurfacelogvel.datatime");

				for(j=0;j<num_cfsurfacelogvels;j++){

					int obs_vector_type=0;
					if ((cfsurfacelogvel_observation_M[j]==iomodel->numberofvertices) || (cfsurfacelogvel_observation_M[j]==iomodel->numberofvertices+1)){
						obs_vector_type=1;
					}
					else if ((cfsurfacelogvel_observation_M[j]==iomodel->numberofelements) || (cfsurfacelogvel_observation_M[j]==iomodel->numberofelements+1)){
						obs_vector_type=2;
					}
					else
					 _error_("cfsurfacelogvel observation size not supported yet");

					int weight_vector_type=0;
					if ((cfsurfacelogvel_weights_M[j]==iomodel->numberofvertices) || (cfsurfacelogvel_weights_M[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((cfsurfacelogvel_weights_M[j]==iomodel->numberofelements) || (cfsurfacelogvel_weights_M[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("cfsurfacelogvel weight size not supported yet");

					/*First create a cfsurfacelogvel object for that specific string (cfsurfacelogvel_modeltring[j]):*/
					output_definitions->AddObject(new Cfsurfacelogvel(cfsurfacelogvel_name[j],StringToEnumx(cfsurfacelogvel_definitionstring[j]),cfsurfacelogvel_datatime[j]));

					/*Now, for this particular cfsurfacelogvel object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){

						Element* element=xDynamicCast<Element*>(object);

						element->DatasetInputAdd(StringToEnumx(cfsurfacelogvel_definitionstring[j]),cfsurfacelogvel_vxobs[j],inputs,iomodel,cfsurfacelogvel_observation_M[j],cfsurfacelogvel_observation_N[j],obs_vector_type,StringToEnumx(cfsurfacelogvel_vxobs_string[j]),VxObsEnum);
							element->DatasetInputAdd(StringToEnumx(cfsurfacelogvel_definitionstring[j]),cfsurfacelogvel_vyobs[j],inputs,iomodel,cfsurfacelogvel_observation_M[j],cfsurfacelogvel_observation_N[j],obs_vector_type,StringToEnumx(cfsurfacelogvel_vyobs_string[j]),VyObsEnum);
						element->DatasetInputAdd(StringToEnumx(cfsurfacelogvel_definitionstring[j]),cfsurfacelogvel_weights[j],inputs,iomodel,cfsurfacelogvel_weights_M[j],cfsurfacelogvel_weights_N[j],weight_vector_type,StringToEnumx(cfsurfacelogvel_weightstring[j]),WeightsSurfaceObservationEnum);

					}

				}

				/*Free resources:*/
				for(j=0;j<num_cfsurfacelogvels;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = cfsurfacelogvel_definitionstring[j];		xDelete<char>(string);
					string = cfsurfacelogvel_vxobs_string[j];	xDelete<char>(string);
					string = cfsurfacelogvel_vyobs_string[j];	xDelete<char>(string);
					string = cfsurfacelogvel_weightstring[j];		xDelete<char>(string);
					string = cfsurfacelogvel_name[j];    xDelete<char>(string);
					matrix = cfsurfacelogvel_weights[j]; xDelete<IssmDouble>(matrix);
					matrix = cfsurfacelogvel_vxobs[j]; xDelete<IssmDouble>(matrix);
					matrix = cfsurfacelogvel_vyobs[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cfsurfacelogvel_name);
				xDelete<char*>(cfsurfacelogvel_definitionstring);
				xDelete<int>(cfsurfacelogvel_observation_M);
				xDelete<IssmDouble*>(cfsurfacelogvel_vxobs);
				xDelete<IssmDouble*>(cfsurfacelogvel_vyobs);
				xDelete<char*>(cfsurfacelogvel_vxobs_string);
				xDelete<char*>(cfsurfacelogvel_vyobs_string);
				xDelete<int>(cfsurfacelogvel_observation_N);
				xDelete<IssmDouble*>(cfsurfacelogvel_weights);
				xDelete<int>(cfsurfacelogvel_weights_M);
				xDelete<int>(cfsurfacelogvel_weights_N);
				xDelete<char*>(cfsurfacelogvel_weightstring);
				xDelete<IssmDouble>(cfsurfacelogvel_datatime);
				/*}}}*/
			}
			else if (output_definition_enums[i]==CflevelsetmisfitEnum){
				/*Deal with cflevelsetmisfit: {{{*/

				/*cflevelsetmisfit variables: */
				int          num_cflevelsetmisfits;
				char**       cflevelsetmisfit_name_s						= NULL;    
				char**		 cflevelsetmisfit_definitionstring_s		= NULL;    
				char**       cflevelsetmisfit_model_string_s			= NULL;
				IssmDouble** cflevelsetmisfit_observation_s			= NULL;
				char**		 cflevelsetmisfit_observation_string_s	= NULL;
				int*         cflevelsetmisfit_observation_M_s			= NULL;
				int*         cflevelsetmisfit_observation_N_s			= NULL;
				IssmDouble** cflevelsetmisfit_weights_s					= NULL;
				int*         cflevelsetmisfit_weights_M_s				= NULL;
				int*         cflevelsetmisfit_weights_N_s				= NULL;
				char**       cflevelsetmisfit_weights_string_s		= NULL;
				IssmDouble*	 cflevelsetmisfit_datatime_s				= NULL;

				/*Fetch name, model_string, observation, observation_string, etc ... (see src/m/classes/cflevelsetmisfit.m): */
				iomodel->FetchMultipleData(&cflevelsetmisfit_name_s,&num_cflevelsetmisfits,                                                        "md.cflevelsetmisfit.name");
				iomodel->FetchMultipleData(&cflevelsetmisfit_definitionstring_s,&num_cflevelsetmisfits,                                            "md.cflevelsetmisfit.definitionstring");
				iomodel->FetchMultipleData(&cflevelsetmisfit_model_string_s,&num_cflevelsetmisfits,                                                "md.cflevelsetmisfit.model_string");
				iomodel->FetchMultipleData(&cflevelsetmisfit_observation_s,&cflevelsetmisfit_observation_M_s,&cflevelsetmisfit_observation_N_s,&num_cflevelsetmisfits, "md.cflevelsetmisfit.observation");
				iomodel->FetchMultipleData(&cflevelsetmisfit_observation_string_s,&num_cflevelsetmisfits,                                          "md.cflevelsetmisfit.observation_string");
				iomodel->FetchMultipleData(&cflevelsetmisfit_weights_s,&cflevelsetmisfit_weights_M_s,&cflevelsetmisfit_weights_N_s,&num_cflevelsetmisfits,             "md.cflevelsetmisfit.weights");
				iomodel->FetchMultipleData(&cflevelsetmisfit_weights_string_s,&num_cflevelsetmisfits,                                              "md.cflevelsetmisfit.weights_string");
				iomodel->FetchMultipleData(&cflevelsetmisfit_datatime_s,&num_cflevelsetmisfits,																	 "md.cflevelsetmisfit.datatime");

				for(j=0;j<num_cflevelsetmisfits;j++){
					int obs_vector_type=0;
					if ((cflevelsetmisfit_observation_M_s[j]==iomodel->numberofvertices) || (cflevelsetmisfit_observation_M_s[j]==iomodel->numberofvertices+1)){
						obs_vector_type=1;
					}
					else if ((cflevelsetmisfit_observation_M_s[j]==iomodel->numberofelements) || (cflevelsetmisfit_observation_M_s[j]==iomodel->numberofelements+1)){
						obs_vector_type=2;
					}
					else
					 _error_("cflevelsetmisfit observation size not supported yet");

					int weight_vector_type=0;
					if ((cflevelsetmisfit_weights_M_s[j]==iomodel->numberofvertices) || (cflevelsetmisfit_weights_M_s[j]==iomodel->numberofvertices+1)){
						weight_vector_type=1;
					}
					else if ((cflevelsetmisfit_weights_M_s[j]==iomodel->numberofelements) || (cflevelsetmisfit_weights_M_s[j]==iomodel->numberofelements+1)){
						weight_vector_type=2;
					}
					else
					 _error_("cflevelsetmisfit weight size not supported yet");

					/*First create a cflevelsetmisfit object for that specific string (cflevelsetmisfit_model_string_s[j]):*/
					output_definitions->AddObject(new Cflevelsetmisfit(cflevelsetmisfit_name_s[j],StringToEnumx(cflevelsetmisfit_definitionstring_s[j]),StringToEnumx(cflevelsetmisfit_model_string_s[j]),cflevelsetmisfit_datatime_s[j]));

					/*Now, for this particular cflevelsetmisfit object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputAdd(StringToEnumx(cflevelsetmisfit_definitionstring_s[j]),cflevelsetmisfit_observation_s[j],inputs,iomodel,cflevelsetmisfit_observation_M_s[j],cflevelsetmisfit_observation_N_s[j],obs_vector_type,StringToEnumx(cflevelsetmisfit_observation_string_s[j]),LevelsetObservationEnum);
						element->DatasetInputAdd(StringToEnumx(cflevelsetmisfit_definitionstring_s[j]),cflevelsetmisfit_weights_s[j],inputs,iomodel,cflevelsetmisfit_weights_M_s[j],cflevelsetmisfit_weights_N_s[j],weight_vector_type,StringToEnumx(cflevelsetmisfit_weights_string_s[j]),WeightsLevelsetObservationEnum);
					}
				}

				  /*Free resources:*/
				for(j=0;j<num_cflevelsetmisfits;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = cflevelsetmisfit_definitionstring_s[j];		xDelete<char>(string);
					string = cflevelsetmisfit_observation_string_s[j];	xDelete<char>(string);
					string = cflevelsetmisfit_model_string_s[j];			xDelete<char>(string);
					string = cflevelsetmisfit_weights_string_s[j];		xDelete<char>(string);
					string = cflevelsetmisfit_name_s[j];    xDelete<char>(string);
					matrix = cflevelsetmisfit_observation_s[j]; xDelete<IssmDouble>(matrix);
					matrix = cflevelsetmisfit_weights_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(cflevelsetmisfit_name_s);
				xDelete<char*>(cflevelsetmisfit_model_string_s);
				xDelete<char*>(cflevelsetmisfit_definitionstring_s);
				xDelete<IssmDouble*>(cflevelsetmisfit_observation_s);
				xDelete<char*>(cflevelsetmisfit_observation_string_s);
				xDelete<int>(cflevelsetmisfit_observation_M_s);
				xDelete<int>(cflevelsetmisfit_observation_N_s);
				xDelete<IssmDouble*>(cflevelsetmisfit_weights_s);
				xDelete<int>(cflevelsetmisfit_weights_M_s);
				xDelete<int>(cflevelsetmisfit_weights_N_s);
				xDelete<char*>(cflevelsetmisfit_weights_string_s);
				xDelete<IssmDouble>(cflevelsetmisfit_datatime_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==NodalvalueEnum){
				/*Deal with nodal values: {{{*/

				/*nodal value variables: */
				int          numnodalvalues;
				char**       nodalvalue_name_s             = NULL;    
				char**       nodalvalue_definitionstrings             = NULL;    
				char**       nodalvalue_modelstrings        = NULL;
				int*         nodalvalue_node_s = NULL;

				/*Fetch name, model_enum, etc ... (see src/m/classes/nodalvalue.m): */
				iomodel->FetchMultipleData(&nodalvalue_name_s,&numnodalvalues,            "md.nodalvalue.name");
				iomodel->FetchMultipleData(&nodalvalue_definitionstrings,&numnodalvalues, "md.nodalvalue.definitionenum");
				iomodel->FetchMultipleData(&nodalvalue_modelstrings,&numnodalvalues,      "md.nodalvalue.model_enum");
				iomodel->FetchMultipleData(&nodalvalue_node_s,&numnodalvalues,            "md.nodalvalue.node");

				for(j=0;j<numnodalvalues;j++){

					/*First create a nodalvalue object for that specific enum (nodalvalue_model_enum_s[j]):*/
					output_definitions->AddObject(new Nodalvalue(nodalvalue_name_s[j],StringToEnumx(nodalvalue_definitionstrings[j]),StringToEnumx(nodalvalue_modelstrings[j]),nodalvalue_node_s[j]-1)); //-1 because matlab to c indexing.
				}

				/*Free resources:*/
				for(j=0;j<numnodalvalues;j++){
					char* string=NULL;
					string = nodalvalue_name_s[j];    xDelete<char>(string);
				}
				xDelete<char*>(nodalvalue_name_s);
				xDelete<char*>(nodalvalue_modelstrings);
				xDelete<char*>(nodalvalue_definitionstrings);
				xDelete<int>(nodalvalue_node_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==MassconEnum){
				/*Deal with masscons: {{{*/

				/*masscon variables: */
				int          nummasscons;
				char**       masscon_name_s               = NULL;
				char**       masscon_definitionstring_s   = NULL;
				IssmDouble** masscon_levelset_s           = NULL;
				int*         masscon_levelset_M_s         = NULL;
				int*         masscon_levelset_N_s         = NULL;

				/*Fetch name and levelset, etc ... (see src/m/classes/masscon.m): */
				iomodel->FetchMultipleData(&masscon_name_s,&nummasscons,                                                "md.masscon.name");
				iomodel->FetchMultipleData(&masscon_definitionstring_s,&nummasscons,                                    "md.masscon.definitionstring");
				iomodel->FetchMultipleData(&masscon_levelset_s,&masscon_levelset_M_s,&masscon_levelset_N_s,&nummasscons,"md.masscon.levelset");

				for(j=0;j<nummasscons;j++){

					/*Create a masscon object: */
					output_definitions->AddObject(new Masscon(masscon_name_s[j],StringToEnumx(masscon_definitionstring_s[j]),masscon_levelset_s[j],masscon_levelset_M_s[j]));

				}

				/*Free resources:*/
				for(j=0;j<nummasscons;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = masscon_name_s[j];    xDelete<char>(string);
					string = masscon_definitionstring_s[j];    xDelete<char>(string);
					matrix = masscon_levelset_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(masscon_name_s);
				xDelete<IssmDouble*>(masscon_levelset_s);
				xDelete<int>(masscon_levelset_M_s);
				xDelete<int>(masscon_levelset_N_s);
				xDelete<char*>(masscon_definitionstring_s);

				/*}}}*/
			}
			else if (output_definition_enums[i]==MassconaxpbyEnum){
				/*Deal with masscon combinations: {{{*/

				/*masscon variables: */
				char**       masscon_name_s             = NULL;    
				char**		 masscon_definitionstring_s		= NULL;    
				char**       masscon_namex_s             = NULL;    
				char**       masscon_namey_s             = NULL;    
				IssmDouble*  masscon_alpha_s     = NULL;
				IssmDouble*  masscon_beta_s     = NULL;
				int          num;

				/*Fetch names and multiplicators, etc ... (see src/m/classes/masscon_axpby.m): */
				iomodel->FetchMultipleData(&masscon_name_s,&num,          "md.massconaxpby.name");
				iomodel->FetchMultipleData(&masscon_definitionstring_s,&num,"md.massconaxpby.definitionstring");
				iomodel->FetchMultipleData(&masscon_namex_s,&num,         "md.massconaxpby.namex");
				iomodel->FetchMultipleData(&masscon_namey_s,&num,         "md.massconaxpby.namey");
				iomodel->FetchMultipleData(&masscon_alpha_s,&num,         "md.massconaxpby.alpha");
				iomodel->FetchMultipleData(&masscon_beta_s,&num,          "md.massconaxpby.beta");
				for(j=0;j<num;j++){

					/*Create a masscon axpyb object: */
					output_definitions->AddObject(new Massconaxpby(masscon_name_s[j],StringToEnumx(masscon_definitionstring_s[j]),masscon_namex_s[j],masscon_namey_s[j],masscon_alpha_s[j],masscon_beta_s[j]));

				}

				/*Free resources:*/
				for(j=0;j<num;j++){
					char* string=NULL;
					string = masscon_definitionstring_s[j];    xDelete<char>(string);
					string = masscon_name_s[j];    xDelete<char>(string);
					string = masscon_namex_s[j];    xDelete<char>(string);
					string = masscon_namey_s[j];    xDelete<char>(string);
				}
				xDelete<char*>(masscon_definitionstring_s);
				xDelete<char*>(masscon_name_s);
				xDelete<char*>(masscon_namex_s);
				xDelete<char*>(masscon_namey_s);
				xDelete<IssmDouble>(masscon_alpha_s);
				xDelete<IssmDouble>(masscon_beta_s);
				/*}}}*/
			}
			else if (output_definition_enums[i]==RegionaloutputEnum){
				/*Deal with regional output: {{{*/

				/*regional output variables: */
				int          numout;
				char**       reg_name_s               = NULL;
				char**       reg_definitionstring_s   = NULL;
				char**       reg_outputnamestring_s   = NULL;
				IssmDouble** reg_mask_s               = NULL;
				int*         reg_mask_M_s             = NULL;
				int*         reg_mask_N_s             = NULL;

				/*Fetch name and mask, etc ... (see src/m/classes/regionaloutput.m): */
				iomodel->FetchMultipleData(&reg_name_s,&numout,                                                "md.regionaloutput.name");
				iomodel->FetchMultipleData(&reg_definitionstring_s,&numout,                                    "md.regionaloutput.definitionstring");
				iomodel->FetchMultipleData(&reg_outputnamestring_s,&numout,                                    "md.regionaloutput.outputnamestring");
				iomodel->FetchMultipleData(&reg_mask_s,&reg_mask_M_s,&reg_mask_N_s,&numout,                    "md.regionaloutput.mask");
				for(j=0;j<numout;j++){

					/*Create a regional output object: */
					output_definitions->AddObject(new Regionaloutput(reg_name_s[j],StringToEnumx(reg_definitionstring_s[j]),reg_outputnamestring_s[j],reg_mask_s[j],reg_mask_M_s[j]));

				}

				/*Free resources:*/
				for(j=0;j<numout;j++){
					char* string=NULL;
					IssmDouble* matrix = NULL;

					string = reg_name_s[j];    xDelete<char>(string);
					string = reg_definitionstring_s[j];    xDelete<char>(string);
					string = reg_outputnamestring_s[j];    xDelete<char>(string);
					matrix = reg_mask_s[j]; xDelete<IssmDouble>(matrix);
				}
				xDelete<char*>(reg_name_s);
				xDelete<IssmDouble*>(reg_mask_s);
				xDelete<int>(reg_mask_M_s);
				xDelete<int>(reg_mask_N_s);
				xDelete<char*>(reg_outputnamestring_s);
				xDelete<char*>(reg_definitionstring_s);
			/*}}}*/
			}
			else if (output_definition_enums[i]==NumberedcostfunctionEnum){
				/*Deal with numbered cost function: {{{*/

				/*Intermediary*/
				int          numout,numout2;
				char       **ncf_name_s             = NULL;
				char       **ncf_definitionstring_s = NULL;
				char       **cost_functions         = NULL;
				IssmDouble **cost_functions_weights = NULL;
				int*         cost_functions_weights_M = NULL;
				int*         cost_functions_weights_N = NULL;
				int          cost_function,domaintype;
				int          num_cost_functions;

				/*Process cost functions and convert from string to enums*/
				iomodel->FindConstant(&num_cost_functions,"md.numberedcostfunction.num_cost_functions");
				iomodel->FindConstant(&cost_functions,&num_cost_functions,"md.numberedcostfunction.cost_functions");
				if(num_cost_functions<1) _error_("No cost functions found");
				int* cost_function_enums=xNew<int>(num_cost_functions);
				for(int i=0;i<num_cost_functions;++i){
					cost_function_enums[i]=StringToEnumx(cost_functions[i]);
				}

				iomodel->FetchMultipleData(&ncf_name_s,&numout,"md.numberedcostfunction.name");
				iomodel->FetchMultipleData(&ncf_definitionstring_s,&numout2,"md.numberedcostfunction.definitionstring"); _assert_(numout2 == numout); 
				iomodel->FetchMultipleData(&cost_functions_weights,&cost_functions_weights_M,&cost_functions_weights_N,&numout2,"md.numberedcostfunction.cost_functions_coefficients");  _assert_(numout2 == numout);
				if(numout!=1) _error_("not implemented yet, check code here");

				/*Fetch Observations */
				iomodel->FindConstant(&domaintype,"md.mesh.domain_type");
				for(int i=0;i<num_cost_functions;i++){
					cost_function=cost_function_enums[i];
					if(     cost_function==ThicknessAbsMisfitEnum) iomodel->FetchDataToInput(inputs,elements,"md.numberedcostfunction.thickness_obs",InversionThicknessObsEnum);
					else if(cost_function==SurfaceAbsMisfitEnum)   iomodel->FetchDataToInput(inputs,elements,"md.numberedcostfunction.surface_obs",InversionSurfaceObsEnum);
					else if(cost_function==SurfaceAbsVelMisfitEnum
							|| cost_function==SurfaceRelVelMisfitEnum
							|| cost_function==SurfaceLogVelMisfitEnum
							|| cost_function==SurfaceLogVxVyMisfitEnum
							|| cost_function==SurfaceAverageVelMisfitEnum){
						iomodel->FetchDataToInput(inputs,elements,"md.numberedcostfunction.vx_obs",InversionVxObsEnum);
						if(domaintype!=Domain2DverticalEnum) iomodel->FetchDataToInput(inputs,elements,"md.numberedcostfunction.vy_obs",InversionVyObsEnum);
					}
				}

				for(j=0;j<numout;j++){

					/*Now, for this particular misfit object, make sure we plug into the elements: the observation, and the weights.*/
					for(Object* & object : elements->objects){
						Element* element=xDynamicCast<Element*>(object);
						element->DatasetInputCreate(cost_functions_weights[j],cost_functions_weights_M[j],cost_functions_weights_N[j],cost_function_enums,num_cost_functions,inputs,iomodel,InversionCostFunctionsCoefficientsEnum);
					}
					output_definitions->AddObject(new Numberedcostfunction(ncf_name_s[j],StringToEnumx(ncf_definitionstring_s[j]),num_cost_functions,cost_function_enums));
				}

				/*Free data: */
				iomodel->DeleteData(2,"md.numberedcostfunction.name","md.numberedcostfunction.definitionstring");
				xDelete<int>(cost_function_enums);
				for(int i=0;i<num_cost_functions;i++) xDelete<char>(cost_functions[i]);
				xDelete<char*>(cost_functions);

				/*Free resources:*/
				for(j=0;j<numout;j++){
					xDelete<char>(ncf_name_s[j]);
					xDelete<char>(ncf_definitionstring_s[j]);
					xDelete<IssmDouble>(cost_functions_weights[j]);
				}
				xDelete<char*>(ncf_name_s);
				xDelete<char*>(ncf_definitionstring_s);
				xDelete<int>(cost_functions_weights_M);
				xDelete<int>(cost_functions_weights_N);
				xDelete<IssmDouble*>(cost_functions_weights);

			/*}}}*/
			}
			else if (output_definition_enums[i]==RadarEnum){		
				/*Deal with radar: {{{*/
				int    numout;
				char **radar_name_s             = NULL;
				char **radar_definitionstring_s = NULL;
				int  **radar_ice_period_s       = NULL;

				/*Fetch name and definition, etc ... (see src/m/classes/radar.m): */
				iomodel->FetchMultipleData(&radar_definitionstring_s,&numout,"md.radar.definitionstring");
				iomodel->FetchMultipleData(&radar_name_s,&numout,"md.radar.name");
				if(numout>1) _error_("not suppored yet"); 
				/*Fetch necessary inputs for calculation*/
				//iomodel->FetchDataToInput(elements,"md.ice_period",RadarIcePeriodEnum);

				/*Add to output definitions*/
				output_definitions->AddObject(new Radar(radar_name_s[0],StringToEnumx(radar_definitionstring_s[0])));
				/*}}}*/ 
			}
		else _error_("output definition enum " << EnumToStringx(output_definition_enums[i]) << " not supported yet!");
		}		
	}
	parameters->AddObject(new DataSetParam(OutputdefinitionEnum,output_definitions));

	/*Free resources:*/
	delete output_definitions;
	xDelete<int>(output_definition_enums);

}
