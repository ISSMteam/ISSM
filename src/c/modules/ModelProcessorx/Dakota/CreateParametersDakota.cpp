/*!\file: CreateParametersDakota.cpp
 * \brief general driver for creating parameters dataset
 */ 

#include "../../../toolkits/toolkits.h"
#include "../../../classes/classes.h"
#include "../../../shared/shared.h"
#include "../../MeshPartitionx/MeshPartitionx.h"
#include "../ModelProcessorx.h"

void CreateParametersDakota(Parameters* parameters,IoModel* iomodel,char* rootpath){

	/*variable declarations*/
	int          i;
	char       **responsedescriptors    = NULL;
	int          numresponsedescriptors;
	char       **variabledescriptors    = NULL;
	int          numvariabledescriptors;
	char        *descriptor             = NULL;
	double      *dakota_parameter       = NULL;

	//qmu files
	char *qmuinname  = NULL;
	char *qmuerrname = NULL;
	char *qmuoutname = NULL;

	//descriptors:
	char tag[50];

	bool  dakota_analysis   = false;
	char *name              = NULL;
	int   numberofresponses;
	int   nrows,ncols;

	//variable partitions: 
	IssmDouble **array                      = NULL;
	int         *mdims_array                = NULL;
	int         *ndims_array                = NULL;
	int          num_partitions;
	int*         intarray = NULL;
	int M,N;

	//qmu statistics
	bool statistics    = false;
	int  numdirectories = 0;
	int  nfilesperdirectory = 0;

	/*recover parameters: */
	iomodel->FindConstant(&dakota_analysis,"md.qmu.isdakota");

	if(dakota_analysis){

		parameters->AddObject(iomodel->CopyConstantObject("md.qmu.output",QmuOutputEnum));

		iomodel->FindConstant(&name,"md.miscellaneous.name");
		iomodel->FindConstant(&numberofresponses,"md.qmu.numberofresponses");

		/*name of qmu input, error and output files*/
		qmuinname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.in")+1));
		sprintf(qmuinname,"%s%s%s",rootpath,name,".qmu.in");
		parameters->AddObject(new   StringParam(QmuInNameEnum,qmuinname));

		qmuoutname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.out")+1));
		sprintf(qmuoutname,"%s%s%s",rootpath,name,".qmu.out");
		parameters->AddObject(new   StringParam(QmuOutNameEnum,qmuoutname));

		qmuerrname=xNew<char>((strlen(rootpath)+strlen(name)+strlen(".qmu.err")+1));
		sprintf(qmuerrname,"%s%s%s",rootpath,name,".qmu.err");
		parameters->AddObject(new   StringParam(QmuErrNameEnum,qmuerrname));

		/*Fetch variable descriptors*/
		iomodel->FindConstant(&variabledescriptors,&numvariabledescriptors,"md.qmu.variabledescriptors");

		/*Fetch response descriptors*/
		iomodel->FindConstant(&responsedescriptors,&numresponsedescriptors,"md.qmu.responsedescriptors");

		/*Ok, we have all the response descriptors. Build a parameter with it: */
		parameters->AddObject(new StringArrayParam(QmuResponsedescriptorsEnum,responsedescriptors,numresponsedescriptors));

		/*Deal with statistics: */
		iomodel->FindConstant(&statistics,"md.qmu.statistics");
		parameters->AddObject(new BoolParam(QmuStatisticsEnum,statistics));
		if(statistics){
			iomodel->FindConstant(&numdirectories,"md.qmu.statistics.ndirectories");
			parameters->AddObject(new IntParam(QmuNdirectoriesEnum,numdirectories));

			iomodel->FindConstant(&nfilesperdirectory,"md.qmu.statistics.nfiles_per_directory");
			parameters->AddObject(new IntParam(QmuNfilesPerDirectoryEnum,nfilesperdirectory));
		}

		/*Load partitioning vectors specific to variables:*/
		iomodel->FetchData(&array,&mdims_array,&ndims_array,&num_partitions,"md.qmu.variablepartitions");
		parameters->AddObject(new DoubleMatArrayParam(QmuVariablePartitionsEnum,array,num_partitions,mdims_array,ndims_array));
		iomodel->FetchData(&intarray,&M,&N,"md.qmu.variablepartitions_npart");
		parameters->AddObject(new IntMatParam(QmuVariablePartitionsNpartEnum,intarray,M,N));
		xDelete<int>(intarray); iomodel->FetchData(&intarray,&M,&N,"md.qmu.variablepartitions_nt");
		parameters->AddObject(new IntMatParam(QmuVariablePartitionsNtEnum,intarray,M,N));

		/*free arrays: {{{*/
		for(i=0;i<num_partitions;i++){
			IssmDouble* matrix=array[i];
			xDelete<IssmDouble>(matrix);
		}
		xDelete<int>(mdims_array); 
		xDelete<int>(ndims_array);
		xDelete<IssmDouble*>(array);
		xDelete<int>(intarray);
		/*}}}*/

		/*Load partitioning vectors specific to responses:*/
		iomodel->FetchData(&array,&mdims_array,&ndims_array,&num_partitions,"md.qmu.responsepartitions");
		parameters->AddObject(new DoubleMatArrayParam(QmuResponsePartitionsEnum,array,num_partitions,mdims_array,ndims_array));
		iomodel->FetchData(&intarray,&M,&N,"md.qmu.responsepartitions_npart");
		parameters->AddObject(new IntMatParam(QmuResponsePartitionsNpartEnum,intarray,M,N));

		/*free arrays: {{{*/
		for(i=0;i<num_partitions;i++){
			IssmDouble* matrix=array[i];
			xDelete<IssmDouble>(matrix);
		}
		xDelete<int>(mdims_array); 
		xDelete<int>(ndims_array);
		xDelete<IssmDouble*>(array);
		xDelete<int>(intarray);
		/*}}}*/

		/*Deal with data needed because of qmu variables*/
		DataSet* dataset_variable_descriptors = new DataSet(QmuVariableDescriptorsEnum);
		for(i=0;i<numvariabledescriptors;i++){
			if (strncmp(variabledescriptors[i],"scaled_",7)==0){
				int code;

				/*Ok, we are dealing with a variable that is distributed over nodes or elements. Recover the name of the variable (ex: scaled_Thickness): */
				sscanf(variabledescriptors[i],"scaled_%s",tag);

				/*Get field name and input enum from tag*/
				char* fieldname  = NULL;
				int   param_enum = -1;
				FieldAndEnumFromCode(&param_enum,&fieldname,tag);

				iomodel->SetFilePointerToData(&code,NULL,fieldname);
				if(code==8) dataset_variable_descriptors->AddObject(new DoubleParam(param_enum,8)); //skip MatArray inputs, as we don't know which input will be scaled yet!
				else{ 
					/*recover more classic data, arrays and scalar mainly:*/
					iomodel->FetchData(&dakota_parameter,&nrows,&ncols,fieldname);
					if(nrows==iomodel->numberofvertices || nrows==iomodel->numberofelements){
						dataset_variable_descriptors->AddObject(new DoubleMatParam(param_enum,dakota_parameter,nrows,ncols));
					}
					else{
						dataset_variable_descriptors->AddObject(new DoubleTransientMatParam(param_enum,dakota_parameter,nrows,ncols));
					}
					xDelete<double>(dakota_parameter);
				}
				xDelete<char>(fieldname);
			}
		}
		parameters->AddObject(new DataSetParam(QmuVariableDescriptorsEnum,dataset_variable_descriptors));
		delete dataset_variable_descriptors;

		/*clean-up {{{*/
		for(i=0;i<numresponsedescriptors;i++){
			descriptor=responsedescriptors[i];
			xDelete<char>(descriptor);
		}
		xDelete<char*>(responsedescriptors);
		for(i=0;i<numvariabledescriptors;i++){
			descriptor=variabledescriptors[i];
			xDelete<char>(descriptor);
		}
		xDelete<char*>(variabledescriptors);
		xDelete<char>(qmuinname);
		xDelete<char>(qmuerrname);
		xDelete<char>(qmuoutname);
		xDelete<char>(name);
		/*}}}*/
	}

}
