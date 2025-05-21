/*!\file FemModel.cpp
 * \brief: implementation of the FemModel object
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <math.h>
#include <float.h>/*  DBL_EPSILON  */
#include "../cores/cores.h"
#include "../shared/io/io.h"
#include "./classes.h"
#include "./Inputs/TriaInput.h"
#include "./modules/modules.h"
#include "../shared/Enum/Enum.h"
#include "../analyses/analyses.h"
#include "./Inputs/DatasetInput.h"
#include "./Inputs/ElementInput.h"
#include "./Inputs/TransientInput.h"

#include "../toolkits/codipack/CoDiPackGlobal.h"

#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
#include <TPZRefPatternDataBase.h>
#endif

/*module includes: {{{*/
#include "../modules/ModelProcessorx/ModelProcessorx.h"
#include "../modules/SpcNodesx/SpcNodesx.h"
#include "../modules/ConfigureObjectsx/ConfigureObjectsx.h"
#include "../modules/ParseToolkitsOptionsx/ParseToolkitsOptionsx.h"
#include "../modules/GetVectorFromInputsx/GetVectorFromInputsx.h"
#include "../modules/InputUpdateFromVectorx/InputUpdateFromVectorx.h"
#include "../modules/NodesDofx/NodesDofx.h"
#include "../modules/SurfaceAbsVelMisfitx/SurfaceAbsVelMisfitx.h"
#include "../modules/SurfaceRelVelMisfitx/SurfaceRelVelMisfitx.h"
#include "../modules/SurfaceLogVelMisfitx/SurfaceLogVelMisfitx.h"
#include "../modules/SurfaceLogVxVyMisfitx/SurfaceLogVxVyMisfitx.h"
#include "../modules/SurfaceAverageVelMisfitx/SurfaceAverageVelMisfitx.h"
#include "../modules/ThicknessAbsMisfitx/ThicknessAbsMisfitx.h"
#include "../modules/ThicknessAlongGradientx/ThicknessAlongGradientx.h"
#include "../modules/ThicknessAcrossGradientx/ThicknessAcrossGradientx.h"
#include "../modules/RheologyBbarAbsGradientx/RheologyBbarAbsGradientx.h"
#include "../modules/DragCoefficientAbsGradientx/DragCoefficientAbsGradientx.h"
#include "../modules/NodalValuex/NodalValuex.h"
#include "../modules/AverageOntoPartitionx/AverageOntoPartitionx.h"
/*}}}*/

/*Object constructors and destructor*/
FemModel::FemModel(void){ /*{{{*/
	/*do nothing:*/
} /*}}}*/
FemModel::FemModel(int argc,char** argv,ISSM_MPI_Comm incomm,bool trace){/*{{{*/

	/*configuration: */
	int  solution_type,amrtype,amr_frequency;
	int  ierr;

	/*File names*/
	char *lockfilename   = NULL;
	char *binfilename    = NULL;
	char *outbinfilename = NULL;
	char *petscfilename  = NULL;
	char *restartfilename  = NULL;
	char *rootpath       = NULL;
	char *modelname       = NULL;

	/*First things first, store the communicator, and set it as a global variable: */
	IssmComm::SetComm(incomm);

	/*Now, initialize PETSC: */
	#ifdef _HAVE_PETSC_
	PETSC_COMM_WORLD=incomm;
	ierr=PetscInitialize(&argc,&argv,(char*)0,"");  if(ierr) _error_("Could not initialize Petsc");
	#endif

	/*Start profiler: */
	this->profiler=new Profiler();
	profiler->Start(TOTAL);

	/*From command line arguments, retrieve different filenames needed to create the FemModel: */
	ProcessArguments(&solution_type,&binfilename,&outbinfilename,&petscfilename,&lockfilename,&restartfilename,&rootpath,&modelname,argc,argv);

	/*Create femmodel from input files: */
	profiler->Start(MPROCESSOR);
	this->InitFromFiles(rootpath,binfilename,outbinfilename,petscfilename,lockfilename,restartfilename, modelname, solution_type,trace,NULL);
	profiler->Stop(MPROCESSOR);

	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));

   /*AMR stuff*/
	this->parameters->FindParam(&amr_frequency,TransientAmrFrequencyEnum);
	this->parameters->FindParam(&amr_frequency,TransientAmrFrequencyEnum);
	#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
	this->amr = NULL;
	#endif
	#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
	this->amrbamg = NULL;
	#endif
	#if !defined(_HAVE_AD_)
	if(amr_frequency && solution_type==TransientSolutionEnum){
		/*Verifications. AMR supports SSA, P1 and horizontal 2D domain*/
		bool isSSA;
		int domaintype,element_type;
		this->analysis_counter=-1;
		this->parameters->FindParam(&isSSA,FlowequationIsSSAEnum);
		this->parameters->FindParam(&domaintype,DomainTypeEnum);
		for(int i=0;i<this->nummodels;i++) {
			if(this->analysis_type_list[i]==StressbalanceAnalysisEnum){
				analysis_counter=i;
				break;
			}
		}
		if(analysis_counter==-1) _error_("Could not find alias for analysis_type StressbalanceAnalysisEnum in list of FemModel analyses\n");
		for(Object* & object : this->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element_type		= element->element_type_list[analysis_counter];
			if(element_type!=P1Enum) _error_("Element type "<<EnumToStringx(element_type)<<" not supported with AMR yet!\n");
		}
		if(!isSSA) _error_("Flow equation not supported with AMR yet!\n ");
		if(domaintype!=Domain2DhorizontalEnum) _error_("Domain "<<EnumToStringx(domaintype)<<" not supported with AMR yet!\n");

		this->parameters->FindParam(&amrtype,AmrTypeEnum);
		switch(amrtype){

			#if defined(_HAVE_NEOPZ_)
			case AmrNeopzEnum: this->InitializeAdaptiveRefinementNeopz(); break;
			#endif

			#if defined(_HAVE_BAMG_)
			case AmrBamgEnum: this->InitializeAdaptiveRefinementBamg(); break;
			#endif

			default: _error_("not implemented yet");
		}
	}
	#endif

	/*Free resources */
	xDelete<char>(lockfilename);
	xDelete<char>(binfilename);
	xDelete<char>(outbinfilename);
	xDelete<char>(petscfilename);
	xDelete<char>(restartfilename);
	xDelete<char>(modelname);
	xDelete<char>(rootpath);

}
/*}}}*/
FemModel::FemModel(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, char* modelname, ISSM_MPI_Comm incomm, int solution_type,IssmPDouble* X){ /*{{{*/

	bool traceon=true;
	this->profiler=NULL; /*avoid leak, as we are not using the profiler ever in ad control run. */

	/*Store the communicator, but do not set it as a global variable, as this has already
	 * been done by the FemModel that called this copy constructor: */
	IssmComm::SetComm(incomm);

	/*Create femmodel from input files, with trace activated: */
	profiler->Start(MPROCESSOR);
	this->InitFromFiles(rootpath,inputfilename,outputfilename,toolkitsfilename,lockfilename,restartfilename, modelname,solution_type,traceon,X);
	profiler->Stop(MPROCESSOR);

	#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
	this->amr = NULL;
	#endif
	#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
	this->amrbamg = NULL;
	#endif

	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));
}
/*}}}*/
FemModel::~FemModel(){/*{{{*/

	/*Intermediary*/
	FILE *output_fid;
	char *outbinfilename = NULL;
	char *lockfilename   = NULL;

	#ifndef _HAVE_JAVASCRIPT_
	if(this->parameters->Exist(OutputFileNameEnum)) this->parameters->FindParam(&outbinfilename,OutputFileNameEnum);
	if(this->parameters->Exist(LockFileNameEnum)) this->parameters->FindParam(&lockfilename,LockFileNameEnum);
	#endif

	/*Delete all the datasets: */
	if(analysis_type_list)xDelete<int>(analysis_type_list);
	if(outbinfilename)xDelete<char>(outbinfilename);
	if(lockfilename)xDelete<char>(lockfilename);
	if(elements)delete elements;
	if(vertices)delete vertices;
	if(this->constraints_list && this->nummodels){
		for(int i=0;i<this->nummodels;i++) delete this->constraints_list[i];
		xDelete<Constraints*>(constraints_list);
	}
	if(this->loads_list && this->nummodels){
		for(int i=0;i<this->nummodels;i++) delete this->loads_list[i];
		xDelete<Loads*>(loads_list);
	}
	if(this->nodes_list && this->nummodels){
		for(int i=0;i<this->nummodels;i++) delete this->nodes_list[i];
		xDelete<Nodes*>(nodes_list);
	}
	if(materials)delete materials;
	if(parameters)delete parameters;
	if(inputs)delete inputs;
	if(results)delete results;

	#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
	if(amr)delete amr;
	#endif

	#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
	if(amrbamg)delete amrbamg;
	#endif

	/*Now delete: */
	if(profiler)delete profiler;
}/*}}}*/

/*Object management*/
int FemModel::AnalysisIndex(int analysis_enum){/*{{{*/

	/*Checks in debugging mode*/
	_assert_(this->analysis_type_list);

	/*Find analysis in list*/
	for(int i=0;i<this->nummodels;i++){
		if(this->analysis_type_list[i]==analysis_enum){
			return i;
			break;
		}
	}

	/*If you reach this point, analysis has not been found*/
	_error_("Could not find index of analysis " << EnumToStringx(analysis_enum) << " in list of FemModel analyses");

}/*}}}*/
void FemModel::CheckPoint(void){/*{{{*/

	FILE* restartfid=NULL;
	char* restartfilename = NULL;
	int   femmodel_size;
	char* femmodel_buffer=NULL;
	char* femmodel_buffer_ini=NULL;

	/*First, recover the name of the restart file: */
	parameters->FindParam(&restartfilename,RestartFileNameEnum);

	/*Open file for writing: */
	restartfid=pfopen(restartfilename,"wb");

	/*Initialize: */
	femmodel_size=this->Size();
	_assert_(femmodel_size);

	/*Create buffer to hold marshalled femmodel: */
	femmodel_buffer=xNew<char>(femmodel_size);

	/*Keep track of initial position of femmodel_buffer: */
	femmodel_buffer_ini=femmodel_buffer;

	/*Marshall:*/
   WriteCheckpointFunctor* marshallhandle = new WriteCheckpointFunctor(&femmodel_buffer);
   this->Marshall(marshallhandle);
	delete marshallhandle;

	/*Reset position of buffer: */
	femmodel_buffer=femmodel_buffer_ini;

	/*write buffer: */
	fwrite(femmodel_buffer,femmodel_size,sizeof(char),restartfid);

	/*Done, close file :*/
	pfclose(restartfid,restartfilename);

	/*Free resources: */
	xDelete<char>(femmodel_buffer);
	xDelete<char>(restartfilename);
}
/*}}}*/
void FemModel::CheckPointAD(int step){/*{{{*/

	/*Get rank*/
	int my_rank = IssmComm::GetRank();

	/*Get string sizes*/
	int rank_length = (my_rank == 0 ? 1 : int(log10(static_cast<double>(my_rank))+1));
	int step_length = (step    == 0 ? 1 : int(log10(static_cast<double>(step))   +1));

	/*Create restart file*/
	int restartfilename_len = strlen("AD_step_")+step_length+strlen("_rank_")+rank_length+strlen(".ckpt")+1;
	char* restartfilename = xNew<char>(restartfilename_len);
	snprintf(restartfilename, restartfilename_len, "%s%i%s%i%s","AD_step_",step,"_rank_",my_rank,".ckpt");
	this->parameters->AddObject(new StringParam(RestartFileNameEnum,restartfilename));

	/*Write files*/
	this->CheckPoint();

	/*Clean up and return*/
	xDelete<char>(restartfilename);

}/*}}}*/
void FemModel::CleanUp(void){/*{{{*/

	/*Intermediary*/
	char *lockfilename   = NULL;
	bool  waitonlock     = false;

	/*Write lock file if requested: */
	this->parameters->FindParam(&waitonlock,SettingsWaitonlockEnum);
	this->parameters->FindParam(&lockfilename,LockFileNameEnum);
	if(waitonlock){
		_printf0_("write lock file:\n");
		WriteLockFile(lockfilename);
	}

	/*Before we delete the profiler, report statistics for this run: */
	profiler->Stop(TOTAL);  //final tagging

	_printf0_("\n");
	_printf0_("   "<<setw(40)<<left<<"FemModel initialization elapsed time:"<<setw(7)<<profiler->TotalTime(MPROCESSOR) << "\n");
	/*Total times*/
	_printf0_("   "<<setw(40)<<left<<"Total Core solution elapsed time:"<<setw(7)<<profiler->TotalTime(CORE) << "\n");

	/*Linear solver only*/
	_printf0_("   "<<setw(40)<<left<<"Linear solver elapsed time:"<<setw(7)<<profiler->TotalTime(SOLVER) << " ("<<setprecision(2)<<profiler->TotalTime(SOLVER)/profiler->TotalTime(CORE)*100.<<"%)\n");
	_printf0_("\n");
	_printf0_("   Total elapsed time: "
				<<profiler->TotalTimeModHour(TOTAL)<<" hrs "
				<<profiler->TotalTimeModMin(TOTAL)<<" min "
				<<profiler->TotalTimeModSec(TOTAL)<<" sec"
				);
	_printf0_("\n");

	/*Finalize PETSC for this model: */
	#ifdef _HAVE_PETSC_
	//_printf0_("closing PETSc\n");
	PetscFinalize();
	#endif

	/*Cleanup toolkit*/
	ToolkitOptions::Delete();

	/*Clean up*/
	xDelete<char>(lockfilename);
} /*}}}*/
FemModel* FemModel::copy(void){/*{{{*/

	FemModel* output=NULL;
	int       i;
	int       analysis_type;

	output=new FemModel(*this); //Use default copy constructor.

	output->nummodels = this->nummodels;
	output->solution_type = this->solution_type;
	output->analysis_counter = this->analysis_counter;

	/*Now, deep copy arrays: */
	output->analysis_type_list=xNew<int>(nummodels);
	xMemCpy<int>(output->analysis_type_list,this->analysis_type_list,this->nummodels);

	/*Analysis dependent arrays*/
	output->constraints_list=xNew<Constraints*>(this->nummodels);
	output->loads_list=xNew<Loads*>(this->nummodels);
	output->nodes_list=xNew<Nodes*>(this->nummodels);

	output->profiler=static_cast<Profiler*>(this->profiler->copy());

	output->materials=static_cast<Materials*>(this->materials->Copy());
	output->parameters=static_cast<Parameters*>(this->parameters->Copy());
	output->inputs=static_cast<Inputs*>(this->inputs->Copy());
	output->results=static_cast<Results*>(this->results->Copy());
	output->vertices=static_cast<Vertices*>(this->vertices->Copy());
	output->elements=static_cast<Elements*>(this->elements->Copy());

	/*reset hooks for elements, loads and nodes: */
	output->elements->ResetHooks();
	output->loads->ResetHooks();
	output->materials->ResetHooks();

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	for(i=0;i<nummodels;i++){
		output->constraints_list[i] = static_cast<Constraints*>(this->constraints_list[i]->Copy());
		output->loads_list[i] = static_cast<Loads*>(this->loads_list[i]->Copy());
		output->nodes_list[i] = static_cast<Nodes*>(this->nodes_list[i]->Copy());
		analysis_type=output->analysis_type_list[i];
		output->SetCurrentConfiguration(analysis_type);
		SpcNodesx(output->nodes_list[i],output->constraints_list[i],output->parameters);
		NodesDofx(output->nodes_list[i],output->parameters);
		ConfigureObjectsx(output->elements,output->loads_list[i],output->nodes_list[i],output->vertices,output->materials,output->parameters,output->inputs);
	}

	/*AMR, no copy for now*/
	#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
	this->amr = NULL;
	#endif
	#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
	this->amrbamg = NULL;
	#endif

	/*Reset current configuration: */
	analysis_type=output->analysis_type_list[analysis_counter];
	output->SetCurrentConfiguration(analysis_type);

	return output;
}
/*}}}*/
void FemModel::Echo(void){/*{{{*/

	_printf_("FemModel echo: \n");
	_printf_("   number of fem models: " << nummodels << "\n");
	_printf_("   analysis_type_list: \n");
	for(int i=0;i<nummodels;i++)_printf_("     " << i << ": " << EnumToStringx(analysis_type_list[i]) << "\n");
	_printf_("   current analysis_type: \n");
	_printf_("     " << analysis_counter << ": " << EnumToStringx(analysis_type_list[analysis_counter]) << "\n");

}
/*}}}*/
void FemModel::InitFromFiles(char* rootpath, char* inputfilename, char* outputfilename, char* toolkitsfilename, char* lockfilename, char* restartfilename, char* modelname, const int in_solution_type,bool trace,IssmPDouble* X){/*{{{*/

	/*intermediary*/
	FILE *IOMODEL            = NULL;
	FILE *toolkitsoptionsfid = NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Open input file descriptor on cpu 0: */
	if(my_rank==0) IOMODEL = pfopen0(inputfilename ,"rb");

	/*Open toolkits file: */
	toolkitsoptionsfid=pfopen(toolkitsfilename,"r");

	/*Now, go create FemModel:*/
	this->InitFromFids(rootpath,IOMODEL,toolkitsoptionsfid,in_solution_type,trace,X);

	/*Close input file and toolkits file descriptors: */
	if(my_rank==0) pfclose(IOMODEL,inputfilename);
	pfclose(toolkitsoptionsfid,toolkitsfilename);

	/*Now save all of these file names into parameters, you never know when you might need them: */
	this->parameters->AddObject(new StringParam(ToolkitsFileNameEnum,toolkitsfilename));
	this->parameters->AddObject(new StringParam(ModelnameEnum,modelname));
	this->parameters->AddObject(new StringParam(RootPathEnum,rootpath));
	this->parameters->AddObject(new StringParam(InputFileNameEnum,inputfilename));
	this->parameters->AddObject(new StringParam(OutputFileNameEnum,outputfilename));
	this->parameters->AddObject(new StringParam(LockFileNameEnum,lockfilename));
	this->parameters->AddObject(new StringParam(RestartFileNameEnum,restartfilename));

}/*}}}*/
void FemModel::InitFromFids(char* rootpath, FILE* IOMODEL, FILE* toolkitsoptionsfid, int in_solution_type, bool trace, IssmPDouble* X){/*{{{*/

	/*Initialize internal data: */
	this->solution_type    = in_solution_type;
	this->analysis_counter = -1;
	this->results          = new Results(); //not initialized by CreateDataSets

	/*create IoModel */
	IoModel* iomodel = new IoModel(IOMODEL,in_solution_type,trace,X);

	/*Figure out what analyses are activated for this solution*/
	SolutionAnalysesList(&this->analysis_type_list,&this->nummodels,iomodel,this->solution_type);

	/*create datasets for all analyses*/
	ModelProcessorx(&this->elements,&this->nodes_list,&this->vertices,&this->materials,&this->constraints_list,&this->loads_list,&this->parameters,&this->inputs,iomodel,toolkitsoptionsfid,rootpath,this->solution_type,this->nummodels,this->analysis_type_list);

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	for(int i=0;i<nummodels;i++){

		if(VerboseMProcessor()) _printf0_("   Processing finite element model of analysis " << EnumToStringx(analysis_type_list[i]) << ":\n");
		this->SetCurrentConfiguration(analysis_type_list[i]);

		if(VerboseMProcessor()) _printf0_("      configuring element and loads\n");
		ConfigureObjectsx(this->elements,this->loads,this->nodes,this->vertices,this->materials,this->parameters,this->inputs);

		if(i==0){
			if(VerboseMProcessor()) _printf0_("      detecting active vertices\n");
			GetMaskOfIceVerticesLSMx0(this);
		}

		if(VerboseMProcessor()) _printf0_("      resolving node constraints\n");
		SpcNodesx(nodes,this->constraints,parameters);

		if(VerboseMProcessor()) _printf0_("      creating nodal degrees of freedom\n");
		NodesDofx(nodes,parameters);
	}

	/*Clean up*/
	delete iomodel;
}/*}}}*/
void FemModel::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	/*Allocate new fields if necessary*/
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		delete this->materials;
		delete this->parameters;
		delete this->inputs;
		if(this->constraints_list && this->nummodels){
			for(int i=0;i<this->nummodels;i++) delete this->constraints_list[i];
			xDelete<Constraints*>(constraints_list);
		}
		if(this->loads_list && this->nummodels){
			for(int i=0;i<this->nummodels;i++) delete this->loads_list[i];
			xDelete<Loads*>(loads_list);
		}
		if(this->nodes_list && this->nummodels){
			for(int i=0;i<this->nummodels;i++) delete this->nodes_list[i];
			xDelete<Nodes*>(nodes_list);
		}
		delete this->results;
		delete this->vertices;
		delete this->elements;
		xDelete<int>(this->analysis_type_list);

		this->materials   = new Materials();
		this->parameters  = new Parameters();
		this->inputs      = new Inputs();
		this->results     = new Results();
		this->vertices    = new Vertices();
		this->elements    = new Elements();
	}

	int obj_enum = FemModelEnum;
	marshallhandle->call(obj_enum);

	marshallhandle->call(this->solution_type);
	marshallhandle->call(this->analysis_counter);
	marshallhandle->call(this->nummodels);
	marshallhandle->call(this->analysis_type_list,nummodels);

	this->materials->Marshall(marshallhandle);
	this->parameters->Marshall(marshallhandle);
	this->inputs->Marshall(marshallhandle);
	this->results->Marshall(marshallhandle);
	this->vertices->Marshall(marshallhandle);
	this->elements->Marshall(marshallhandle);

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		this->constraints_list = xNew<Constraints*>(this->nummodels);
		for(int i=0;i<nummodels;i++) this->constraints_list[i] = new Constraints();
		this->loads_list = xNew<Loads*>(this->nummodels);
		for(int i=0;i<nummodels;i++) this->loads_list[i] = new Loads();
		this->nodes_list = xNew<Nodes*>(this->nummodels);
		for(int i=0;i<nummodels;i++) this->nodes_list[i] = new Nodes();
	}

	for(int i=0;i<nummodels;i++){
		this->constraints_list[i]->Marshall(marshallhandle);
		this->loads_list[i]->Marshall(marshallhandle);
		this->nodes_list[i]->Marshall(marshallhandle);
	}
	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD){
		/*reset hooks for elements, loads and nodes:*/
		this->elements->ResetHooks();
		this->materials->ResetHooks();

		/*do the post-processing of the datasets to get an FemModel that can actually run analyses:*/
		for(int i=0;i<nummodels;i++){
			this->loads_list[i]->ResetHooks();
			int analysis_type=this->analysis_type_list[i];
			SetCurrentConfiguration(analysis_type);
			SpcNodesx(this->nodes_list[i],this->constraints_list[i],this->parameters);
			NodesDofx(this->nodes_list[i],this->parameters);
			ConfigureObjectsx(this->elements,this->loads_list[i],this->nodes_list[i],this->vertices,this->materials,this->parameters,this->inputs);
		}

		/*Reset current configuration*/
		SetCurrentConfiguration(this->analysis_type_list[this->analysis_counter]);
	}
}
/*}}}*/
void FemModel::Restart(int verboselevel){ /*{{{*/

	FILE *restartfid          = NULL;
	char *restartfilename     = NULL;
	int   femmodel_size       = 0;
	int   fread_return        = 0;
	char *femmodel_buffer     = NULL;
	char *femmodel_buffer_ini = NULL;

	/*First, recover the name of the restart file: */
	parameters->FindParam(&restartfilename,RestartFileNameEnum);

	/*Now, figure out whether this file actually exists!: */
	restartfid=pfopen(restartfilename,"r",false);

	if(restartfid==NULL){
		xDelete<char>(restartfilename);
		return; //could not find the file, so no restart possible.
	}

	/*Print banner*/
   if(verboselevel>1){
      _printf0_("                                                                    \n");
      _printf0_("====================================================================\n");
      _printf0_(" RESTART DETECTED: "<<restartfilename<<                            "\n");
      _printf0_("====================================================================\n");
      _printf0_("                                                                    \n");
   }
   else if(verboselevel==1){
      _printf0_("    == restarting from "<<restartfilename<<"\n");
   }
   else{
      /*Do not print anything*/
   }

	/*Figure out size of buffer to be read: */
	fseek(restartfid, 0L, SEEK_END);
	femmodel_size = ftell(restartfid);
	fseek(restartfid, 0L, SEEK_SET);

	/*Allocate buffer: */
	femmodel_buffer=xNew<char>(femmodel_size);

	/*Read buffer from file: */
	fread_return=fread(femmodel_buffer,femmodel_size,sizeof(char),restartfid); if(fread_return!=1)_error_("error reading the buffer from marshalled file!");
	femmodel_buffer_ini=femmodel_buffer; //keep track of the initial position, so as to free later.

	/*Create new FemModel by demarshalling the buffer: */
   LoadCheckpointFunctor* marshallhandle = new LoadCheckpointFunctor(&femmodel_buffer);
   this->Marshall(marshallhandle);
	delete marshallhandle;

	/*Reset position of buffer: */
	femmodel_buffer=femmodel_buffer_ini;

	/*Done, close file :*/
	pfclose(restartfid,restartfilename);

	/*Free resources: */
	xDelete<char>(restartfilename);
	xDelete<char>(femmodel_buffer);
}/*}}}*/
void FemModel::RestartAD(int step){ /*{{{*/

	/*Get rank*/
	int my_rank = IssmComm::GetRank();

	/*Get string sizes*/
	int rank_length = (my_rank == 0 ? 1 : int(log10(static_cast<double>(my_rank))+1));
	int step_length = (step    == 0 ? 1 : int(log10(static_cast<double>(step))   +1));

	/*Create restart file*/
	int   restartfilename_len = strlen("AD_step_")+step_length+strlen("_rank_")+rank_length+strlen(".ckpt")+1;
	char* restartfilename  = xNew<char>(restartfilename_len);
	snprintf(restartfilename, restartfilename_len,"%s%i%s%i%s","AD_step_",step,"_rank_",my_rank,".ckpt");
	this->parameters->AddObject(new StringParam(RestartFileNameEnum,restartfilename));

	/*Read files*/
	this->Restart(1);

	/*Delete checkpoint file to save disk space*/
	/*_printf0_("    == deleting  file  "<<restartfilename<<"\n");*/
	std::remove(restartfilename);

	/*Clean up and return*/
	xDelete<char>(restartfilename);
}/*}}}*/
void FemModel::SetCurrentConfiguration(int configuration_type,int analysis_type){/*{{{*/

	/*Use configuration_type to setup the analysis counter, the configurations of objects etc ... but use
	 * analysis_type to drive the element numerics. This allows for use of 1 configuration_type for several
	 * analyses. For example: do a SurfaceSlopeX, SurfaceSlopeY, BedSlopeX and BedSlopeY analysis using the
	 * Slope configuration.*/
	int index = AnalysisIndex(configuration_type);

	/*If we already have the right analysis, return*/
	//if(this->analysis_counter==index) return;
	this->analysis_counter=index;

	/*Now, plug analysis_counter and analysis_type inside the parameters: */
	this->parameters->SetParam(analysis_counter,AnalysisCounterEnum);
	this->parameters->SetParam(analysis_type,AnalysisTypeEnum);
	this->parameters->SetParam(configuration_type,ConfigurationTypeEnum);

	/*configure elements, loads and nodes, for this new analysis: */
	this->loads = this->loads_list[this->analysis_counter];
	this->constraints = this->constraints_list[this->analysis_counter];
	this->nodes = this->nodes_list[this->analysis_counter];
	this->loads->SetCurrentConfiguration(elements, loads, nodes,vertices, materials,parameters);
	this->elements->SetCurrentConfiguration(elements,loads, nodes,vertices, materials,parameters);

	/*take care of toolkits options, that depend on this analysis type (present only after model processor)*/
	if(this->parameters->Exist(ToolkitsOptionsStringsEnum)){
		ToolkitsOptionsFromAnalysis(this->parameters,analysis_type);
		if(VerboseSolver()) _printf0_("      toolkits Options set for analysis: " << EnumToStringx(analysis_type) << "\n");
	}

}/*}}}*/
void FemModel::SetCurrentConfiguration(int configuration_type){/*{{{*/
	this->SetCurrentConfiguration(configuration_type,configuration_type);
}
/*}}}*/
int  FemModel::Size(){ /*{{{*/

	SizeCheckpointFunctor* marshallhandle = new SizeCheckpointFunctor();
	this->Marshall(marshallhandle);
	int femmodel_size = marshallhandle->MarshalledSize();

	/*Cleanup and return*/
	delete marshallhandle;
	return femmodel_size;
}
/*}}}*/
void FemModel::SolutionAnalysesList(int** panalyses,int* pnumanalyses,IoModel* iomodel,int solutiontype){/*{{{*/

	/*output: */
	int  numanalyses = 0;
	int *analyses    = NULL;

	/*Intermediaries*/
	const int MAXANALYSES = 30;
	int   analyses_temp[MAXANALYSES];

	/*Analyses lists*/
	switch(solutiontype){

		case StressbalanceSolutionEnum:{
			bool isSIA,isFS;
			int  fe_FS;
			iomodel->FindConstant(&fe_FS,"md.flowequation.fe_FS");
			iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
			iomodel->FindConstant(&isFS,"md.flowequation.isFS");
			analyses_temp[numanalyses++]=StressbalanceAnalysisEnum;
			analyses_temp[numanalyses++]=StressbalanceVerticalAnalysisEnum;
			if(isSIA){
				analyses_temp[numanalyses++]=StressbalanceSIAAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;
			analyses_temp[numanalyses++]=DepthAverageAnalysisEnum;
			if(fe_FS==LATaylorHoodEnum || fe_FS==LACrouzeixRaviartEnum){
				analyses_temp[numanalyses++]=UzawaPressureAnalysisEnum;
			}
			}
			break;

		case SteadystateSolutionEnum:{
			bool isSIA,isenthalpy;
			iomodel->FindConstant(&isSIA,"md.flowequation.isSIA");
			iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
			analyses_temp[numanalyses++]=StressbalanceAnalysisEnum;
			analyses_temp[numanalyses++]=StressbalanceVerticalAnalysisEnum;
			if(isSIA){
				analyses_temp[numanalyses++]=StressbalanceSIAAnalysisEnum;
			}
			if(isenthalpy){
				analyses_temp[numanalyses++]=EnthalpyAnalysisEnum;
			}
			else{
				analyses_temp[numanalyses++]=ThermalAnalysisEnum;
				analyses_temp[numanalyses++]=MeltingAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			}
			break;

		case ThermalSolutionEnum:{
			bool isenthalpy;
			iomodel->FindConstant(&isenthalpy,"md.thermal.isenthalpy");
			if(isenthalpy){
				analyses_temp[numanalyses++]=EnthalpyAnalysisEnum;
			}
			else{
				analyses_temp[numanalyses++]=ThermalAnalysisEnum;
				analyses_temp[numanalyses++]=MeltingAnalysisEnum;
			}
			}
			break;

		case HydrologySolutionEnum:{
			int hydrology_model;
			iomodel->FindConstant(&hydrology_model,"md.hydrology.model");
			if(hydrology_model==HydrologyshreveEnum){
				analyses_temp[numanalyses++]=HydrologyShreveAnalysisEnum;
			}
			else if(hydrology_model==HydrologyGlaDSEnum){
				analyses_temp[numanalyses++]=HydrologyGlaDSAnalysisEnum;
			}
			if(hydrology_model==HydrologyshaktiEnum){
				analyses_temp[numanalyses++]=HydrologyShaktiAnalysisEnum;
			}
			if(hydrology_model==HydrologypismEnum){
				analyses_temp[numanalyses++]=HydrologyPismAnalysisEnum;
			}
			if(hydrology_model==HydrologyTwsEnum){
				analyses_temp[numanalyses++]=HydrologyTwsAnalysisEnum;
			}
			if(hydrology_model==HydrologydcEnum){
				analyses_temp[numanalyses++]=HydrologyDCInefficientAnalysisEnum;
				analyses_temp[numanalyses++]=HydrologyDCEfficientAnalysisEnum;
				analyses_temp[numanalyses++]=L2ProjectionEPLAnalysisEnum;
				analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			}
			if(hydrology_model==HydrologyarmapwEnum){
				analyses_temp[numanalyses++]=HydrologyArmapwAnalysisEnum;
			}
		}
			break;

		case MasstransportSolutionEnum:
			analyses_temp[numanalyses++]=DepthAverageAnalysisEnum;
			analyses_temp[numanalyses++]=SmbAnalysisEnum;
			analyses_temp[numanalyses++]=MasstransportAnalysisEnum;
			analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;

			break;

		case BalancethicknessSolutionEnum:
			analyses_temp[numanalyses++]=BalancethicknessAnalysisEnum;
			break;

		case Balancethickness2SolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			analyses_temp[numanalyses++]=SmoothAnalysisEnum;
			analyses_temp[numanalyses++]=Balancethickness2AnalysisEnum;
			break;

		case BalancethicknessSoftSolutionEnum:
			analyses_temp[numanalyses++]=BalancethicknessAnalysisEnum;
			break;

		case BalancevelocitySolutionEnum:
			analyses_temp[numanalyses++]=BalancevelocityAnalysisEnum;
			analyses_temp[numanalyses++]=SmoothAnalysisEnum;
			break;

		case SurfaceSlopeSolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			break;

		case BedSlopeSolutionEnum:
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			break;

		case LoveSolutionEnum:
			analyses_temp[numanalyses++]=LoveAnalysisEnum;
			break;

		case EsaSolutionEnum:
			analyses_temp[numanalyses++]=EsaAnalysisEnum;
			break;

		case SamplingSolutionEnum:
			analyses_temp[numanalyses++]=SamplingAnalysisEnum;
			break;

		case SmbSolutionEnum:
			analyses_temp[numanalyses++]=SmbAnalysisEnum;
			break;

		case DamageEvolutionSolutionEnum:
			analyses_temp[numanalyses++]=DamageEvolutionAnalysisEnum;
			break;

		case TransientSolutionEnum:{
			/*We have multiple analyses here, process one by one*/
			bool isSIA,isFS,isthermal,isenthalpy,ismasstransport,ismmemasstransport,isoceantransport,isgroundingline,isstressbalance,ismovingfront,ishydrology,isdamage,issmb,isslc,isesa,isdebris,issampling,isfreesurface;
			iomodel->FindConstant(&isthermal,"md.transient.isthermal");
			iomodel->FindConstant(&ismovingfront,"md.transient.ismovingfront");
			iomodel->FindConstant(&ismasstransport,"md.transient.ismasstransport");
			iomodel->FindConstant(&ismmemasstransport,"md.transient.ismmemasstransport");
			iomodel->FindConstant(&isoceantransport,"md.transient.isoceantransport");
			iomodel->FindConstant(&isstressbalance,"md.transient.isstressbalance");
			iomodel->FindConstant(&isgroundingline,"md.transient.isgroundingline");
			iomodel->FindConstant(&isdamage,"md.transient.isdamageevolution");
			iomodel->FindConstant(&ishydrology,"md.transient.ishydrology");
			iomodel->FindConstant(&issmb,"md.transient.issmb");
			iomodel->FindConstant(&isfreesurface,"md.masstransport.isfreesurface");
			iomodel->FindConstant(&isslc,"md.transient.isslc");
			iomodel->FindConstant(&isesa,"md.transient.isesa");
			iomodel->FindConstant(&isdebris,"md.transient.isdebris");
			iomodel->FindConstant(&issampling,"md.transient.issampling");
      int* analyses_iter     = NULL;
      int  num_analyses_iter = 0;
			if(isstressbalance){
				SolutionAnalysesList(&analyses_iter,&num_analyses_iter,iomodel,StressbalanceSolutionEnum);
            xMemCpy<int>(&analyses_temp[numanalyses],analyses_iter,num_analyses_iter);
				numanalyses+=num_analyses_iter; xDelete<int>(analyses_iter);
         }
			if(isthermal && iomodel->domaintype==Domain3DEnum){
				SolutionAnalysesList(&analyses_iter,&num_analyses_iter,iomodel,ThermalSolutionEnum);
            xMemCpy<int>(&analyses_temp[numanalyses],analyses_iter,num_analyses_iter);
				numanalyses+=num_analyses_iter; xDelete<int>(analyses_iter);
			}
			if(ismasstransport || isgroundingline){
				analyses_temp[numanalyses++]=MasstransportAnalysisEnum;
				if(isfreesurface){
					analyses_temp[numanalyses++]=FreeSurfaceBaseAnalysisEnum;
					analyses_temp[numanalyses++]=FreeSurfaceTopAnalysisEnum;
				}
				int  basalforcing_model;
				iomodel->FindConstant(&basalforcing_model,"md.basalforcings.model");
				if(basalforcing_model==BasalforcingsPicoEnum){
					bool isplume;
					iomodel->FindConstant(&isplume,"md.basalforcings.isplume");
					if(isplume){
						analyses_temp[numanalyses++]=GLheightadvectionAnalysisEnum;
					}
				}
			}
			if(issmb) analyses_temp[numanalyses++]=SmbAnalysisEnum;
			if(ismovingfront){
				analyses_temp[numanalyses++]=ExtrapolationAnalysisEnum;
				analyses_temp[numanalyses++]=LevelsetAnalysisEnum;
			}
			if(ishydrology){
				SolutionAnalysesList(&analyses_iter,&num_analyses_iter,iomodel,HydrologySolutionEnum);
            xMemCpy<int>(&analyses_temp[numanalyses],analyses_iter,num_analyses_iter);
				numanalyses+=num_analyses_iter; xDelete<int>(analyses_iter);
			}
			if(isdamage){
				analyses_temp[numanalyses++]=DamageEvolutionAnalysisEnum;
			}
			if(isoceantransport){
				analyses_temp[numanalyses++]=OceantransportAnalysisEnum;
			}
			if(ismmemasstransport){
				analyses_temp[numanalyses++]=MmemasstransportAnalysisEnum;
			}
			if(isslc){
				analyses_temp[numanalyses++]=SealevelchangeAnalysisEnum;
			}
			if(isesa){
				analyses_temp[numanalyses++]=EsaAnalysisEnum;
			}
			if(isdebris){
				analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
				analyses_temp[numanalyses++]=SmbAnalysisEnum;
				analyses_temp[numanalyses++]=ExtrudeFromTopAnalysisEnum;
				analyses_temp[numanalyses++]=DebrisAnalysisEnum;
			}
			if(issampling){
				analyses_temp[numanalyses++]=SamplingAnalysisEnum;
			}

			if(iomodel->domaintype==Domain2DverticalEnum || iomodel->domaintype==Domain3DEnum){
				analyses_temp[numanalyses++]=ExtrudeFromBaseAnalysisEnum;
				analyses_temp[numanalyses++]=ExtrudeFromTopAnalysisEnum;
			}
			analyses_temp[numanalyses++]=L2ProjectionBaseAnalysisEnum;
			}
			break;

		default:
			_error_("solution type: " << EnumToStringx(solutiontype) << " not supported yet!");
			break;
	}

	/*Copy analyses from temp to output*/
	_assert_(numanalyses<MAXANALYSES);
	analyses=xNew<int>(numanalyses);
	for(int i=0;i<numanalyses;i++) analyses[i]=analyses_temp[i];

	/*Assign output pointers:*/
	if(pnumanalyses) *pnumanalyses=numanalyses;
	if(panalyses)    *panalyses=analyses;
	else              xDelete<int>(analyses);
}/*}}}*/
void FemModel::Solve(void){/*{{{*/

	/*profiling: */
	bool profiling = false;
	IssmDouble solution_time;
	IssmDouble solution_flops;
	IssmDouble solution_memory;

	/*solution: */
	int solution_type;
	void (*solutioncore)(FemModel*)=NULL; //core solution function pointer

	_printf0_("call computational core:\n");

	/*Retrieve solution_type from parameters: */
	parameters->FindParam(&solution_type,SolutionTypeEnum);

	/*Figure out which solution core we are going to run with the current solution type: */
	WrapperCorePointerFromSolutionEnum(&solutioncore,this->parameters,solution_type);

	/*run solution core: */
	profiler->Start(CORE);
	solutioncore(this);
	profiler->Stop(CORE);

	/*run AD core if needed: */
	profiler->Start(ADCORE);
	ad_core(this);
	profiler->Stop(ADCORE);

	/*some profiling results for the core: */
	parameters->FindParam(&profiling,DebugProfilingEnum);
	if(profiling){

		solution_time=profiler->TotalTime(CORE);
		solution_flops=profiler->TotalFlops(CORE);
		solution_memory=profiler->Memory(CORE);

		_printf0_("\n======================================================\n");
		_printf0_("Core solution profiling\n");
		_printf0_("   elapsed time    : " << solution_time   << " Seconds\n");
		_printf0_("   number of flops : " << solution_flops  << " Flops\n");
		_printf0_("   memory used     : " << solution_memory << " Bytes\n");

		/*Individual cores*/
		_printf0_("\nIndividual core profiling\n");
		if(profiler->Used(THERMALCORE)) _printf0_("   "<<setw(40)<<left<<"Thermal core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(THERMALCORE) << " sec\n");
		if(profiler->Used(HYDROLOGYCORE)) _printf0_("   "<<setw(40)<<left<<"Hydrology core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(HYDROLOGYCORE) << " sec\n");
		if(profiler->Used(STRESSBALANCECORE)) _printf0_("   "<<setw(40)<<left<<"Stress balance core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(STRESSBALANCECORE) << " sec\n");
		if(profiler->Used(DAMAGECORE)) _printf0_("   "<<setw(40)<<left<<"Damage core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(DAMAGECORE) << " sec\n");
		if(profiler->Used(MOVINGFRONTCORE)) _printf0_("   "<<setw(40)<<left<<"Moving front core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(MOVINGFRONTCORE) << " sec\n");
		if(profiler->Used(MASSTRANSPORTCORE)) _printf0_("   "<<setw(40)<<left<<"Mass transport core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(MASSTRANSPORTCORE) << " sec\n");
		if(profiler->Used(SMBCORE)) _printf0_("   "<<setw(40)<<left<<"SMB core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SMBCORE) << " sec\n");
		if(profiler->Used(GROUNDINGLINECORE)) _printf0_("   "<<setw(40)<<left<<"Groundingline migration core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(GROUNDINGLINECORE) << " sec\n");
		if(profiler->Used(ESACORE)) _printf0_("   "<<setw(40)<<left<<"ESA core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(ESACORE) << " sec\n");
		if(profiler->Used(DEBRISCORE)) _printf0_("   "<<setw(40)<<left<<"DEBRIS core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SLRCORE) << " sec\n");
		if(profiler->Used(SLRCORE)) _printf0_("   "<<setw(40)<<left<<"SLR core elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SLRCORE) << " sec\n");
		if(profiler->Used(MPISERIAL)) _printf0_("   "<<setw(40)<<left<<"MPISERIAL elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(MPISERIAL) << " sec\n");

		if(profiler->Used(SEDLOOP)) _printf0_("   "<<setw(40)<<left<<"SedimentLoop elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SEDLOOP) << " sec\n");
		if(profiler->Used(SEDMatrix)) _printf0_("   "<<setw(40)<<left<<"SedimentMatrix elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SEDMatrix) << " sec\n");
		if(profiler->Used(SEDUpdate)) _printf0_("   "<<setw(40)<<left<<"SedimentUpdate elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(SEDUpdate) << " sec\n");
		if(profiler->Used(EPLLOOP)) _printf0_("   "<<setw(40)<<left<<"EplLoop elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(EPLLOOP) << " sec\n");
		if(profiler->Used(EPLMasking)) _printf0_("   "<<setw(40)<<left<<"EPL masking elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(EPLMasking) << " sec\n");
		if(profiler->Used(EPLMatrices)) _printf0_("   "<<setw(40)<<left<<"EPLMatrices elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(EPLMatrices) << " sec\n");
		if(profiler->Used(EPLUpdate)) _printf0_("   "<<setw(40)<<left<<"EPLUpdate elapsed time:"<<setw(7)<<setprecision(6)<<profiler->TotalTime(EPLUpdate) << " sec\n");

		/*Add to results: */
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingSolutionTimeEnum,  solution_time));
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingCurrentMemEnum,  solution_memory));
		results->AddObject(new GenericExternalResult<IssmDouble>(results->Size()+1, ProfilingCurrentFlopsEnum, solution_flops));

		#ifdef _HAVE_AD_
		solution_time   = profiler->TotalTime(ADCORE);
		solution_flops  = profiler->TotalFlops(ADCORE);
		solution_memory = profiler->Memory(ADCORE);

		_printf0_("AD profiling\n");
		_printf0_("   elapsed time    : " << solution_time   << " Seconds\n");
		_printf0_("   number of flops : " << solution_flops  << " Flops\n");
		_printf0_("   memory used     : " << solution_memory << " Bytes\n");
		#endif
		_printf0_("======================================================\n");
	}
}
/*}}}*/

/*Modules:*/
void FemModel::BalancethicknessMisfitx(IssmDouble* presponse){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  weight,vx,vy,H,dvx[2],dvy[2],dH[2];
	IssmDouble  temp,Jdet,dhdt,groundedice_melting,surface_mass_balance;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);
		DatasetInput* weights_input                   = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum);   _assert_(weights_input);
		Input* thickness_input                 = element->GetInput(ThicknessEnum); _assert_(thickness_input);
		Input* vx_input                        = element->GetInput(VxEnum);                                  _assert_(vx_input);
		Input* vy_input                        = element->GetInput(VyEnum);                                  _assert_(vy_input);
		Input* surface_mass_balance_input      = element->GetInput(SmbMassBalanceEnum);          _assert_(surface_mass_balance_input);
		Input* groundedice_melting_input       = element->GetInput(BasalforcingsGroundediceMeltingRateEnum); _assert_(groundedice_melting_input);
		Input* dhdt_input                      = element->GetInput(BalancethicknessThickeningRateEnum);      _assert_(dhdt_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,BalancethicknessMisfitEnum);
			thickness_input->GetInputValue(&H, gauss);
			thickness_input->GetInputDerivativeValue(&dH[0],xyz_list,gauss);
			surface_mass_balance_input->GetInputValue(&surface_mass_balance,gauss);
			groundedice_melting_input->GetInputValue(&groundedice_melting,gauss);
			dhdt_input->GetInputValue(&dhdt,gauss);
			vx_input->GetInputValue(&vx,gauss);
			vx_input->GetInputDerivativeValue(&dvx[0],xyz_list,gauss);
			vy_input->GetInputValue(&vy,gauss);
			vy_input->GetInputDerivativeValue(&dvy[0],xyz_list,gauss);

			/*Balance thickness soft constraint J = 1/2 (div(Hv)-a)^2*/
			temp  = vx*dH[0]+vy*dH[1]+H*(dvx[0]+dvy[1]) - (surface_mass_balance-groundedice_melting-dhdt);
			J    +=weight*1/2*temp*temp*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*presponse=J;

}/*}}}*/
void FemModel::CalvingRateVonmisesx(){/*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->CalvingRateVonmises();
	}
}
/*}}}*/
void FemModel::CalvingRateLevermannx(){/*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->CalvingRateLevermann();
	}
}
/*}}}*/
void FemModel::CalvingFluxLevelsetx(){/*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->CalvingFluxLevelset();
	}
}
/*}}}*/
void FemModel::CalvingMeltingFluxLevelsetx(){/*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->CalvingMeltingFluxLevelset();
	}
}
/*}}}*/
void FemModel::CostFunctionx(IssmDouble* pJ,IssmDouble** pJlist,int* pn){/*{{{*/

	/*Intermediary*/
	int      num_responses;
	int     *responses      = NULL;
	Results *cost_functions = NULL;

	/*Recover parameters*/
	parameters->FindParam(&num_responses,InversionNumCostFunctionsEnum);
	parameters->FindParam(&responses,NULL,InversionCostFunctionsEnum);

	/*Get the value of all cost functions*/
	this->RequestedOutputsx(&cost_functions,responses,num_responses);

	/*Get and add all contributions one by one*/
	IssmDouble  J=0;
	IssmDouble* Jlist = xNew<IssmDouble>(num_responses);
	for(int i=0;i<num_responses;i++){
		ExternalResult* result=(ExternalResult*)cost_functions->GetObjectByOffset(i);
		Jlist[i] = reCast<IssmDouble>(result->GetValue());
		J       += Jlist[i];
	}
	_assert_(cost_functions->Size()==num_responses);

	/*Assign output pointers: */
	delete cost_functions;
	xDelete<int>(responses);
	if(pJ)     *pJ     = J;
	if(pJlist) *pJlist = Jlist;
	else        xDelete<IssmDouble>(Jlist);
	if(pn)     *pn     = num_responses;
}
/*}}}*/
void FemModel::DeviatoricStressx(){/*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->ComputeDeviatoricStressTensor();
	}
}
/*}}}*/
void FemModel::DistanceToFieldValue(int fieldenum,IssmDouble fieldvalue,int distanceenum){/*{{{*/

	/*recover my_rank:*/
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Get domain type (2d or 3d)*/
	int domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);

	/*1: go throug all elements of this partition and figure out how many
	 * segments we have (corresopnding to field = value)*/
	DataSet* segments=new DataSet();
	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		if(!element->IsOnBase()) continue;
		Element* basalelement = element->SpawnBasalElement();
		basalelement->WriteFieldIsovalueSegment(segments,fieldenum,fieldvalue);
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	}

	/*2: now get the segments from all partitions*/
	int  segcount=segments->Size();
	int* allsegcount=xNew<int>(num_procs);
	ISSM_MPI_Gather(&segcount,1,ISSM_MPI_INT,allsegcount,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(allsegcount,num_procs,ISSM_MPI_INT,0,IssmComm::GetComm());

	/* Every cpu should start its own dof count at the end of the dofcount from cpu-1*/
	int numseg_offset=0;
	int numseg=0;
	for(int i=0;i<my_rank;  i++) numseg_offset+=allsegcount[i];
	for(int i=0;i<num_procs;i++) numseg+=allsegcount[i];
	IssmDouble* segmentlist    = xNewZeroInit<IssmDouble>(4*numseg);
	IssmDouble* allsegmentlist = xNewZeroInit<IssmDouble>(4*numseg);
	int i=0;
	for(Object* & object : segments->objects){
		Contour<IssmDouble>* segment=(Contour<IssmDouble>*)object;
		_assert_(segment->nods == 2);
		segmentlist[(numseg_offset+i)*4 + 0] = segment->x[0];
		segmentlist[(numseg_offset+i)*4 + 1] = segment->y[0];
		segmentlist[(numseg_offset+i)*4 + 2] = segment->x[1];
		segmentlist[(numseg_offset+i)*4 + 3] = segment->y[1];
		i++;
	}

	ISSM_MPI_Allreduce((void*)segmentlist,(void*)allsegmentlist,4*numseg,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	delete segments;
	xDelete<IssmDouble>(segmentlist);
	xDelete<int>(allsegcount);

	/*3: Add distance input to all elements*/
	IssmDouble* distances = xNew<IssmDouble>(vertices->Size());
	IssmDouble  d,xn,yn,dmin;
	int         last = -1;
	for(Object* & object : this->vertices->objects){
      Vertex* vertex= xDynamicCast<Vertex*>(object);
		IssmDouble x = vertex->x;
		IssmDouble y = vertex->y;

		/*Most of the time the last checked segment is the closest so start with that one*/
		if(last>0){
			dmin = pow(allsegmentlist[4*last+0] - x,2) + pow(y-allsegmentlist[4*last+1],2);
		}
		else{
			dmin = 1.e+50;
		}

		for(int i=0;i<numseg;i++){

			/*Skip if tip is more than 10xdmin away*/
			if(pow(allsegmentlist[4*i+0] - x,2)>10*dmin) continue;
			if(pow(allsegmentlist[4*i+1] - y,2)>10*dmin) continue;

			IssmDouble l2 = (allsegmentlist[4*i+2]-allsegmentlist[4*i+0])*(allsegmentlist[4*i+2]-allsegmentlist[4*i+0]) + (allsegmentlist[4*i+3]-allsegmentlist[4*i+1])*(allsegmentlist[4*i+3]-allsegmentlist[4*i+1]);

			/*Segment has a length of 0*/
			if(l2==0.){
				d = (x-allsegmentlist[4*i+0])*(x-allsegmentlist[4*i+0])+(y-allsegmentlist[4*i+1])*(y-allsegmentlist[4*i+1]);
				if(d<dmin){
					dmin = d;
					last = i;
				}
				continue;
			}

			/*Consider the line extending the segment, parameterized as v + t (w - v).
			 *We find projection of point p onto the line.
			 *It falls where t = [(p-v) . (w-v)] / |w-v|^2*/
			IssmDouble t = ((x-allsegmentlist[4*i+0])*(allsegmentlist[4*i+2]-allsegmentlist[4*i+0]) + (y-allsegmentlist[4*i+1])*(allsegmentlist[4*i+3]-allsegmentlist[4*i+1]))/l2;
			if(t < 0.0){
				// Beyond the 'v' end of the segment
				d = (x-allsegmentlist[4*i+0])*(x-allsegmentlist[4*i+0])+(y-allsegmentlist[4*i+1])*(y-allsegmentlist[4*i+1]);
			}
			else if (t > 1.0){
				// Beyond the 'w' end of the segment
				d = (x-allsegmentlist[4*i+2])*(x-allsegmentlist[4*i+2])+(y-allsegmentlist[4*i+3])*(y-allsegmentlist[4*i+3]);
			}
			else{
				// Projection falls on the segment
				xn = allsegmentlist[4*i+0] + t * (allsegmentlist[4*i+2] - allsegmentlist[4*i+0]);
				yn = allsegmentlist[4*i+1] + t * (allsegmentlist[4*i+3] - allsegmentlist[4*i+1]);
				d = (x-xn)*(x-xn)+(y-yn)*(y-yn);
			}

			if(d<dmin){
				dmin = d;
				last = i;
			}
		}

		/*Update signed distance*/
		_assert_(vertex->lid<vertices->Size());
		distances[vertex->lid] = sqrt(dmin);
	}

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		element->CreateDistanceInputFromSegmentlist(distances,distanceenum);
	}
	//InputUpdateFromVectorx(this,distances,distanceenum,VertexLIdEnum);

	/*Clean up and return*/
	xDelete<IssmDouble>(distances);
	xDelete<IssmDouble>(allsegmentlist);
}/*}}}*/
void FemModel::Divergencex(IssmDouble* pdiv){/*{{{*/

	IssmDouble local_divergence=0;
	IssmDouble total_divergence;

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		local_divergence+=element->Divergence();
	}
	ISSM_MPI_Reduce(&local_divergence,&total_divergence,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_divergence,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pdiv=total_divergence;

}/*}}}*/
void FemModel::ElementOperationx(void (Element::*function)(void)){ /*{{{*/

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		(element->*function)();
	}

}
/*}}}*/
void FemModel::ElementResponsex(IssmDouble* presponse,int response_enum){/*{{{*/

	int found=0;
	int sumfound=0;
	int cpu_found=-1;
	int index;
	IssmDouble response;
	Element* element=NULL;

	/*retrieve element we are interested in: */
	this->parameters->FindParam(&index,IndexEnum);
	int my_rank=IssmComm::GetRank();

	/*now, go through our elements, and retrieve the one with this id: index: */
	for(Object* & object : this->elements->objects){
      element= xDynamicCast<Element*>(object);
		if (element->Id()==index){
			found=1;
			cpu_found=my_rank;
			break;
		}
	}

	/*Broadcast whether we found the element: */
	ISSM_MPI_Allreduce ( &found,&sumfound,1,ISSM_MPI_INT,ISSM_MPI_SUM,IssmComm::GetComm());
	if(!sumfound)_error_("could not find material with id" << index << " to compute ElementResponse");

	/*Ok, we found the element, compute responseocity: */
	if(my_rank==cpu_found){
		element->ElementResponse(&response,response_enum);
	}

	/*Broadcast and plug into response: */
	ISSM_MPI_Allreduce ( &cpu_found,&cpu_found,1,ISSM_MPI_INT,ISSM_MPI_MAX,IssmComm::GetComm());
	ISSM_MPI_Bcast(&response,1,ISSM_MPI_DOUBLE,cpu_found,IssmComm::GetComm());

	/*Assign output pointers: */
	*presponse=response;

}/*}}}*/
void FemModel::FloatingAreax(IssmDouble* pV, bool scaled){/*{{{*/

	IssmDouble local_floating_area= 0;
	IssmDouble total_floating_area;

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		local_floating_area+=element->FloatingArea(scaled);
	}
	ISSM_MPI_Reduce(&local_floating_area,&total_floating_area,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_floating_area,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_floating_area;

}/*}}}*/
void FemModel::GetInputLocalMinMaxOnNodesx(IssmDouble** pmin,IssmDouble** pmax,IssmDouble* ug){/*{{{*/

	/*Initialize output vectors*/
	int numnodes = this->nodes->NumberOfNodes();
	IssmDouble* uLmin_local = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmax_local = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmin = xNew<IssmDouble>(numnodes);
	IssmDouble* uLmax = xNew<IssmDouble>(numnodes);
	for(int i=0;i<numnodes;i++){
		uLmin_local[i] = +1.e+50;
		uLmax_local[i] = -1.e+50;
	}

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		element->GetInputLocalMinMaxOnNodes(uLmin_local,uLmax_local,ug);
	}

	/*Synchronize all CPUs*/
	ISSM_MPI_Allreduce((void*)uLmin_local,(void*)uLmin,numnodes,ISSM_MPI_DOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());
	ISSM_MPI_Allreduce((void*)uLmax_local,(void*)uLmax,numnodes,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
	xDelete<IssmDouble>(uLmin_local);
	xDelete<IssmDouble>(uLmax_local);

	/*Assign output pointers: */
	*pmin=uLmin;
	*pmax=uLmax;

}/*}}}*/
void FemModel::GetLocalVectorWithClonesGset(IssmDouble** plocal_ug,Vector<IssmDouble> *ug){/*{{{*/

	this->nodes->GetLocalVectorWithClonesGset(plocal_ug,ug);

}/*}}}*/
void FemModel::GetLocalVectorWithClonesVertices(IssmDouble** plocal_vector,Vector<IssmDouble> *vector){/*{{{*/

	/*retrieve vertex info*/
	int localsize         = this->vertices->NumberOfVerticesLocalAll();
	int localsize_masters = this->vertices->NumberOfVerticesLocal();

	/*Get local vector of vector*/
	int        *indices_vector_masters = NULL;
	IssmDouble *local_vector_masters   = NULL;
	vector->GetLocalVector(&local_vector_masters,&indices_vector_masters);
	_assert_(localsize_masters==indices_vector_masters[localsize_masters-1] - indices_vector_masters[0]+1);
	xDelete<int>(indices_vector_masters);

	/*Now, extend vectors to account for clones (make vectors longer, for clones at the end)*/
	IssmDouble *local_vector  = xNew<IssmDouble>(localsize);
	xMemCpy<IssmDouble>(local_vector,local_vector_masters,localsize_masters);
	xDelete<IssmDouble>(local_vector_masters);

	/*Now send and receive vector for vertices on partition edge*/
	SyncLocalVectorWithClonesVertices(local_vector);

	/*Assign output pointer*/
	*plocal_vector = local_vector;
}/*}}}*/
void FemModel::SyncLocalVectorWithClonesVertices(IssmDouble* local_vector){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Now send and receive vector for vertices on partition edge*/
	IssmDouble **send_buffers = xNewZeroInit<IssmDouble*>(num_procs);
	IssmDouble  *recv_buffer  = xNewZeroInit<IssmDouble>(this->vertices->Size());
	ISSM_MPI_Request  *send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for (int rank = 0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;

	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_send[rank]){
			int  numids = this->vertices->common_send[rank];
			send_buffers[rank] = xNew<IssmDouble>(numids,"t"); //only one alloc, "t" is required by adolc
			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_send_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(!vertex->clone);
            send_buffers[rank][i] = local_vector[vertex->lid];
			}
         ISSM_MPI_Isend(send_buffers[rank],numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_recv[rank]){
			int  numids = this->vertices->common_recv[rank];
         ISSM_MPI_Recv(recv_buffer,numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_recv_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(vertex->clone);
            local_vector[vertex->lid] = recv_buffer[i];
			}
		}
	}
   xDelete<IssmDouble>(recv_buffer);
   for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<IssmDouble>(send_buffers[rank]);
   }
   xDelete<IssmDouble*>(send_buffers);
   xDelete<ISSM_MPI_Request>(send_requests);
}/*}}}*/
void FemModel::SyncLocalVectorWithClonesVerticesAdd(IssmDouble* local_vector){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*Now send and receive vector for vertices on partition edge*/
	IssmDouble **send_buffers = xNewZeroInit<IssmDouble*>(num_procs);
	IssmDouble  *recv_buffer  = xNewZeroInit<IssmDouble>(this->vertices->Size());
	ISSM_MPI_Request  *send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for (int rank = 0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;

	/*1st: add slaves to master values (reverse of what we usually do)*/
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_recv[rank]){
			int  numids = this->vertices->common_recv[rank];
         send_buffers[rank] = xNew<IssmDouble>(numids,"t"); //only one alloc, "t" is required by adolc
			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_recv_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(vertex->clone);
				send_buffers[rank][i] = local_vector[vertex->lid];
			}
			ISSM_MPI_Isend(send_buffers[rank],numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_send[rank]){
			int  numids = this->vertices->common_send[rank];
			ISSM_MPI_Recv(recv_buffer,numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_send_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(!vertex->clone);
				local_vector[vertex->lid] += recv_buffer[i];
			}
		}
	}

	/*Wait until MPI is done*/
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_recv[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
	}

	/*Now sync masters across partitions*/
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_send[rank]){
			int  numids = this->vertices->common_send[rank];
			xDelete<IssmDouble>(send_buffers[rank]);
			send_buffers[rank] = xNew<IssmDouble>(numids,"t"); //only one alloc, "t" is required by adolc
			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_send_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(!vertex->clone);
				send_buffers[rank][i] = local_vector[vertex->lid];
			}
			ISSM_MPI_Isend(send_buffers[rank],numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_recv[rank]){
			int  numids = this->vertices->common_recv[rank];
			ISSM_MPI_Recv(recv_buffer,numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&status);

			for(int i=0;i<numids;i++){
				int   master_lid = this->vertices->common_recv_ids[rank][i];
				Vertex* vertex=xDynamicCast<Vertex*>(this->vertices->GetObjectByOffset(master_lid));
				_assert_(vertex->clone);
				local_vector[vertex->lid] = recv_buffer[i];
			}
		}
	}
	xDelete<IssmDouble>(recv_buffer);
	for(int rank=0;rank<num_procs;rank++){
		if(this->vertices->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<IssmDouble>(send_buffers[rank]);
	}
	xDelete<IssmDouble*>(send_buffers);
	xDelete<ISSM_MPI_Request>(send_requests);
}/*}}}*/
void FemModel::GetLocalVectorWithClonesNodes(IssmDouble** plocal_vector,Vector<IssmDouble> *vector){/*{{{*/

	/*recover my_rank:*/
	ISSM_MPI_Status status;
	int my_rank   = IssmComm::GetRank();
	int num_procs = IssmComm::GetSize();

	/*retrieve vertex info*/
	int localsize         = this->nodes->NumberOfNodesLocalAll();
	int localsize_masters = this->nodes->NumberOfNodesLocal();

	/*Get local vector of vector*/
	int        *indices_vector_masters = NULL;
	IssmDouble *local_vector_masters   = NULL;
	vector->GetLocalVector(&local_vector_masters,&indices_vector_masters);
	_assert_(localsize_masters==indices_vector_masters[localsize_masters-1] - indices_vector_masters[0]+1);
	xDelete<int>(indices_vector_masters);

	/*Now, extend vectors to account for clones (make vectors longer, for clones at the end)*/
	IssmDouble *local_vector  = xNew<IssmDouble>(localsize);
	xMemCpy<IssmDouble>(local_vector,local_vector_masters,localsize_masters);
	xDelete<IssmDouble>(local_vector_masters);

	/*Now send and receive vector for nodes on partition edge*/
	IssmDouble **send_buffers = xNewZeroInit<IssmDouble*>(num_procs);
	IssmDouble  *recv_buffer  = xNewZeroInit<IssmDouble>(this->nodes->Size(),"t"); //only one alloc, "t" is required by adolc
	ISSM_MPI_Request  *send_requests = xNew<ISSM_MPI_Request>(num_procs);
	for (int rank = 0;rank<num_procs;rank++) send_requests[rank] = ISSM_MPI_REQUEST_NULL;

	for(int rank=0;rank<num_procs;rank++){
		if(this->nodes->common_send[rank]){
			int  numids = this->nodes->common_send[rank];
			send_buffers[rank] = xNew<IssmDouble>(numids,"t"); //only one alloc, "t" is required by adolc
			for(int i=0;i<numids;i++){
				int   master_lid = this->nodes->common_send_ids[rank][i];
				Node* vertex=xDynamicCast<Node*>(this->nodes->GetObjectByOffset(master_lid));
				_assert_(!vertex->clone);
				send_buffers[rank][i] = local_vector[vertex->lid];
			}
			ISSM_MPI_Isend(send_buffers[rank],numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&send_requests[rank]);
		}
	}
	for(int rank=0;rank<num_procs;rank++){
		if(this->nodes->common_recv[rank]){
			int  numids = this->nodes->common_recv[rank];
			ISSM_MPI_Recv(recv_buffer,numids,ISSM_MPI_DOUBLE,rank,0,IssmComm::GetComm(),&status);
			for(int i=0;i<numids;i++){
				int   master_lid = this->nodes->common_recv_ids[rank][i];
				Node* vertex=xDynamicCast<Node*>(this->nodes->GetObjectByOffset(master_lid));
				_assert_(vertex->clone);
				local_vector[vertex->lid] = recv_buffer[i];
			}
		}
	}

	xDelete<IssmDouble>(recv_buffer);
	for(int rank=0;rank<num_procs;rank++){
		if(this->nodes->common_send[rank]) ISSM_MPI_Wait(&send_requests[rank],&status);
		xDelete<IssmDouble>(send_buffers[rank]);
	}
	xDelete<IssmDouble*>(send_buffers);
	xDelete<ISSM_MPI_Request>(send_requests);

	/*Assign output pointer*/
	*plocal_vector = local_vector;
}/*}}}*/
void FemModel::GroundedAreax(IssmDouble* pV, bool scaled){/*{{{*/

	IssmDouble local_grounded_area= 0;
	IssmDouble total_grounded_area;

	for(Object* & object : this->elements->objects){
      Element* element= xDynamicCast<Element*>(object);
		local_grounded_area+=element->GroundedArea(scaled);
	}
	ISSM_MPI_Reduce(&local_grounded_area,&total_grounded_area,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_grounded_area,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_grounded_area;

}/*}}}*/
void FemModel::IcefrontAreax(){/*{{{*/

	int numbasins,BasinId;
	this->parameters->FindParam(&numbasins,FrontalForcingsNumberofBasinsEnum);
	IssmDouble* basin_icefront_area = xNewZeroInit<IssmDouble>(numbasins);

	for(int basin=0;basin<numbasins;basin++){
		IssmDouble local_icefront_area = 0;
		IssmDouble total_icefront_area;

		/*Add contribution of each element to icefront area of the basin*/
		for(Object* & object : this->elements->objects){
			Element* element= xDynamicCast<Element*>(object);
			element->GetInputValue(&BasinId,FrontalForcingsBasinIdEnum);
			if(BasinId==basin) local_icefront_area+=element->GetIcefrontArea();
		}
		ISSM_MPI_Reduce(&local_icefront_area,&total_icefront_area,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
		ISSM_MPI_Bcast(&total_icefront_area,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		basin_icefront_area[basin]=total_icefront_area;
	}

	this->parameters->AddObject(new DoubleVecParam(FrontalForcingsBasinIcefrontAreaEnum,basin_icefront_area,numbasins));

	xDelete<IssmDouble>(basin_icefront_area);
}/*}}}*/
void FemModel::IcefrontMassFluxx(IssmDouble* pM, bool scaled){/*{{{*/

	IssmDouble local_mass_flux = 0;
	IssmDouble total_mass_flux;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_mass_flux+=element->IcefrontMassFlux(scaled);
	}
	ISSM_MPI_Reduce(&local_mass_flux,&total_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_mass_flux,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_mass_flux;

}/*}}}*/
void FemModel::IcefrontMassFluxLevelsetx(IssmDouble* pM, bool scaled){/*{{{*/

	IssmDouble local_mass_flux = 0;
	IssmDouble total_mass_flux;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_mass_flux+=element->IcefrontMassFluxLevelset(scaled);
	}
	ISSM_MPI_Reduce(&local_mass_flux,&total_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_mass_flux,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_mass_flux;

}/*}}}*/
void FemModel::InputMakeDiscontinuous(int enum_in){/*{{{*/

	int numvertices  = 6;
	IssmDouble* P1DGlist = xNew<IssmDouble>(numvertices);

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		element->GetInputListOnVertices(P1DGlist,enum_in);
		element->AddInput(DummyEnum,P1DGlist,P1DGEnum);
	}
	xDelete<IssmDouble>(P1DGlist);

	this->inputs->ChangeEnum(DummyEnum,enum_in);
	this->inputs->DeleteInput(DummyEnum);

}/*}}}*/
void FemModel::GroundinglineMassFluxx(IssmDouble* pM, bool scaled){/*{{{*/

	/*First we need to depth average the velocities*/
	int domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);
	this->parameters->SetParam(VxEnum,InputToDepthaverageInEnum);
	this->parameters->SetParam(VxAverageEnum,InputToDepthaverageOutEnum);
	depthaverage_core(this);
	if(domaintype!=Domain2DverticalEnum){
		this->parameters->SetParam(VyEnum,InputToDepthaverageInEnum);
		this->parameters->SetParam(VyAverageEnum,InputToDepthaverageOutEnum);
		depthaverage_core(this);
	}

	IssmDouble local_mass_flux = 0;
	IssmDouble total_mass_flux;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_mass_flux+=element->GroundinglineMassFlux(scaled);
	}
	ISSM_MPI_Reduce(&local_mass_flux,&total_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_mass_flux,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_mass_flux;

}/*}}}*/
void FemModel::IceMassx(IssmDouble* pM, bool scaled){/*{{{*/

	IssmDouble local_ice_mass = 0;
	IssmDouble total_ice_mass;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_ice_mass+=element->IceMass(scaled);
	}
	ISSM_MPI_Reduce(&local_ice_mass,&total_ice_mass,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_mass,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_ice_mass;

}/*}}}*/
void FemModel::IceVolumeAboveFloatationx(IssmDouble* pV, bool scaled){/*{{{*/

	IssmDouble local_ice_volume_af = 0;
	IssmDouble total_ice_volume_af;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_ice_volume_af+=element->IceVolumeAboveFloatation(scaled);
	}
	ISSM_MPI_Reduce(&local_ice_volume_af,&total_ice_volume_af,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_volume_af,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_ice_volume_af;

}/*}}}*/
void FemModel::IceVolumex(IssmDouble* pV, bool scaled){/*{{{*/

	IssmDouble local_ice_volume = 0;
	IssmDouble total_ice_volume;

	for(Object* & object : this->elements->objects){
		Element* element= xDynamicCast<Element*>(object);
		local_ice_volume+=element->IceVolume(scaled);
	}
	ISSM_MPI_Reduce(&local_ice_volume,&total_ice_volume,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_ice_volume,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pV=total_ice_volume;

}/*}}}*/
void FemModel::InputToP0(int inputenum,int outputenum){/*{{{*/

	IssmDouble average;

	/*Collapse input to P0, by doing the average. We need to have the elements 
	 * to do so, so loop onto elements: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(inputenum);
		input->GetInputAverage(&average);
		element->AddInput(outputenum,&average,P0Enum);
	}
}/*}}}*/
void FemModel::MassFluxx(IssmDouble* pmass_flux){/*{{{*/

	int          i,j;
	Element     *element       = NULL;
	int          element_id;
	bool         ispresent     = false;
	IssmDouble   mass_flux     = 0;
	IssmDouble   all_mass_flux = 0;
	int          counter;
	IssmDouble **array         = NULL;
	int          M;
	int         *mdims_array   = NULL;
	int         *ndims_array   = NULL;
	IssmDouble  *segments      = NULL;
	int          num_segments;

	/*First, figure out which segment to compute our mass flux on. Start with retrieving qmu_mass_flux_segments: */
	this->parameters->FindParam(&ispresent,MassFluxSegmentsPresentEnum);
	if(!ispresent)_error_("no mass flux segments available!");
	this->parameters->FindParam(&array,&M,&mdims_array,&ndims_array,MassFluxSegmentsEnum);

	/*Retrieve index of segments being used for MassFlux computation: */
	parameters->FindParam(&counter,IndexEnum);

	/*retrieve segments from array: */
	segments     = array[counter-1]; //matlab to "C" indexing
	num_segments = mdims_array[counter-1];

	/*Go through segments, and then elements, and figure out which elements belong to a segment.
	 * When we find one, use the element to compute the mass flux on the segment: */
	for(i=0;i<num_segments;i++){
		element_id=reCast<int,IssmDouble>(*(segments+5*i+4));
		for(Object* & object : this->elements->objects){
			element = xDynamicCast<Element*>(object);
			if (element->Id()==element_id){
				/*We found the element which owns this segment, use it to compute the mass flux: */
				mass_flux+=element->MassFlux(segments+5*i+0);
				break;
			}
		}
	}

	ISSM_MPI_Allreduce ( (void*)&mass_flux,(void*)&all_mass_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
	mass_flux=all_mass_flux;

	/*Free resources:*/
	for(i=0;i<M;i++){
		IssmDouble* matrix=array[i];
		xDelete<IssmDouble>(matrix);
	}
	xDelete<int>(mdims_array);
	xDelete<int>(ndims_array);
	xDelete<IssmDouble*>(array);

	/*Assign output pointers: */
	*pmass_flux=mass_flux;

}/*}}}*/
void FemModel::MaxAbsVxx(IssmDouble* pmaxabsvx){/*{{{*/

	int i;
	IssmDouble maxabsvx;
	IssmDouble node_maxabsvx;
	IssmDouble element_maxabsvx;

	/*Go through elements, and request velocity: */
	maxabsvx=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VxEnum);
		element_maxabsvx=input->GetInputMaxAbs();
		if(element_maxabsvx>maxabsvx) maxabsvx=element_maxabsvx;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvx,&node_maxabsvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxabsvx=node_maxabsvx;

	/*Assign output pointers:*/
	*pmaxabsvx=maxabsvx;

}/*}}}*/
void FemModel::MaxAbsVyx(IssmDouble* pmaxabsvy){/*{{{*/

	int i;
	IssmDouble maxabsvy;
	IssmDouble node_maxabsvy;
	IssmDouble element_maxabsvy;

	/*Go through elements, and request velocity: */
	maxabsvy=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VyEnum);
		element_maxabsvy=input->GetInputMaxAbs();
		if(element_maxabsvy>maxabsvy) maxabsvy=element_maxabsvy;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvy,&node_maxabsvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxabsvy=node_maxabsvy;

	/*Assign output pointers:*/
	*pmaxabsvy=maxabsvy;

}/*}}}*/
void FemModel::MaxAbsVzx(IssmDouble* pmaxabsvz){/*{{{*/

	int i;
	IssmDouble maxabsvz;
	IssmDouble node_maxabsvz;
	IssmDouble element_maxabsvz;

	/*Go through elements, and request velocity: */
	maxabsvz=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VzEnum);
		element_maxabsvz=input->GetInputMaxAbs();
		if(element_maxabsvz>maxabsvz) maxabsvz=element_maxabsvz;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxabsvz,&node_maxabsvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxabsvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxabsvz=node_maxabsvz;

	/*Assign output pointers:*/
	*pmaxabsvz=maxabsvz;

}/*}}}*/
void FemModel::MaxDivergencex(IssmDouble* pdiv){/*{{{*/

	IssmDouble local_divergence;
	IssmDouble node_max_divergence;
	IssmDouble max_divergence = -INFINITY;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_divergence=element->Divergence();
		if(fabs(local_divergence)>max_divergence) max_divergence=fabs(local_divergence);
	}
	ISSM_MPI_Reduce(&max_divergence,&node_max_divergence,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_max_divergence,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	max_divergence=node_max_divergence;

	/*Assign output pointers: */
	*pdiv=max_divergence;

}/*}}}*/
void FemModel::MaxVelx(IssmDouble* pmaxvel){/*{{{*/

	int i;
	IssmDouble maxvel;
	IssmDouble node_maxvel;
	IssmDouble element_maxvel;

	/*Go through elements, and request velocity: */
	maxvel=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input* vel_input = element->GetInput(VelEnum); _assert_(vel_input);
		element_maxvel = vel_input->GetInputMax();
		if(element_maxvel>maxvel) maxvel=element_maxvel;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvel,&node_maxvel,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvel,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxvel=node_maxvel;

	/*Assign output pointers:*/
	*pmaxvel=maxvel;

}/*}}}*/
void FemModel::MaxVxx(IssmDouble* pmaxvx){/*{{{*/

	int i;
	IssmDouble maxvx;
	IssmDouble node_maxvx;
	IssmDouble element_maxvx;

	/*Go through elements, and request velocity: */
	maxvx=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input* vx_input = element->GetInput(VxEnum); _assert_(vx_input);
		element_maxvx = vx_input->GetInputMax();
		if(element_maxvx>maxvx) maxvx=element_maxvx;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvx,&node_maxvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxvx=node_maxvx;

	/*Assign output pointers:*/
	*pmaxvx=maxvx;

}/*}}}*/
void FemModel::MaxVyx(IssmDouble* pmaxvy){/*{{{*/

	int i;
	IssmDouble maxvy;
	IssmDouble node_maxvy;
	IssmDouble element_maxvy;

	/*Go through elements, and request velocity: */
	maxvy=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input* vy_input = element->GetInput(VyEnum); _assert_(vy_input);
		element_maxvy = vy_input->GetInputMax();
		if(element_maxvy>maxvy) maxvy=element_maxvy;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvy,&node_maxvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxvy=node_maxvy;

	/*Assign output pointers:*/
	*pmaxvy=maxvy;

}/*}}}*/
void FemModel::MaxVzx(IssmDouble* pmaxvz){/*{{{*/

	int i;
	IssmDouble maxvz;
	IssmDouble node_maxvz;
	IssmDouble element_maxvz;

	/*Go through elements, and request velocity: */
	maxvz=-INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input* vz_input = element->GetInput(VzEnum); _assert_(vz_input);
		element_maxvz = vz_input->GetInputMax();
		if(element_maxvz>maxvz) maxvz=element_maxvz;
	}

	/*Figure out maximum across the cluster: */
	ISSM_MPI_Reduce(&maxvz,&node_maxvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_maxvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	maxvz=node_maxvz;

	/*Assign output pointers:*/
	*pmaxvz=maxvz;

}/*}}}*/
void FemModel::MinVelx(IssmDouble* pminvel){/*{{{*/

	int i;
	IssmDouble minvel;
	IssmDouble node_minvel;
	IssmDouble element_minvel;

	/*Go through elements, and request velocity: */
	minvel=INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VelEnum);
		element_minvel =input->GetInputMin();
		if(element_minvel<minvel) minvel=element_minvel;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvel,&node_minvel,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvel,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	minvel=node_minvel;

	/*Assign output pointers:*/
	*pminvel=minvel;

}/*}}}*/
void FemModel::MinVxx(IssmDouble* pminvx){/*{{{*/

	int i;
	IssmDouble minvx;
	IssmDouble node_minvx;
	IssmDouble element_minvx;

	/*Go through elements, and request velocity: */
	minvx=INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VxEnum);
		element_minvx =input->GetInputMin();
		if(element_minvx<minvx) minvx=element_minvx;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvx,&node_minvx,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvx,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	minvx=node_minvx;

	/*Assign output pointers:*/
	*pminvx=minvx;

}/*}}}*/
void FemModel::MinVyx(IssmDouble* pminvy){/*{{{*/

	int i;
	IssmDouble minvy;
	IssmDouble node_minvy;
	IssmDouble element_minvy;

	/*Go through elements, and request velocity: */
	minvy=INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VyEnum);
		element_minvy =input->GetInputMin();
		if(element_minvy<minvy) minvy=element_minvy;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvy,&node_minvy,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvy,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	minvy=node_minvy;

	/*Assign output pointers:*/
	*pminvy=minvy;

}/*}}}*/
void FemModel::MinVzx(IssmDouble* pminvz){/*{{{*/

	int i;
	IssmDouble minvz;
	IssmDouble node_minvz;
	IssmDouble element_minvz;

	/*Go through elements, and request velocity: */
	minvz=INFINITY;
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		Input*  input = element->GetInput(VzEnum);
		element_minvz =input->GetInputMin();
		if(element_minvz<minvz) minvz=element_minvz;
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce(&minvz,&node_minvz,1,ISSM_MPI_DOUBLE,ISSM_MPI_MAX,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_minvz,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	minvz=node_minvz;

	/*Assign output pointers:*/
	*pminvz=minvz;

}/*}}}*/
void FemModel::MmeToInputFromId(int id, int rootenum, int interpolationenum){ /*{{{*/

	MmeToInputFromIdx(this->inputs,this->elements,this->parameters,id,rootenum,interpolationenum);

}	//}}}
void FemModel::OmegaAbsGradientx( IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  omega,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		DatasetInput* weights_input = element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* omega_input   = element->GetInput(BalancethicknessOmegaEnum);              _assert_(omega_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,OmegaAbsGradientEnum);
			omega_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

			/*Tikhonov regularization: J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */
			//J+=weight*1/2*(dp[0]*dp[0]+dp[1]*dp[1])*Jdet*gauss->weight;
			J+=weight*1/2*pow(dp[0]*dp[0]+dp[1]*dp[1],2)*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::EtaDiffx( IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  omega,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  p,p0;

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		DatasetInput* weights_input =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* omega_input   =element->GetInput(BalancethicknessOmegaEnum);              _assert_(omega_input);
		Input* omega0_input  =element->GetInput(BalancethicknessOmega0Enum);             _assert_(omega0_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,EtaDiffEnum);
			omega_input->GetInputValue(&p,gauss);
			omega0_input->GetInputValue(&p0,gauss);

			/*Tikhonov regularization: J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */
			//J+=weight*1/2*(dp[0]*dp[0]+dp[1]*dp[1])*Jdet*gauss->weight;
			J+=weight*1/2*pow(p - p0,2)*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::OutputControlsx(Results **presults){/*{{{*/

	/*parameters: */
	int         num_controls,step;
	IssmDouble  time;
	int        *control_type   = NULL;
	int        *control_interp = NULL;
	int        *M = NULL;
	int        *N = NULL;

	/*recover results*/
	Results* results = *presults;
	if(!results) results = new Results();

	/*Get list of Controls*/
	this->parameters->FindParam(&num_controls,InversionNumControlParametersEnum);
	this->parameters->FindParam(&control_type,NULL,InversionControlParametersEnum);
	this->parameters->FindParam(&M,NULL,ControlInputSizeMEnum);
	this->parameters->FindParam(&N,NULL,ControlInputSizeNEnum);
	this->parameters->FindParam(&control_interp,NULL,ControlInputInterpolationEnum);
	this->parameters->FindParam(&step,StepEnum);
	this->parameters->FindParam(&time,TimeEnum);

	for(int i=0;i<num_controls;i++){

		int control_enum = control_type[i];
		int gradient_enum;

		switch(i){
			case 0: gradient_enum = Gradient1Enum; break;
			case 1: gradient_enum = Gradient2Enum; break;
			case 2: gradient_enum = Gradient3Enum; break;
			case 3: gradient_enum = Gradient4Enum; break;
			default: _error_("more than 4 controls not implemented yet");
		}

		/*Allocate vector*/
		Vector<IssmPDouble> *vector_control  = new Vector<IssmPDouble>(M[i]*N[i]);
		Vector<IssmPDouble> *vector_gradient = new Vector<IssmPDouble>(M[i]*N[i]);

		/*Fill in vector*/
		for(Object* & object : this->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->ControlToVectors(vector_control,vector_gradient,control_enum,control_interp[i]);
		}
		vector_control->Assemble();
		vector_gradient->Assemble();

		results->AddResult(new GenericExternalResult<Vector<IssmPDouble>*>(results->Size()+1,control_enum,vector_control ,step,time));
		results->AddResult(new GenericExternalResult<Vector<IssmPDouble>*>(results->Size()+1,gradient_enum,vector_gradient,step,time));
	}

	/*Clean up and return*/
	xDelete<int>(control_type);
	xDelete<int>(control_interp);
	xDelete<int>(M);
	xDelete<int>(N);
}
/*}}}*/
void FemModel::RequestedDependentsx(void){/*{{{*/

	/*AD mode on?: */
	bool isautodiff;
	parameters->FindParam(&isautodiff,AutodiffIsautodiffEnum);

	if(isautodiff){
		#ifdef _HAVE_AD_
		int      num_dependents;
		DataSet* dependent_objects=NULL;
		parameters->FindParam(&num_dependents,AutodiffNumDependentsEnum);
		parameters->FindParam(&dependent_objects,AutodiffDependentObjectsEnum);
		if(num_dependents){
			IssmPDouble* dependents=xNew<IssmPDouble>(num_dependents);

			/*Go through our dependent variables, and compute the response:*/
			int my_rank=IssmComm::GetRank();
			int i = 0;
			for(Object* & object : dependent_objects->objects){
				DependentObject* dep=(DependentObject*)object;
				dep->RecordResponsex(this);
				IssmDouble output_value = dep->GetValue();
				if (my_rank==0) {
					#if defined(_HAVE_CODIPACK_)
						codi_global.registerOutput(output_value);
					#else
						output_value>>=dependents[i];
					#endif
				}
				i++;
			}
			xDelete<IssmPDouble>(dependents);
		}
		delete dependent_objects;
		#else
		_error_("Should not be requesting dependents when an AD library is not available!");
		#endif
	}
}
/*}}}*/
void FemModel::RequestedOutputsx(Results **presults,char** requested_outputs, int numoutputs, bool save_results){/*{{{*/

	/*Intermediaries*/
	bool        isvec,results_on_nodes;
	int         step,output_enum,numonnodes,ierr;
	IssmDouble  time;
	IssmDouble  double_result;
	const char *output_string = NULL;
	char**      resultsonnodes = NULL;

	/*recover results*/
	Results* results = *presults;
	if(!results) results = new Results();

	/*Get time and step*/
	parameters->FindParam(&step,StepEnum);
	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&numonnodes,SettingsNumResultsOnNodesEnum);
	if(numonnodes) parameters->FindParam(&resultsonnodes,&numonnodes,SettingsResultsOnNodesEnum);

	/*Go through all requested output*/
	for(int i=0;i<numoutputs;i++){
		output_string = requested_outputs[i];
		output_enum   = StringToEnumx(output_string,false);
		isvec         = false;

		/*If string is not an enum, it is defined in output definitions*/
		if(output_enum<0){
			ierr = OutputDefinitionsResponsex(&double_result, this,output_string);
			if(save_results){
				if(!ierr){
					results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
				}
				continue;
			}
		}
		else{
			/*last chance for the output definition, if the enum is one of Outputdefinition[1-10]Enum:*/
			if(output_enum>=Outputdefinition1Enum && output_enum <=Outputdefinition100Enum){
				ierr = OutputDefinitionsResponsex(&double_result, this, output_enum);
				if(save_results){
					if(!ierr){
						results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
					}
					continue;
				}
			}
			else{
				switch(output_enum){

					/*Scalar output*/
					case DivergenceEnum:                     this->Divergencex(&double_result);                     break;
					case MaxDivergenceEnum:                  this->MaxDivergencex(&double_result);                  break;
					case IceMassEnum:                        this->IceMassx(&double_result,false);                  break;
					case IcefrontMassFluxEnum:               this->IcefrontMassFluxx(&double_result,false);         break;
					case IcefrontMassFluxLevelsetEnum:       this->IcefrontMassFluxLevelsetx(&double_result,false);         break;
					case IceMassScaledEnum:                  this->IceMassx(&double_result,true);                   break;
					case IceVolumeEnum:                      this->IceVolumex(&double_result,false);                break;
					case IceVolumeScaledEnum:                this->IceVolumex(&double_result,true);                 break;
					case IceVolumeAboveFloatationEnum:       this->IceVolumeAboveFloatationx(&double_result,false); break;
					case IceVolumeAboveFloatationScaledEnum: this->IceVolumeAboveFloatationx(&double_result,true);  break;
					case GroundedAreaEnum:                   this->GroundedAreax(&double_result,false);             break;
					case GroundedAreaScaledEnum:             this->GroundedAreax(&double_result,true);              break;
					case GroundinglineMassFluxEnum:          this->GroundinglineMassFluxx(&double_result,false);    break;
					case FloatingAreaEnum:                   this->FloatingAreax(&double_result,false);             break;
					case FloatingAreaScaledEnum:             this->FloatingAreax(&double_result,true);              break;
					case MinVelEnum:                         this->MinVelx(&double_result);                         break;
					case MaxVelEnum:                         this->MaxVelx(&double_result);                         break;
					case MinVxEnum:                          this->MinVxx(&double_result);                          break;
					case MaxVxEnum:                          this->MaxVxx(&double_result);                          break;
					case MaxAbsVxEnum:                       this->MaxAbsVxx(&double_result);                       break;
					case MinVyEnum:                          this->MinVyx(&double_result);                          break;
					case MaxVyEnum:                          this->MaxVyx(&double_result);                          break;
					case MaxAbsVyEnum:                       this->MaxAbsVyx(&double_result);                       break;
					case MinVzEnum:                          this->MinVzx(&double_result);                          break;
					case MaxVzEnum:                          this->MaxVzx(&double_result);                          break;
					case MaxAbsVzEnum:                       this->MaxAbsVzx(&double_result);                       break;
					case MassFluxEnum:                       this->MassFluxx(&double_result);                       break;
					case TotalCalvingFluxLevelsetEnum:		  this->TotalCalvingFluxLevelsetx(&double_result,false);          break;
					case TotalCalvingMeltingFluxLevelsetEnum:this->TotalCalvingMeltingFluxLevelsetx(&double_result,false);          break;
					case TotalFloatingBmbEnum:               this->TotalFloatingBmbx(&double_result,false);         break;
					case TotalFloatingBmbScaledEnum:         this->TotalFloatingBmbx(&double_result,true);          break;
					case TotalGroundedBmbEnum:               this->TotalGroundedBmbx(&double_result,false);         break;
					case TotalGroundedBmbScaledEnum:         this->TotalGroundedBmbx(&double_result,true);          break;
					case TotalSmbEnum:                       this->TotalSmbx(&double_result,false);                 break;
					case TotalSmbMeltEnum:                   this->TotalSmbMeltx(&double_result,false);             break;
					case TotalSmbRefreezeEnum:               this->TotalSmbRefreezex(&double_result,false);         break;
					case TotalSmbScaledEnum:                 this->TotalSmbx(&double_result,true);                  break;

					/*Scalar control output*/
					case SurfaceAbsVelMisfitEnum:       SurfaceAbsVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
					case SurfaceRelVelMisfitEnum:       SurfaceRelVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
					case SurfaceLogVelMisfitEnum:       SurfaceLogVelMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);        break;
					case SurfaceLogVxVyMisfitEnum:      SurfaceLogVxVyMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);       break;
					case SurfaceAverageVelMisfitEnum:   SurfaceAverageVelMisfitx(&double_result,this);                                                  break;
					case ThicknessAbsMisfitEnum:        ThicknessAbsMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);         break;
					case ThicknessAbsGradientEnum:      this->ThicknessAbsGradientx(&double_result);                                                    break;
					case ThicknessAlongGradientEnum:    ThicknessAlongGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);     break;
					case ThicknessAcrossGradientEnum:   ThicknessAcrossGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);    break;
					case ThicknessPositiveEnum:         this->ThicknessPositivex(&double_result);                                                       break;
					case RheologyBbarAbsGradientEnum:   RheologyBbarAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);    break;
					case RheologyBAbsGradientEnum:      RheologyBAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters);       break;
					case RheologyBInitialguessMisfitEnum:  RheologyBInitialguessMisfitx(&double_result,elements,nodes,vertices,loads,materials,parameters);       break;
					case DragCoefficientAbsGradientEnum:DragCoefficientAbsGradientx(&double_result,elements,nodes,vertices,loads,materials,parameters); break;
					case BalancethicknessMisfitEnum:    BalancethicknessMisfitx(&double_result);                                                        break;
					case SurfaceAbsMisfitEnum:          SurfaceAbsMisfitx(&double_result); break;
					case OmegaAbsGradientEnum:          OmegaAbsGradientx(&double_result); break;
					case EtaDiffEnum:                   EtaDiffx(&double_result); break;

					/*Vector special case (maybe should go to specific analysis?)*/
					case ChannelAreaEnum:
					case ChannelDischargeEnum:{

							/*Get Number of Channels*/
							int numchannels_local=0,numchannels;
							for(int j=0;j<this->loads->Size();j++){
								if(this->loads->GetEnum(j)==ChannelEnum) numchannels_local++;
							}
							ISSM_MPI_Reduce(&numchannels_local,&numchannels,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
							ISSM_MPI_Bcast(&numchannels,1,ISSM_MPI_INT,0,IssmComm::GetComm());

							IssmPDouble* values    = xNewZeroInit<IssmPDouble>(numchannels);
							IssmPDouble* allvalues = xNew<IssmPDouble>(numchannels);

							/*Fill-in vector*/
							for(Object* & object : this->loads->objects){
								if(object->ObjectEnum()==ChannelEnum){
									Channel* channel=(Channel*)object;
									if(output_enum==ChannelAreaEnum){
										channel->WriteChannelCrossSection(values);
									}
									else if(output_enum==ChannelDischargeEnum){
										channel->WriteChannelDischarge(values);
									}
									else{
										_error_("not supported");
									}
								}
							}

							/*Gather from all cpus*/
							ISSM_MPI_Allreduce((void*)values,(void*)allvalues,numchannels,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
							xDelete<IssmPDouble>(values);

							if(save_results)results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,output_enum,allvalues,numchannels,1,step,time));
							xDelete<IssmPDouble>(allvalues);

							isvec = true;
					}
					break;

				   /*Default is always Vector */
					default:

					/*Some preliminary calculation may be required (use similar syntax for other inputs)*/
						if(output_enum==NewDamageEnum){
							InputDuplicatex(this,DamageDEnum,DamageDOldEnum);
							InputDuplicatex(this,DamageDbarEnum,DamageDbarOldEnum);
							this->ElementOperationx(&Element::ComputeNewDamage);
						}
						else if(output_enum==FrictionAlpha2Enum){
							for(Object* & object : this->elements->objects){
								Element* element=xDynamicCast<Element*>(object);
								element->SetElementInput(FrictionAlpha2Enum,0.,P1Enum);
							}
							this->ElementOperationx(&Element::FrictionAlpha2CreateInput);
						}

						/*Vector layout*/
						if(!IsInputEnum(output_enum)) _error_("Cannot output \""<<EnumToStringx(output_enum)<<"\" because it is not an input");
						int interpolation,nodesperelement,size,nlines,ncols,array_size;
						int rank_interpolation=-1,rank_nodesperelement=-1,rank_arraysize=-1,max_rank_arraysize=0;
						bool isarray=false;

						/*Get interpolation (and compute input if necessary)*/
						for(Object* & object : this->elements->objects){
							Element* element = xDynamicCast<Element*>(object);
							element->ResultInterpolation(&rank_interpolation,&rank_nodesperelement,&rank_arraysize,output_enum);
							if(rank_arraysize>max_rank_arraysize)max_rank_arraysize=rank_arraysize;
						}
						rank_arraysize=max_rank_arraysize;

						/*Broadcast for cpus that do not have any elements*/
						ISSM_MPI_Reduce(&rank_interpolation,&interpolation,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Reduce(&rank_nodesperelement,&nodesperelement,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Reduce(&rank_arraysize,&array_size,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&interpolation,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&nodesperelement,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						ISSM_MPI_Bcast(&array_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

						results_on_nodes=false;
						/*Loop to see if this output was requested on nodes*/
						for(int j=0;j<numonnodes & results_on_nodes==false;j++){
							if(strcmp(resultsonnodes[j],output_string) == 0 || strcmp(resultsonnodes[j],"all") == 0) results_on_nodes=true;
						}

						if(results_on_nodes){

							/*Allocate matrices*/
							int         nbe       = this->elements->NumberOfElements();
							IssmDouble* values    = xNewZeroInit<IssmDouble>(nbe*nodesperelement);
							IssmDouble* allvalues = xNew<IssmDouble>(nbe*nodesperelement);

							/*Fill-in matrix*/
							for(Object* & object : this->elements->objects){
								Element* element = xDynamicCast<Element*>(object);
								element->ResultToPatch(values,nodesperelement,output_enum);
							}

							/*Gather from all cpus*/
							ISSM_MPI_Allreduce((void*)values,(void*)allvalues,nbe*nodesperelement,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
							xDelete<IssmDouble>(values);

							if(save_results)results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,output_enum,allvalues,nbe,nodesperelement,step,time));
							xDelete<IssmDouble>(allvalues);

						}
						else{

							/*Allocate vector depending on interpolation*/
							switch(interpolation){
								case P0Enum: isarray = false; size = this->elements->NumberOfElements(); break;
								case P1Enum: isarray = false; size = this->vertices->NumberOfVertices(); break;
								case P0ArrayEnum: isarray = true; nlines = this->elements->NumberOfElements(); ncols= array_size; break;
								default:     _error_("Interpolation "<<EnumToStringx(interpolation)<<" not supported yet");

							}
							if(!isarray){
								Vector<IssmDouble> *vector_result = new Vector<IssmDouble>(size);

								/*Fill in vector*/
								for(Object* & object : this->elements->objects){
									Element* element = xDynamicCast<Element*>(object);
									element->ResultToVector(vector_result,output_enum);
								}
								vector_result->Assemble();

								if(save_results){
									results->AddResult(new GenericExternalResult<Vector<IssmDouble>*>(results->Size()+1,output_enum,vector_result,step,time));
									/*We do not do a copy for Vectors, so don't delete*/
								}
								else{
									delete vector_result;
								}
							}
							else{
								IssmDouble* values    = xNewZeroInit<IssmDouble>(nlines*ncols);
								IssmDouble* allvalues = xNew<IssmDouble>(nlines*ncols);

								/*Fill-in matrix*/
								for(Object* & object : this->elements->objects){
									Element* element = xDynamicCast<Element*>(object);
									element->ResultToMatrix(values,ncols,output_enum);
								}
								/*Gather from all cpus*/
								ISSM_MPI_Allreduce((void*)values,(void*)allvalues,ncols*nlines,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
								xDelete<IssmDouble>(values);

								if(save_results)results->AddResult(new GenericExternalResult<IssmDouble*>(results->Size()+1,output_enum,allvalues,nlines,ncols,step,time));
								xDelete<IssmDouble>(allvalues);
							}
						}
						isvec = true;
						break;
				}
			}

		}

		/*Add result to Results*/
		if(!isvec && save_results){
			results->AddResult(new GenericExternalResult<IssmPDouble>(results->Size()+1,output_string,reCast<IssmPDouble>(double_result),step,time));
		}
	}

	/*Clean up*/
	for(int i=0;i<numonnodes;i++) xDelete<char>(resultsonnodes[i]);
	xDelete<char*>(resultsonnodes);

	/*Assign pointer and clean up*/
	*presults = results;
}
/*}}}*/
void FemModel::RequestedOutputsx(Results **presults,int* requested_outputs, int numoutputs,bool save_results){/*{{{*/

	/*Convert list of enums to list of string*/
	char** enumlist = xNew<char*>(numoutputs);
	for(int i=0;i<numoutputs;i++) EnumToStringx(&enumlist[i],requested_outputs[i]);

	/*Call main module*/
	this->RequestedOutputsx(presults,enumlist,numoutputs,save_results);

	/*clean up and return*/
	for(int i=0;i<numoutputs;i++) xDelete<char>(enumlist[i]);
	xDelete<char*>(enumlist);
}
/*}}}*/
void FemModel::ResetLevelset(void){/*{{{*/

	this->DistanceToFieldValue(MaskIceLevelsetEnum,0.,MaskIceLevelsetEnum);

}
/*}}}*/
void FemModel::Responsex(IssmDouble* responses,const char* response_descriptor){/*{{{*/

	int response_descriptor_enum=StringToEnumx(response_descriptor);
	this->Responsex(responses, response_descriptor_enum);

}
/*}}}*/
void FemModel::Responsex(IssmDouble* responses,int response_descriptor_enum){/*{{{*/

	switch (response_descriptor_enum){

		case DivergenceEnum:                     this->Divergencex(responses); break;
		case MaxDivergenceEnum:                  this->MaxDivergencex(responses); break;
		case IceMassEnum:                        this->IceMassx(responses, false); break;
		case IceMassScaledEnum:                  this->IceMassx(responses, true); break;
		case IceVolumeEnum:                      this->IceVolumex(responses, false); break;
		case IceVolumeScaledEnum:                this->IceVolumex(responses, true); break;
		case IceVolumeAboveFloatationEnum:       this->IceVolumeAboveFloatationx(responses, false); break;
		case IceVolumeAboveFloatationScaledEnum: this->IceVolumeAboveFloatationx(responses, true); break;
		case IcefrontMassFluxEnum:               this->IcefrontMassFluxx(responses, false); break;
		case IcefrontMassFluxLevelsetEnum:       this->IcefrontMassFluxLevelsetx(responses, false); break;
		case GroundedAreaEnum:                   this->GroundedAreax(responses, false); break;
		case GroundedAreaScaledEnum:             this->GroundedAreax(responses, true); break;
		case FloatingAreaEnum:                   this->FloatingAreax(responses, false); break;
		case FloatingAreaScaledEnum:             this->FloatingAreax(responses, true); break;
		case MinVelEnum:                         this->MinVelx(responses); break;
		case MaxVelEnum:                         this->MaxVelx(responses); break;
		case MinVxEnum:                          this->MinVxx(responses); break;
		case MaxVxEnum:                          this->MaxVxx(responses); break;
		case MaxAbsVxEnum:                       this->MaxAbsVxx(responses); break;
		case MinVyEnum:                          this->MinVyx(responses); break;
		case MaxVyEnum:                          this->MaxVyx(responses); break;
		case MaxAbsVyEnum:                       this->MaxAbsVyx(responses); break;
		case MinVzEnum:                          this->MinVzx(responses); break;
		case MaxVzEnum:                          this->MaxVzx(responses); break;
		case MaxAbsVzEnum:                       this->MaxAbsVzx(responses); break;
		case MassFluxEnum:                       this->MassFluxx(responses); break;
		case SurfaceAbsVelMisfitEnum:            SurfaceAbsVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceRelVelMisfitEnum:            SurfaceRelVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceLogVelMisfitEnum:            SurfaceLogVelMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceLogVxVyMisfitEnum:           SurfaceLogVxVyMisfitx(responses, elements,nodes, vertices, loads, materials,parameters); break;
		case SurfaceAverageVelMisfitEnum:        SurfaceAverageVelMisfitx(responses,this); break;
		case ThicknessAbsMisfitEnum:             ThicknessAbsMisfitx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case ThicknessAbsGradientEnum:           this->ThicknessAbsGradientx(responses); break;
		case ThicknessAlongGradientEnum:         ThicknessAlongGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case ThicknessAcrossGradientEnum:        ThicknessAcrossGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case RheologyBbarAbsGradientEnum:        RheologyBbarAbsGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case DragCoefficientAbsGradientEnum:     DragCoefficientAbsGradientx(responses, elements,nodes, vertices, loads, materials, parameters); break;
		case BalancethicknessMisfitEnum:         BalancethicknessMisfitx(responses); break;
		case TotalCalvingFluxLevelsetEnum:		  this->TotalCalvingFluxLevelsetx(responses, false); break;
		case TotalCalvingMeltingFluxLevelsetEnum:this->TotalCalvingMeltingFluxLevelsetx(responses, false); break;
		case TotalFloatingBmbEnum:			        this->TotalFloatingBmbx(responses, false); break;
		case TotalFloatingBmbScaledEnum:			  this->TotalFloatingBmbx(responses, true); break;
		case TotalGroundedBmbEnum:			        this->TotalGroundedBmbx(responses, false); break;
		case TotalGroundedBmbScaledEnum:			  this->TotalGroundedBmbx(responses, true); break;
		case TotalSmbEnum:					        this->TotalSmbx(responses, false); break;
		case TotalSmbMeltEnum:					     this->TotalSmbMeltx(responses, false); break;
		case TotalSmbRefreezeEnum:					  this->TotalSmbRefreezex(responses, false); break;
		case TotalSmbScaledEnum:					  this->TotalSmbx(responses, true); break;
		case MaterialsRheologyBbarEnum:          this->ElementResponsex(responses,MaterialsRheologyBbarEnum); break;
		case VelEnum:                            this->ElementResponsex(responses,VelEnum); break;
		case FrictionCoefficientEnum:            NodalValuex(responses, FrictionCoefficientEnum,elements,nodes, vertices, loads, materials, parameters); break;
		case GroundinglineMassFluxEnum:          this->GroundinglineMassFluxx(responses, false);    break;
		default:
			if(response_descriptor_enum>=Outputdefinition1Enum && response_descriptor_enum <=Outputdefinition2000Enum){
				int ierr = OutputDefinitionsResponsex(responses, this,response_descriptor_enum);
				if(ierr) _error_("could not evaluate response");
			}
			else _error_("response descriptor \"" << EnumToStringx(response_descriptor_enum) << "\" not supported yet!");
			break;
	}

}
/*}}}*/
void FemModel::RignotMeltParameterizationx(){/*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->RignotMeltParameterization();
	}
}
/*}}}*/
void FemModel::StrainRateparallelx(){/*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->StrainRateparallel();
	}
}
/*}}}*/
void FemModel::StrainRateperpendicularx(){/*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->StrainRateperpendicular();
	}
}
/*}}}*/
void FemModel::StrainRateeffectivex(){/*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->ComputeStrainRate();
	}
}
/*}}}*/
void FemModel::StressIntensityFactorx(){/*{{{*/

	/*Update input for basal element only*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->StressIntensityFactor();
	}
}
	/*}}}*/
void FemModel::SurfaceAbsMisfitx(IssmDouble* presponse){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  surface,surfaceobs,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		 /*Retrieve all inputs we will be needing: */
		 DatasetInput* weights_input   =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		 Input* surface_input   =element->GetInput(SurfaceEnum);                            _assert_(surface_input);
		 Input* surfaceobs_input=element->GetInput(InversionSurfaceObsEnum);                _assert_(surfaceobs_input);

		 /* Start  looping on the number of gaussian points: */
		 Gauss* gauss=element->NewGauss(2);
		 while(gauss->next()){

			 /* Get Jacobian determinant: */
			 element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			 /*Get all parameters at gaussian point*/
			 weights_input->GetInputValue(&weight,gauss,SurfaceAbsMisfitEnum);
			 surface_input->GetInputValue(&surface,gauss);
			 surfaceobs_input->GetInputValue(&surfaceobs,gauss);

			 /*Compute SurfaceAbsMisfitEnum*/
			 J+=0.5*(surface-surfaceobs)*(surface-surfaceobs)*weight*Jdet*gauss->weight;
		 }
		 delete gauss;
		 xDelete<IssmDouble>(xyz_list);
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*presponse=J;

}/*}}}*/
void FemModel::ThicknessAbsGradientx( IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  thickness,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  dp[3];

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		DatasetInput* weights_input   =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* thickness_input =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,ThicknessAbsGradientEnum);
			thickness_input->GetInputDerivativeValue(&dp[0],xyz_list,gauss);

			/*Tikhonov regularization: J = 1/2 ((dp/dx)^2 + (dp/dy)^2) */
			J+=weight*1/2*(dp[0]*dp[0]+dp[1]*dp[1])*Jdet*gauss->weight;
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::ThicknessAverage(){/*{{{*/

	int elementswidth                   = this->GetElementsWidth();//just 2D mesh, tria elements
   int numberofvertices                = this->vertices->NumberOfVertices();//total number of vertices

   IssmDouble weight                   = 0.;
   IssmDouble* totalweight             = NULL;
	IssmDouble* Hserial						= NULL;
   IssmDouble* H                       = xNew<IssmDouble>(elementswidth);
   Vector<IssmDouble>* vecH				= new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vectotalweight  = new Vector<IssmDouble>(numberofvertices);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*check if there is ice in this element*/
		if(!element->IsIceInElement()) continue;

		/*get H on the vertices*/
		element->GetInputListOnVertices(H,ThicknessEnum);

      /*weight to calculate the smoothed H*/
      weight=1.;//simple average

		/*add in the serial vector*/
      vecH->SetValue(element->vertices[0]->Sid(),weight*H[0],ADD_VAL);
      vecH->SetValue(element->vertices[1]->Sid(),weight*H[1],ADD_VAL);
      vecH->SetValue(element->vertices[2]->Sid(),weight*H[2],ADD_VAL);
      /*total weight*/
      vectotalweight->SetValue(element->vertices[0]->Sid(),weight,ADD_VAL);
      vectotalweight->SetValue(element->vertices[1]->Sid(),weight,ADD_VAL);
      vectotalweight->SetValue(element->vertices[2]->Sid(),weight,ADD_VAL);
   }

   /*Assemble and serialize*/
   vecH->Assemble();
   vectotalweight->Assemble();
   Hserial=vecH->ToMPISerial();
   totalweight=vectotalweight->ToMPISerial();

   /*Divide for the total weight*/
   for(int i=0;i<numberofvertices;i++){
      _assert_(totalweight[i]>0);
      Hserial[i]=Hserial[i]/totalweight[i];
   }

   /*Set element inputs*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		H[0]=Hserial[element->vertices[0]->Sid()];
		H[1]=Hserial[element->vertices[1]->Sid()];
		H[2]=Hserial[element->vertices[2]->Sid()];
		element->AddInput(ThicknessEnum,H,P1Enum);
	}

 	/*Cleanup*/
   delete vecH;
   delete vectotalweight;
   xDelete<IssmDouble>(H);
   xDelete<IssmDouble>(Hserial);
   xDelete<IssmDouble>(totalweight);
}
/*}}}*/
void FemModel::ThicknessPositivex(IssmDouble* pJ){/*{{{*/

	/*output: */
	IssmDouble J=0.;
	IssmDouble J_sum;

	IssmDouble  thickness,weight;
	IssmDouble  Jdet;
	IssmDouble* xyz_list = NULL;
	IssmDouble  H;

	/*Compute Misfit: */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*If on water, return 0: */
		if(!element->IsIceInElement()) continue;

		/* Get node coordinates*/
		element->GetVerticesCoordinates(&xyz_list);

		/*Retrieve all inputs we will be needing: */
		DatasetInput* weights_input   =element->GetDatasetInput(InversionCostFunctionsCoefficientsEnum); _assert_(weights_input);
		Input* thickness_input =element->GetInput(ThicknessEnum);                          _assert_(thickness_input);

		/* Start  looping on the number of gaussian points: */
		Gauss* gauss=element->NewGauss(2);
		while(gauss->next()){

			/* Get Jacobian determinant: */
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);

			/*Get all parameters at gaussian point*/
			weights_input->GetInputValue(&weight,gauss,ThicknessPositiveEnum);
			thickness_input->GetInputValue(&H,gauss);

			/*int min(H,0)^2 */
			if(H<=0){
				J+=weight*H*H*Jdet*gauss->weight;
			}
		}

		/*clean up and Return: */
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Sum all J from all cpus of the cluster:*/
	ISSM_MPI_Reduce (&J,&J_sum,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&J_sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	J=J_sum;

	/*Assign output pointers: */
	*pJ=J;
}
/*}}}*/
void FemModel::TimeAdaptx(IssmDouble* pdt){/*{{{*/

	/*intermediary: */
	IssmDouble   min_dt      = 1.e+50;
	IssmDouble   node_min_dt = 0.;

	/*Go through elements, and figure out the minimum of the time steps for each element (using CFL criterion): */
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		if(element->IsIceInElement()){ /*verify if there is ice in the element*/
			IssmDouble dt=element->TimeAdapt();
			if(dt<min_dt)min_dt=dt;
		}
	}

	/*Figure out minimum across the cluster: */
	ISSM_MPI_Reduce (&min_dt,&node_min_dt,1,ISSM_MPI_DOUBLE,ISSM_MPI_MIN,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&node_min_dt,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	min_dt=node_min_dt;

	/*Constrain dt */
	IssmDouble dt_low,dt_high;
	parameters->FindParam(&dt_low,TimesteppingTimeStepMinEnum);
	parameters->FindParam(&dt_high,TimesteppingTimeStepMaxEnum);
	if(min_dt<dt_low)  min_dt = dt_low;
	if(min_dt>dt_high) min_dt = dt_high;

	/*Assign output pointers:*/
	*pdt=min_dt;
}/*}}}*/
void FemModel::TotalCalvingFluxLevelsetx(IssmDouble* pM, bool scaled){/*{{{*/

	IssmDouble local_calving_flux = 0;
	IssmDouble total_calving_flux;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_calving_flux+=element->TotalCalvingFluxLevelset(scaled);
	}
	ISSM_MPI_Reduce(&local_calving_flux,&total_calving_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_calving_flux,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_calving_flux;

}/*}}}*/
void FemModel::TotalCalvingMeltingFluxLevelsetx(IssmDouble* pM, bool scaled){/*{{{*/

	IssmDouble local_calving_flux = 0;
	IssmDouble total_calving_flux;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_calving_flux+=element->TotalCalvingMeltingFluxLevelset(scaled);
	}
	ISSM_MPI_Reduce(&local_calving_flux,&total_calving_flux,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_calving_flux,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pM=total_calving_flux;

}/*}}}*/
void FemModel::TotalFloatingBmbx(IssmDouble* pFbmb, bool scaled){/*{{{*/

	IssmDouble local_fbmb = 0;
	IssmDouble total_fbmb;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_fbmb+=element->TotalFloatingBmb(scaled);
	}
	ISSM_MPI_Reduce(&local_fbmb,&total_fbmb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_fbmb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pFbmb=total_fbmb;

}/*}}}*/
void FemModel::TotalGroundedBmbx(IssmDouble* pGbmb, bool scaled){/*{{{*/

	IssmDouble local_gbmb = 0;
	IssmDouble total_gbmb;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_gbmb+=element->TotalGroundedBmb(scaled);
	}
	ISSM_MPI_Reduce(&local_gbmb,&total_gbmb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_gbmb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pGbmb=total_gbmb;

}/*}}}*/
void FemModel::TotalSmbx(IssmDouble* pSmb, bool scaled){/*{{{*/

	IssmDouble local_smb = 0;
	IssmDouble total_smb;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_smb+=element->TotalSmb(scaled);
	}
	ISSM_MPI_Reduce(&local_smb,&total_smb,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_smb,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pSmb=total_smb;

}/*}}}*/
void FemModel::TotalSmbMeltx(IssmDouble* pSmbMelt, bool scaled){/*{{{*/

	IssmDouble local_smbmelt = 0;
	IssmDouble total_smbmelt;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_smbmelt+=element->TotalSmbMelt(scaled);
	}
	ISSM_MPI_Reduce(&local_smbmelt,&total_smbmelt,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_smbmelt,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pSmbMelt=total_smbmelt;

}/*}}}*/
void FemModel::TotalSmbRefreezex(IssmDouble* pSmbRefreeze, bool scaled){/*{{{*/

	IssmDouble local_smbrefreeze = 0;
	IssmDouble total_smbrefreeze;

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		local_smbrefreeze+=element->TotalSmbRefreeze(scaled);
	}
	ISSM_MPI_Reduce(&local_smbrefreeze,&total_smbrefreeze,1,ISSM_MPI_DOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&total_smbrefreeze,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pSmbRefreeze=total_smbrefreeze;

}/*}}}*/
void FemModel::UpdateConstraintsExtrudeFromBasex(void){ /*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->UpdateConstraintsExtrudeFromBase();
	}

}
/*}}}*/
void FemModel::UpdateConstraintsExtrudeFromTopx(void){ /*{{{*/

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->UpdateConstraintsExtrudeFromTop();
	}

}
/*}}}*/
void FemModel::UpdateConstraintsx(void){ /*{{{*/

	IssmDouble time,yts;
	int        analysis_type,config_type;

	/*retrieve parameters: */
	parameters->FindParam(&analysis_type,AnalysisTypeEnum);
	parameters->FindParam(&config_type,ConfigurationTypeEnum);
	parameters->FindParam(&time,TimeEnum);
	parameters->FindParam(&yts,ConstantsYtsEnum);

	int index=AnalysisIndex(config_type);
	_assert_(this->analysis_type_list[index]==config_type);

	/*start module: */
	if(VerboseModule()) _printf0_("   Updating constraints and active domain of analysis " << EnumToStringx(analysis_type)  << " for time: " << time/yts << "\n");

	Analysis* analysis= EnumToAnalysis(analysis_type);
	analysis->UpdateConstraints(this);
	delete analysis;

	/*Second, constraints might be time dependent: */
	SpcNodesx(nodes,constraints,parameters);

	/*Now, update degrees of freedoms: */
	NodesDofx(nodes,parameters);

	/*Update FileInputs if need be*/
	if(this->inputs->IsFileInputUpdate(time)){
		_error_("not implemented yet");
	}

}/*}}}*/
int  FemModel::UpdateVertexPositionsx(void){ /*{{{*/

	IssmDouble         *surface = NULL;
	IssmDouble         *bed     = NULL;

	if(VerboseSolution()) _printf0_("   updating vertices positions\n");

	/*get vertex vectors for bed and thickness: */
	GetVectorFromInputsx(&surface  ,this, SurfaceEnum,VertexPIdEnum);
	GetVectorFromInputsx(&bed      ,this, BaseEnum,   VertexPIdEnum);

	/*Allocate vector*/
	int numvert       = vertices->NumberOfVertices();
	int numvert_local = vertices->NumberOfVerticesLocal();
	Vector<IssmDouble> *vx=new Vector<IssmDouble>(numvert_local,numvert);
	Vector<IssmDouble> *vy=new Vector<IssmDouble>(numvert_local,numvert);
	Vector<IssmDouble> *vz=new Vector<IssmDouble>(numvert_local,numvert);

	/*Update verices new geometry: */
	for(Object* & object : this->vertices->objects){
		Vertex* vertex = xDynamicCast<Vertex*>(object);
		vertex->UpdatePosition(vx,vy,vz,parameters,surface,bed);
	}

	/*Assemble mesh velocity*/
	vx->Assemble();
	vy->Assemble();
	vz->Assemble();

	/*Update element inputs*/
	InputUpdateFromVectorx(this,vx,VxMeshEnum,VertexPIdEnum);
	InputUpdateFromVectorx(this,vy,VyMeshEnum,VertexPIdEnum);
	InputUpdateFromVectorx(this,vz,VzMeshEnum,VertexPIdEnum);

	/*Free resources:*/
	delete vx;
	delete vy;
	delete vz;
	xDelete<IssmDouble>(bed);
	xDelete<IssmDouble>(surface);
	return 1;
}
/*}}}*/

/*AMR*/
#ifndef _HAVE_AD_
void FemModel::ReMesh(void){/*{{{*/

	/*Intermediaries*/
	IssmDouble *newx			= NULL;
	IssmDouble *newy			= NULL;
	IssmDouble *newz			= NULL;
	int *newelementslist		= NULL;
	int newnumberofvertices	= -1;
	int newnumberofelements = -1;

	int elementswidth       = this->GetElementsWidth();//just tria elements in this version
	int amrtype,basalforcing_model;
	bool isgroundingline;

	/*Branch to specific amr depending on requested method*/
	this->parameters->FindParam(&amrtype,AmrTypeEnum);
	switch(amrtype){
		#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
		case AmrNeopzEnum: this->ReMeshNeopz(&newnumberofvertices,&newnumberofelements,&newx,&newy,&newz,&newelementslist); break;
		#endif

		#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
		case AmrBamgEnum: this->ReMeshBamg(&newnumberofvertices,&newnumberofelements,&newx,&newy,&newz,&newelementslist); break;
		#endif

		default: _error_("not implemented yet");
	}

	/*Create iomodel for model processing*/
	IoModel* iomodel = new IoModel();
	this->parameters->FindParam(&iomodel->domaintype,DomainTypeEnum);
	this->parameters->FindParam(&iomodel->domaindim ,DomainDimensionEnum);
	this->parameters->FindParam(&iomodel->meshelementtype,MeshElementtypeEnum);
	iomodel->numberofvertices = newnumberofvertices;
	iomodel->numberofelements = newnumberofelements;
	iomodel->elements         = newelementslist;
	iomodel->AddConstant(new IoConstant(0,"md.rifts.numrifts"));
	iomodel->AddConstant(new IoConstant(false,"md.transient.isoceancoupling"));
	bool temp; int tempint;
	this->parameters->FindParam(&temp,FlowequationIsSIAEnum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isSIA"));
	this->parameters->FindParam(&temp,FlowequationIsSSAEnum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isSSA"));
	this->parameters->FindParam(&temp,FlowequationIsL1L2Enum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isL1L2"));
	this->parameters->FindParam(&temp,FlowequationIsMOLHOEnum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isMOLHO"));
	this->parameters->FindParam(&temp,FlowequationIsHOEnum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isHO"));
	this->parameters->FindParam(&temp,FlowequationIsFSEnum); iomodel->AddConstant(new IoConstant(temp,"md.flowequation.isFS"));
	this->parameters->FindParam(&tempint,MasstransportStabilizationEnum); iomodel->AddConstant(new IoConstant(tempint,"md.masstransport.stabilization"));
	iomodel->AddConstant(new IoConstant(P1Enum,"md.flowequation.fe_SSA"));

	/*Partitioning the new mesh. Maybe ElementsAndVerticesPartitioning.cpp could be modified to set this without iomodel.*/
	::ElementsAndVerticesPartitioning(iomodel);

	/*Creating elements*/
	/*Just Tria in this version*/
	Elements* new_elements=new Elements();
	this->CreateElements(newnumberofelements,elementswidth,newelementslist,iomodel->my_elements,new_elements);

	/*Create vertices*/
	Vertices* new_vertices=new Vertices();
	CreateNumberNodeToElementConnectivity(iomodel);
	::CreateVertices(new_elements,new_vertices,iomodel,TransientSolutionEnum,true);
	for(Object* & object : new_vertices->objects){
		Vertex *vertex=(Vertex*)object;
		int     sid = vertex->Sid();
		vertex->x=newx[sid];
		vertex->y=newy[sid];
		vertex->z=newz[sid];
	}

	/*Creating inputs*/
	Inputs* new_inputs=new Inputs(newnumberofelements,newnumberofvertices);

	/*Creating materials*/
	Materials* new_materials=new Materials();
	this->CreateMaterials(newnumberofelements,iomodel->my_elements,new_materials);

	/*Creating nodes and constraints*/
	/*Just SSA (2D) and P1 in this version*/
	Constraints **new_constraints_list = xNew<Constraints*>(this->nummodels);
	Nodes       **new_nodes_list       = xNew<Nodes*>(this->nummodels);

	this->analysis_counter=-1;
	for(int i=0;i<this->nummodels;i++){//create nodes for each analysis in analysis_type_list

		int analysis_enum = this->analysis_type_list[i];
		if(VerboseMProcessor()) _printf0_("   creating datasets for analysis " << EnumToStringx(analysis_enum) << "\n");

		if(this->loads_list[i]->Size()!=0) _error_("not supported yet");
		new_constraints_list[i] = new Constraints();
		new_nodes_list[i] = new Nodes();

		/*As the domain is 2D, it is not necessary to create nodes for this analysis*/
		if(analysis_enum==StressbalanceVerticalAnalysisEnum) continue;
		Analysis* analysis = EnumToAnalysis(analysis_enum);
		analysis->CreateNodes(new_nodes_list[i],iomodel,true);
		delete analysis;
		this->UpdateElements(newnumberofelements,newelementslist,iomodel->my_elements,i,new_elements);
		this->CreateConstraints(new_vertices,analysis_enum,new_constraints_list[i]);

		new_constraints_list[i]->Presort();
		new_nodes_list[i]->Presort();
	}

	new_elements->Presort();
	new_vertices->Presort();
	//this->loads->Presort();
	new_materials->Presort();

	/*reset hooks*/
	new_elements->ResetHooks();
	//this->loads->ResetHooks();
	new_materials->ResetHooks();

	/*do the post-processing of the datasets to get an FemModel that can actually run analyses: */
	int analysis_type;
	for(int i=0;i<this->nummodels;i++){
		analysis_type=this->analysis_type_list[i];
		SetCurrentConfiguration(analysis_type);

		this->analysis_counter=i;
		/*Now, plug analysis_counter and analysis_type inside the parameters: */
		this->parameters->SetParam(this->analysis_counter,AnalysisCounterEnum);
		this->parameters->SetParam(analysis_type,AnalysisTypeEnum);
		this->parameters->SetParam(analysis_type,ConfigurationTypeEnum);

		/*configure elements, loads and nodes, for this new analysis: */
		new_elements->SetCurrentConfiguration(new_elements,this->loads,new_nodes_list[i],new_vertices,new_materials,this->parameters);
		this->loads->SetCurrentConfiguration(new_elements,this->loads,new_nodes_list[i],new_vertices,new_materials,this->parameters);

		/*take care of toolkits options, that depend on this analysis type (present only after model processor)*/
		if(this->parameters->Exist(ToolkitsOptionsStringsEnum)){
			ToolkitsOptionsFromAnalysis(this->parameters,analysis_type);
			if(VerboseSolver()) _printf0_("      toolkits Options set for analysis type: " << EnumToStringx(analysis_type) << "\n");
		}

		ConfigureObjectsx(new_elements,this->loads,new_nodes_list[i],new_vertices,new_materials,this->parameters,new_inputs);
		SpcNodesx(new_nodes_list[i],new_constraints_list[i],this->parameters);
		NodesDofx(new_nodes_list[i],this->parameters);
	}

	/*Interpolate all inputs and insert them into the new elements.*/
	this->InterpolateInputs(new_vertices,new_elements,new_inputs);

	/*Delete old structure and set new pointers*/
	delete this->inputs;   this->inputs = new_inputs;
	delete this->vertices;  this->vertices = new_vertices;
	delete this->elements;  this->elements = new_elements;
	delete this->materials; this->materials = new_materials;
	if(this->constraints_list && this->nummodels){
		for(int i=0;i<this->nummodels;i++) delete this->constraints_list[i];
		xDelete<Constraints*>(this->constraints_list);
	}
	this->constraints_list= new_constraints_list;
	if(this->nodes_list && this->nummodels){
		for(int i=0;i<this->nummodels;i++) delete this->nodes_list[i];
		xDelete<Nodes*>(this->nodes_list);
	}
	this->nodes_list = new_nodes_list;

	/*Reset mask*/
	GetMaskOfIceVerticesLSMx0(this);

	/*Insert MISMIP+ bed topography FIXME it could be stay in another place*/
	this->parameters->FindParam(&basalforcing_model,BasalforcingsEnum);
	if(basalforcing_model==MismipFloatingMeltRateEnum) this->BedrockFromMismipPlus();

	/*Adjust base, thickness and mask grounded ice leve set*/
	this->parameters->FindParam(&isgroundingline,TransientIsgroundinglineEnum);
	if(isgroundingline) this->AdjustBaseThicknessAndMask();

	/*Reset current configuration: */
	analysis_type=this->analysis_type_list[this->analysis_counter];
	SetCurrentConfiguration(analysis_type);

	/*Set the new mesh*/
	this->SetMesh(&newelementslist,&newx,&newy,&newnumberofvertices,&newnumberofelements);

	/*Cleanup*/
	xDelete<IssmDouble>(newz);
	/*Delete iomodel, but make sure to not erase some pointers*/
	iomodel->elements = NULL;
	delete iomodel;
}
/*}}}*/
void FemModel::BedrockFromMismipPlus(void){/*{{{*/

	/*Insert bedrock from mismip+ setup*/
	/*This was used to Misomip project/simulations*/

	if(VerboseSolution())_printf0_("	call Mismip bedrock adjust module\n");

	IssmDouble x,y,bx,by;
	int numvertices 		= this->GetElementsWidth();
	IssmDouble* xyz_list = NULL;
   IssmDouble* r        = xNew<IssmDouble>(numvertices);

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);
		element->GetVerticesCoordinates(&xyz_list);
		for(int i=0;i<numvertices;i++){
			x		= xyz_list[3*i+0];
			y		= xyz_list[3*i+1];
			bx		= -150.-728.8*pow(x/300000.,2)+343.91*pow(x/300000.,4)-50.57*pow(x/300000.,6);
			by		= 500./(1.+exp((-2./4000.)*(y-80000./2.-24000.)))+500./(1.+exp((2./4000.)*(y-80000./2.+24000.)));
			r[i]	= max(bx+by,-720.);
		}
		/*insert new bedrock*/
		element->AddInput(BedEnum,&r[0],P1Enum);
		/*Cleanup*/
		xDelete<IssmDouble>(xyz_list);
	}
   /*Delete*/
   xDelete<IssmDouble>(r);
}
/*}}}*/
void FemModel::AdjustBaseThicknessAndMask(void){/*{{{*/

	if(VerboseSolution())_printf0_("	call adjust base and thickness module\n");

	int     numvertices = this->GetElementsWidth();
   IssmDouble rho_water,rho_ice,density,base_float;
   IssmDouble* phi     = xNew<IssmDouble>(numvertices);
   IssmDouble* h       = xNew<IssmDouble>(numvertices);
   IssmDouble* s       = xNew<IssmDouble>(numvertices);
   IssmDouble* b       = xNew<IssmDouble>(numvertices);
   IssmDouble* r       = xNew<IssmDouble>(numvertices);
   IssmDouble* sl      = xNew<IssmDouble>(numvertices);

	for(Object* & object : this->elements->objects){
      Element* element = xDynamicCast<Element*>(object);

		element->GetInputListOnVertices(&s[0],SurfaceEnum);
		element->GetInputListOnVertices(&r[0],BedEnum);
		element->GetInputListOnVertices(&sl[0],SealevelEnum);
		rho_water   = element->FindParam(MaterialsRhoSeawaterEnum);
		rho_ice     = element->FindParam(MaterialsRhoIceEnum);
		density     = rho_ice/rho_water;

		for(int i=0;i<numvertices;i++){
			/*calculate base floatation (which supports given surface*/
			base_float = rho_ice*s[i]/(rho_ice-rho_water);
			if(r[i]>base_float){
				b[i] = r[i];
			}
			else {
				b[i] = base_float;
			}

			if(fabs(sl[i])>0) _error_("Sea level value "<<sl[i]<<" not supported!");
			/*update thickness and mask grounded ice level set*/
			h[i]	  = s[i]-b[i];
			phi[i]  = h[i]+r[i]/density;
		}

		/*Update inputs*/
		element->AddInput(MaskOceanLevelsetEnum,&phi[0],P1Enum);
		element->AddInput(ThicknessEnum,&h[0],P1Enum);
		element->AddInput(BaseEnum,&b[0],P1Enum);
	}

   /*Delete*/
   xDelete<IssmDouble>(phi);
   xDelete<IssmDouble>(h);
   xDelete<IssmDouble>(s);
   xDelete<IssmDouble>(b);
   xDelete<IssmDouble>(r);
   xDelete<IssmDouble>(sl);
}
/*}}}*/
void FemModel::GetInputs(int* pnumP0inputs,IssmDouble** pP0inputs,int** pP0input_enums,int** pP0input_interp,int* pnumP1inputs,IssmDouble** pP1inputs,int** pP1input_enums,int** pP1input_interp){/*{{{*/

	int numberofvertices = this->vertices->NumberOfVertices();
	int numberofelements = this->elements->NumberOfElements();
	int elementswidth    = this->GetElementsWidth();
	int numinputs,numP0inputs,numP1inputs;
	IssmDouble* P0inputs								= NULL;
	Vector<IssmDouble>* vP0inputs					= NULL;
	int* P0input_enums								= NULL;
	int* P0input_interp 								= NULL;
	IssmDouble* P1inputs								= NULL;
	Vector<IssmDouble>* vP1inputs					= NULL;
	int* P1input_enums  								= NULL;
	int* P1input_interp 								= NULL;
	int* input_interpolations                 = NULL;
	int* input_enums                          = NULL;
   int* pos												= NULL;
	IssmDouble value									= 0;

	/*Figure out how many inputs we have and their respective interpolation*/
	this->inputs->GetInputsInterpolations(&numinputs,&input_interpolations,&input_enums);

	/*Count and get enums of all inputs in old mesh*/
	for(int step=0;step<2;step++){
		if(step){
			P0input_enums  = xNew<int>(numP0inputs);
			P0input_interp = xNew<int>(numP0inputs);
			P1input_enums  = xNew<int>(numP1inputs);
			P1input_interp = xNew<int>(numP1inputs);
		}
		numP0inputs = 0;
		numP1inputs = 0;
		for(int i=0;i<numinputs;i++){
			int inputinterp = input_interpolations[i];
			switch(inputinterp){
				case 0:
					/*Input not found, go to the next*/
					break;
				case P1Enum:
					if(step){
						P1input_enums[numP1inputs]  = input_enums[i];
						P1input_interp[numP1inputs] = inputinterp;
					}
					numP1inputs++;
					break;
				case P0Enum:
				case IntInputEnum:
				case BoolInputEnum:
					if(step){
						P0input_enums[numP0inputs]  = input_enums[i];
						P0input_interp[numP0inputs] = inputinterp;
					}
					numP0inputs++;
					break;
				default:
					_error_(EnumToStringx(inputinterp)<<" ("<<inputinterp<<") Not supported yet");
			}
		}
	}

	/*Get P0 and P1 inputs over the elements*/
	pos		= xNew<int>(elementswidth);
	vP0inputs= new Vector<IssmDouble>(numberofelements*numP0inputs);
	vP1inputs= new Vector<IssmDouble>(numberofvertices*numP1inputs);
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);

		/*Get P0 inputs*/
		for(int j=0;j<numP0inputs;j++){
			switch(P0input_interp[j]){
				case P0Enum:{
					Input* input=element->GetInput(P0input_enums[j]);
					input->GetInputAverage(&value);
				}
							 break;
				case IntInputEnum:{
					int valueint;
					element->GetInputValue(&valueint,P0input_enums[j]);
					value = reCast<IssmDouble>(valueint);
				}
									 break;
				case BoolInputEnum:{
					bool valuebool;
					element->GetInputValue(&valuebool,P0input_enums[j]);
					value = reCast<IssmDouble>(valuebool);
				}
					break;
				default:
					_error_(EnumToStringx(P0input_interp[j])<<" ("<<P0input_interp[j]<<") Not supported yet");
			}
			pos[0]=element->Sid()*numP0inputs+j;
			/*Insert input in the vector*/
			vP0inputs->SetValues(1,pos,&value,INS_VAL);
		}

		/*Get P1 inputs*/
		for(int j=0;j<numP1inputs;j++){
			Input* temp = element->GetInput(P1input_enums[j]); _assert_(temp);
			ElementInput* input=xDynamicCast<ElementInput*>(temp);
			pos[0]=element->vertices[0]->Sid()*numP1inputs+j;
			pos[1]=element->vertices[1]->Sid()*numP1inputs+j;
			pos[2]=element->vertices[2]->Sid()*numP1inputs+j;
			/*Insert input in the vector*/
			vP1inputs->SetValues(elementswidth,pos,input->element_values,INS_VAL);
		}
	}

	/*Assemble and serialize*/
	vP0inputs->Assemble();
	vP1inputs->Assemble();
	P0inputs=vP0inputs->ToMPISerial();
	P1inputs=vP1inputs->ToMPISerial();

	/*Assign pointers*/
	*pnumP0inputs		= numP0inputs;
	*pP0inputs			= P0inputs;
	*pP0input_enums	= P0input_enums;
	*pP0input_interp	= P0input_interp;
	*pnumP1inputs		= numP1inputs;
	*pP1inputs			= P1inputs;
	*pP1input_enums	= P1input_enums;
	*pP1input_interp	= P1input_interp;

	/*Cleanup*/
	delete vP0inputs;
	delete vP1inputs;
	xDelete<int>(input_interpolations);
	xDelete<int>(input_enums);
	xDelete<int>(pos);
}
/*}}}*/
void FemModel::InterpolateInputs(Vertices* newfemmodel_vertices,Elements* newfemmodel_elements,Inputs* newinputs){/*{{{*/

	int numberofelements			= -1;												//global, entire old mesh
	int newnumberofelements		= newfemmodel_elements->Size();			//just on the new partition
	int numberofvertices			= -1;												//global, entire old mesh
	int newnumberofvertices 	= newfemmodel_vertices->Size();			//just on the new partition
	int elementswidth				= this->GetElementsWidth(); //just tria in this version
	int numP0inputs				= -1;
	IssmDouble* P0inputs			= NULL; //global, entire old mesh
	IssmDouble* newP0inputs		= NULL; //just on the new partition
	int* P0input_enums			= NULL;
	int* P0input_interp 			= NULL;
	int numP1inputs				= -1;
	IssmDouble* P1inputs			= NULL; //global, entire old mesh
	IssmDouble* newP1inputs 	= NULL; //just on the new partition
	int* P1input_enums  			= NULL;
	int* P1input_interp 			= NULL;
	IssmDouble* values			= NULL;
   IssmDouble* vector      	= NULL;
	IssmDouble* x					= NULL;//global, entire old mesh
	IssmDouble* y					= NULL;//global, entire old mesh
	int* elementslist				= NULL;//global, entire old mesh
	IssmDouble* newx				= NULL;//just on the new partition
	IssmDouble* newy				= NULL;//just on the new partition
	IssmDouble* newz				= NULL;//just on the new partition
	IssmDouble* newxc				= NULL;//just on the new partition
	IssmDouble* newyc				= NULL;//just on the new partition
	int* newelementslist			= NULL;//just on the new partition
	int* sidtoindex				= NULL;//global vertices sid to partition index

	/*Get old P0 and P1  inputs (entire mesh)*/
	this->GetInputs(&numP0inputs,&P0inputs,&P0input_enums,&P0input_interp,&numP1inputs,&P1inputs,&P1input_enums,&P1input_interp);

	/*Get the old mesh (global, entire mesh)*/
	this->GetMesh(&elementslist,&x,&y,&numberofvertices,&numberofelements);

	/*Get the new mesh (just on the new partition)*/
	this->GetMeshOnPartition(newfemmodel_vertices,newfemmodel_elements,&newx,&newy,&newz,&newelementslist,&sidtoindex);

	/*Calculate the center points xc and xy (new mesh, new partition)*/
	newxc=xNewZeroInit<IssmDouble>(newnumberofelements);
	newyc=xNewZeroInit<IssmDouble>(newnumberofelements);
	for(int i=0;i<newnumberofelements;i++){
		for(int j=0;j<elementswidth;j++){
			int vid = newelementslist[i*elementswidth+j]-1;//Transform to C indexing
			newxc[i]+=newx[vid]/elementswidth;
			newyc[i]+=newy[vid]/elementswidth;
		}
	}

	/*Interplate P0 inputs in the new mesh (just on the new partition)*/
	InterpFromMeshToMesh2dx(&newP0inputs,elementslist,x,y,numberofvertices,numberofelements,
				P0inputs,numberofelements,numP0inputs,
				newxc,newyc,newnumberofelements,NULL);

	/*Interpolate P1 inputs in the new mesh (just on the new partition)*/
	InterpFromMeshToMesh2dx(&newP1inputs,elementslist,x,y,numberofvertices,numberofelements,
				P1inputs,numberofvertices,numP1inputs,
				newx,newy,newnumberofvertices,NULL);

	/*Insert P0 and P1 inputs into the new elements (just on the new partition)*/
	int vertexlids[3];
	values=xNew<IssmDouble>(elementswidth);
	for(int i=0;i<newfemmodel_elements->Size();i++){//just on the new partition
		Element* element=xDynamicCast<Element*>(newfemmodel_elements->GetObjectByOffset(i));
		/*newP0inputs is just on the new partition*/
		for(int j=0;j<numP0inputs;j++){
			switch(P0input_interp[j]){
				case P0Enum:
					element->SetElementInput(newinputs,P0input_enums[j],newP0inputs[i*numP0inputs+j]);
					break;
				case IntInputEnum:
					element->SetIntInput(newinputs,P0input_enums[j],reCast<int>(newP0inputs[i*numP0inputs+j]));
					break;
				case BoolInputEnum:
					element->SetBoolInput(newinputs,P0input_enums[j],reCast<bool>(newP0inputs[i*numP0inputs+j]));
					break;
				default:
					_error_(EnumToStringx(P0input_interp[j])<<" Not supported yet");
			}
		}
		/*newP1inputs is just on the new partition*/
		for(int i=0;i<3;i++) vertexlids[i]=element->vertices[i]->lid;
		for(int j=0;j<numP1inputs;j++){
			values[0]=newP1inputs[sidtoindex[element->vertices[0]->Sid()]*numP1inputs+j];
			values[1]=newP1inputs[sidtoindex[element->vertices[1]->Sid()]*numP1inputs+j];
			values[2]=newP1inputs[sidtoindex[element->vertices[2]->Sid()]*numP1inputs+j];
			newinputs->SetTriaInput(P1input_enums[j],P1Enum,3,vertexlids,values);
		}
	}

	/*Cleanup*/
	xDelete<IssmDouble>(P0inputs);
	xDelete<IssmDouble>(newP0inputs);
	xDelete<int>(P0input_enums);
	xDelete<int>(P0input_interp);
	xDelete<IssmDouble>(P1inputs);
	xDelete<IssmDouble>(newP1inputs);
	xDelete<int>(P1input_enums);
	xDelete<int>(P1input_interp);
	xDelete<IssmDouble>(newx);
	xDelete<IssmDouble>(newy);
	xDelete<IssmDouble>(newz);
	xDelete<IssmDouble>(newxc);
	xDelete<IssmDouble>(newyc);
	xDelete<int>(newelementslist);
	xDelete<int>(sidtoindex);
	xDelete<IssmDouble>(values);
}
/*}}}*/
void FemModel::WriteMeshInResults(void){/*{{{*/

	/*Write the erros estimators*/
	this->WriteErrorEstimatorsInResults();

	int step					= -1;
	int numberofelements = -1;
	int numberofvertices = -1;
	IssmDouble time		= -1;
	IssmDouble* x			= NULL;
	IssmDouble* y			= NULL;
	int* elementslist		= NULL;

	if(!this->elements || !this->vertices || !this->results || !this->parameters) return;

	parameters->FindParam(&step,StepEnum);
	parameters->FindParam(&time,TimeEnum);

	/*Get mesh. Elementslist comes in Matlab indexing*/
	this->GetMesh(&elementslist,&x,&y,&numberofvertices,&numberofelements);

	/*Write mesh in Results*/
	this->results->AddResult(new GenericExternalResult<int*>(this->results->Size()+1,MeshElementsEnum,
					elementslist,numberofelements,this->GetElementsWidth(),step,time));

	this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,MeshXEnum,
					x,numberofvertices,1,step,time));

	this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,MeshYEnum,
					y,numberofvertices,1,step,time));
}
/*}}}*/
void FemModel::WriteErrorEstimatorsInResults(void){/*{{{*/

   int step                   = -1;
   int numberofelements       = -1;
   IssmDouble time            = -1;
   IssmDouble* stresserror    = NULL;
   IssmDouble* thicknesserror = NULL;

   if(!this->elements || !this->vertices || !this->results || !this->parameters) return;

   parameters->FindParam(&step,StepEnum);
   parameters->FindParam(&time,TimeEnum);
   numberofelements=this->elements->NumberOfElements();

   /*Compute the deviatoric stress tensor*/
   this->ZZErrorEstimator(&stresserror);

   /*Compute the thickness error*/
   this->ThicknessZZErrorEstimator(&thicknesserror);

   /*Write error estimators in Results*/
   this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,DeviatoricStressErrorEstimatorEnum,
                                                                  stresserror,numberofelements,1,step,time));

   this->results->AddResult(new GenericExternalResult<IssmDouble*>(this->results->Size()+1,ThicknessErrorEstimatorEnum,
                                                                  thicknesserror,numberofelements,1,step,time));
   /*Cleanup*/
   xDelete<IssmDouble>(stresserror);
   xDelete<IssmDouble>(thicknesserror);

   return;
}
/*}}}*/
void FemModel::CreateElements(int newnumberofelements,int elementswidth,int* newelementslist,bool* my_elements,Elements* elements){/*{{{*/

	/*newlementslist is in Matlab indexing*/
	int lid=0;
	for(int i=0;i<newnumberofelements;i++){
		if(my_elements[i]){
			/*Create element - just tria in this version*/
			Tria *newtria=new Tria();
			newtria->id=i+1;
			newtria->sid=i;
			newtria->lid=lid++;
			newtria->iscollapsed=0;
			newtria->isonsurface = true;
			newtria->isonbase = true;
			newtria->parameters=NULL;
			newtria->inputs=NULL;
			newtria->nodes=NULL;
			newtria->vertices=NULL;
			newtria->material=NULL;
			if(this->nummodels>0){
				newtria->element_type_list=xNew<int>(this->nummodels);
				for(int j=0;j<nummodels;j++) newtria->element_type_list[j]=0;
			}
			else newtria->element_type_list=NULL;

			/*Element hook*/
			int material_id=i+1; // retrieve material_id = i+1;
			/*retrieve vertices ids*/
			int* vertex_ids=xNew<int>(elementswidth);
			for(int j=0;j<elementswidth;j++)	vertex_ids[j]=reCast<int>(newelementslist[elementswidth*i+j]);//this Hook wants Matlab indexing
			/*Setting the hooks*/
			newtria->numanalyses =this->nummodels;
			newtria->hnodes		=new Hook*[this->nummodels];
			newtria->hvertices   =new Hook(&vertex_ids[0],elementswidth);
			newtria->hmaterial   =new Hook(&material_id,1);
			newtria->hneighbors  =NULL;
			/*Initialize hnodes as NULL*/
			for(int j=0;j<this->nummodels;j++) newtria->hnodes[j]=NULL;
			/*Clean up*/
			xDelete<int>(vertex_ids);
			elements->AddObject(newtria);
		}
	}

}
/*}}}*/
void FemModel::CreateMaterials(int newnumberofelements,bool* my_elements,Materials* materials){/*{{{*/

	/*Just Matice in this version*/
	for(int i=0;i<newnumberofelements;i++){
		if(my_elements[i]){
			materials->AddObject(new Matice(i+1,i,MaticeEnum));
		}
	}
}
/*}}}*/
void FemModel::GetMesh(Vertices* femmodel_vertices, Elements* femmodel_elements,IssmDouble** px, IssmDouble** py, int** pelementslist){/*{{{*/

	if(!femmodel_vertices) _error_("GetMesh: vertices are NULL.");
	if(!femmodel_elements) _error_("GetMesh: elements are NULL.");

	int numberofvertices = femmodel_vertices->NumberOfVertices();
	int numberofelements = femmodel_elements->NumberOfElements();
	int elementswidth		= this->GetElementsWidth(); // just 2D mesh in this version (just tria elements)
	IssmDouble* x			= NULL;
	IssmDouble* y			= NULL;
	IssmDouble* z			= NULL;
	int* elementslist 	= NULL;
	int* elem_vertices	= NULL;
	IssmDouble *id1		= NULL;
   IssmDouble *id2 		= NULL;
	IssmDouble *id3 		= NULL;

	/*Get vertices coordinates*/
	VertexCoordinatesx(&x,&y,&z,femmodel_vertices,false) ;

	/*Get element vertices*/
	elem_vertices				= xNew<int>(elementswidth);
	Vector<IssmDouble>* vid1= new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vid2= new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vid3= new Vector<IssmDouble>(numberofelements);

	/*Go through elements, and for each element, get vertices*/
   for(Object* & object : femmodel_elements->objects){
    	Element* element=xDynamicCast<Element*>(object);
    	element->GetVerticesSidList(elem_vertices);
    	vid1->SetValue(element->sid,elem_vertices[0],INS_VAL);
    	vid2->SetValue(element->sid,elem_vertices[1],INS_VAL);
    	vid3->SetValue(element->sid,elem_vertices[2],INS_VAL);
   }

	/*Assemble*/
   vid1->Assemble();
   vid2->Assemble();
   vid3->Assemble();

   /*Serialize*/
	id1 = vid1->ToMPISerial();
   id2 = vid2->ToMPISerial();
	id3 = vid3->ToMPISerial();

	/*Construct elements list*/
	elementslist=xNew<int>(numberofelements*elementswidth);
	if(numberofelements*elementswidth<0) _error_("numberofelements negative.");
	for(int i=0;i<numberofelements;i++){
		elementslist[elementswidth*i+0] = reCast<int>(id1[i])+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+1] = reCast<int>(id2[i])+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+2] = reCast<int>(id3[i])+1; //InterpMesh wants Matlab indexing
	}

	/*Assign pointers*/
	*px				= x;
	*py				= y;
	*pelementslist = elementslist; //Matlab indexing. InterMesh uses this type.

	/*Cleanup*/
	xDelete<int>(elem_vertices);
	xDelete<IssmDouble>(id1);
	xDelete<IssmDouble>(id2);
	xDelete<IssmDouble>(id3);
	xDelete<IssmDouble>(z);
	delete vid1;
	delete vid2;
	delete vid3;
}
/*}}}*/
void FemModel::GetMesh(int** elementslist, IssmDouble** x, IssmDouble** y, int* numberofvertices, int* numberofelements){/*{{{*/

	int amrtype;
	this->parameters->FindParam(&amrtype,AmrTypeEnum);

	switch(amrtype){

      #if defined(_HAVE_NEOPZ_)
      case AmrNeopzEnum: this->amr->GetMesh(elementslist,x,y,numberofvertices,numberofelements); break;
      #endif

      #if defined(_HAVE_BAMG_)
      case AmrBamgEnum: this->amrbamg->GetMesh(elementslist,x,y,numberofvertices,numberofelements); break;
      #endif

      default: _error_("not implemented yet");
   }
}/*}}}*/
void FemModel::SetMesh(int** elementslist, IssmDouble** x, IssmDouble** y, int* numberofvertices, int* numberofelements){/*{{{*/

	int amrtype;
	this->parameters->FindParam(&amrtype,AmrTypeEnum);

	switch(amrtype){

      #if defined(_HAVE_NEOPZ_)
      case AmrNeopzEnum: this->amr->SetMesh(elementslist,x,y,numberofvertices,numberofelements); break;
      #endif

      #if defined(_HAVE_BAMG_)
      case AmrBamgEnum: this->amrbamg->SetMesh(elementslist,x,y,numberofvertices,numberofelements); break;
      #endif

      default: _error_("not implemented yet");
   }
}/*}}}*/
void FemModel::GetMeshOnPartition(Vertices* femmodel_vertices,Elements* femmodel_elements,IssmDouble** px,IssmDouble** py,IssmDouble** pz,int** pelementslist,int** psidtoindex){/*{{{*/

	if(!femmodel_vertices) _error_("GetMesh: vertices are NULL.");
	if(!femmodel_elements) _error_("GetMesh: elements are NULL.");

	int numberofvertices			= femmodel_vertices->Size();	//number of vertices of this partition
	int numbertotalofvertices	= femmodel_vertices->NumberOfVertices();	//number total of vertices (entire mesh)
	int numberofelements			= femmodel_elements->Size();  //number of elements of this partition
	int elementswidth				= this->GetElementsWidth();	//just 2D mesh in this version (just tria elements)
	IssmDouble* x					= NULL;
	IssmDouble* y					= NULL;
	IssmDouble* z					= NULL;
	int* elementslist				= NULL;
	int* sidtoindex				= NULL;
	int* elem_vertices			= NULL;

	/*Get vertices coordinates of this partition*/
	sidtoindex	= xNewZeroInit<int>(numbertotalofvertices);//entire mesh, all vertices
	x				= xNew<IssmDouble>(numberofvertices);//just this partition
	y				= xNew<IssmDouble>(numberofvertices);//just this partitio;
	z				= xNew<IssmDouble>(numberofvertices);//just this partitio;

	/*Go through in this partition (vertices)*/
	for(int i=0;i<numberofvertices;i++){//just this partition
		Vertex* vertex=(Vertex*)femmodel_vertices->GetObjectByOffset(i);
		/*Attention: no spherical coordinates*/
		x[i]=vertex->GetX();
		y[i]=vertex->GetY();
		z[i]=vertex->GetZ();
		/*Keep the index and sid pair*/
		sidtoindex[vertex->Sid()]=i;
	}

	/*Go through in this partition (elements) and build the element list*/
	elem_vertices= xNew<int>(elementswidth);
	elementslist = xNew<int>(numberofelements*elementswidth);
	if(numberofelements*elementswidth<0) _error_("numberofelements negative.");

	for(int i=0;i<numberofelements;i++){//just this partition
    	Element* element=xDynamicCast<Element*>(femmodel_elements->GetObjectByOffset(i));
    	element->GetVerticesSidList(elem_vertices);
		elementslist[elementswidth*i+0] = sidtoindex[elem_vertices[0]]+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+1] = sidtoindex[elem_vertices[1]]+1; //InterpMesh wants Matlab indexing
		elementslist[elementswidth*i+2] = sidtoindex[elem_vertices[2]]+1; //InterpMesh wants Matlab indexing
	}

	/*Assign pointers*/
	*px				= x;
	*py				= y;
	*pz				= z;
	*pelementslist = elementslist; //Matlab indexing. InterMesh uses this type.
	*psidtoindex	= sidtoindex;  //it is ncessary to insert inputs

	/*Cleanup*/
	xDelete<int>(elem_vertices);
}
/*}}}*/
void FemModel::CreateConstraints(Vertices* newfemmodel_vertices,int analysis_enum,Constraints* newfemmodel_constraints){/*{{{*/

	/*ATTENTION: JUST SPCVX AND SPCVY*/
	/*OTHERS CONSTRAINTS MUST BE IMPLEMENTED*/
	if(analysis_enum!=StressbalanceAnalysisEnum) return;
	int analysis_index = AnalysisIndex(analysis_enum);

	int numberofnodes_analysistype= this->nodes_list[analysis_index]->NumberOfNodes();
	int dofpernode						= 2;														//vx and vy
	int numberofcols					= dofpernode*2;										//to keep dofs and flags in the vspc vector
	int numberofvertices				= -1;														//global, entire old mesh
	int numberofelements				= -1;														//global, entire old mesh
	int newnumberofvertices			= newfemmodel_vertices->Size();					//local, just the new partition
	int count							= 0;
	IssmDouble* x						= NULL;													//global, entire old mesh
	IssmDouble* y						= NULL;													//global, entire old mesh
	int*			elementslist		= NULL;													//global, entire old mesh
	IssmDouble* spc					= NULL;													//global, entire old mesh
	IssmDouble* newx					= NULL;													//local, just new partition
	IssmDouble* newy					= NULL;													//local, just new partition
	IssmDouble* newspc				= NULL;													//local, just new partition
	IssmDouble eps						= 1.e-8;
	Vector<IssmDouble>* vspc		= new Vector<IssmDouble>(numberofnodes_analysistype*numberofcols);

	/*Get old mesh (global, entire mesh). Elementslist comes in Matlab indexing*/
	this->GetMesh(&elementslist,&x,&y,&numberofvertices,&numberofelements);

	/*Get vertices coordinates of the new partition*/
	newx=xNew<IssmDouble>(newnumberofvertices);//just the new partition
	newy=xNew<IssmDouble>(newnumberofvertices);//just the new partition
	for(int i=0;i<newnumberofvertices;i++){//just the new partition
		Vertex* vertex=(Vertex*)newfemmodel_vertices->GetObjectByOffset(i);
		/*Attention: no spherical coordinates*/
		newx[i]=vertex->GetX();
		newy[i]=vertex->GetY();
	}

	/*Get spcvx and spcvy of old mesh*/
	for(Object* & object : this->constraints_list[analysis_index]->objects){
		Constraint* constraint=(Constraint*)object;
		SpcStatic* spcstatic = xDynamicCast<SpcStatic*>(constraint);
		int dof					= spcstatic->GetDof();
		int node					= spcstatic->GetNodeId();
		IssmDouble spcvalue	= spcstatic->GetValue();
		int nodeindex			= node-1;

		/*vx and vx flag insertion*/
		if(dof==0) {//vx
			vspc->SetValue(nodeindex*numberofcols,spcvalue,INS_VAL);    //vx
			vspc->SetValue(nodeindex*numberofcols+dofpernode,1,INS_VAL);//vxflag
		}
		/*vy and vy flag insertion*/
		if(dof==1){//vy
			vspc->SetValue(nodeindex*numberofcols+1,spcvalue,INS_VAL);	//vy
			vspc->SetValue(nodeindex*numberofcols+dofpernode+1,1,INS_VAL);//vyflag
		}
	}

	/*Assemble and serialize*/
	vspc->Assemble();
	spc=vspc->ToMPISerial();

	/*Interpolate spc values and flags in the new partition*/
	InterpFromMeshToMesh2dx(&newspc,elementslist,x,y,numberofvertices,numberofelements,
								spc,numberofvertices,numberofcols,
								newx,newy,newnumberofvertices,NULL);

	/*Now, insert the interpolated constraints in the data set (constraints)*/
	count=0;
	for(int i=0;i<newnumberofvertices;i++){//just in the new partition
		Vertex* vertex=(Vertex*)newfemmodel_vertices->GetObjectByOffset(i);
		/*spcvx*/
		if(!xIsNan<IssmDouble>(newspc[i*numberofcols]) && newspc[i*numberofcols+dofpernode]>(1-eps)){
			newfemmodel_constraints->AddObject(new SpcStatic(count+1,vertex->Sid()+1,0,newspc[i*numberofcols],analysis_enum));
			//add count'th spc, on node i+1, setting dof 1 to vx.
			count++;
		}
	}
	count=0;
	for(int i=0;i<newnumberofvertices;i++){//just in the new partition
		Vertex* vertex=(Vertex*)newfemmodel_vertices->GetObjectByOffset(i);
		/*spcvy*/
		if(!xIsNan<IssmDouble>(newspc[i*numberofcols+1]) && newspc[i*numberofcols+dofpernode+1]>(1-eps) ){
			newfemmodel_constraints->AddObject(new SpcStatic(count+1,vertex->Sid()+1,1,newspc[i*numberofcols+1],analysis_enum));
			//add count'th spc, on node i+1, setting dof 1 to vx.
			count++;
		}
	}

	/*Cleanup*/
	xDelete<IssmDouble>(spc);
	xDelete<IssmDouble>(newspc);
	xDelete<IssmDouble>(newx);
	xDelete<IssmDouble>(newy);
	delete vspc;
}
/*}}}*/
void FemModel::UpdateElements(int newnumberofelements,int* newelementslist,bool* my_elements,int analysis_counter,Elements* newelements){/*{{{*/

	/*newelementslist is in Matlab indexing*/

	/*Update elements, set hnode.
	This code is in all analysis */
	int elemcounter=0;
	for(int iel=0;iel<newnumberofelements;iel++){
		if(my_elements[iel]){
			Tria* tria=(Tria*)newelements->GetObjectByOffset(elemcounter);
			//element update
			tria->element_type_list[analysis_counter]=P1Enum;
			int numnodes=3;
         int* tria_node_ids=xNew<int>(numnodes);
         tria_node_ids[0]=newelementslist[3*iel+0]; //matlab indexing
         tria_node_ids[1]=newelementslist[3*iel+1]; //matlab indexing
         tria_node_ids[2]=newelementslist[3*iel+2]; //matlab indexing
			tria->SetHookNodes(tria_node_ids,numnodes,analysis_counter); tria->nodes=NULL;
   		xDelete<int>(tria_node_ids);
			elemcounter++;
		}
	}
	return;
}
/*}}}*/
void FemModel::SmoothedDeviatoricStressTensor(IssmDouble** ptauxx,IssmDouble** ptauyy,IssmDouble** ptauxy){/*{{{*/

	int elementswidth							= this->GetElementsWidth();//just 2D mesh, tria elements
   int numberofvertices						= this->vertices->NumberOfVertices();
   IssmDouble weight 						= 0.;
	IssmDouble*	tauxx							= NULL;
	IssmDouble*	tauyy							= NULL;
	IssmDouble*	tauxy							= NULL;
   IssmDouble* totalweight 				= NULL;
	IssmDouble* deviatoricstressxx 		= xNew<IssmDouble>(elementswidth);
   IssmDouble* deviatoricstressyy 		= xNew<IssmDouble>(elementswidth);
   IssmDouble* deviatoricstressxy 		= xNew<IssmDouble>(elementswidth);
   int* elem_vertices 						= xNew<int>(elementswidth);
   Vector<IssmDouble>* vectauxx			= new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vectauyy			= new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vectauxy			= new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vectotalweight	= new Vector<IssmDouble>(numberofvertices);

	/*Update the Deviatoric Stress tensor over the elements*/
	this->DeviatoricStressx();

   /*Calculate the Smoothed Deviatoric Stress tensor*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      element->GetInputListOnVertices(deviatoricstressxx,DeviatoricStressxxEnum);
      element->GetInputListOnVertices(deviatoricstressyy,DeviatoricStressyyEnum);
      element->GetInputListOnVertices(deviatoricstressxy,DeviatoricStressxyEnum);
      element->GetVerticesSidList(elem_vertices);

		/*weight to calculate the smoothed deviatoric stress*/
		Tria* triaelement = xDynamicCast<Tria*>(element);
		weight				= triaelement->GetArea();//the tria area is a choice for the weight

      /*taux xx*/
		vectauxx->SetValue(elem_vertices[0],weight*deviatoricstressxx[0],ADD_VAL);
      vectauxx->SetValue(elem_vertices[1],weight*deviatoricstressxx[1],ADD_VAL);
      vectauxx->SetValue(elem_vertices[2],weight*deviatoricstressxx[2],ADD_VAL);
      /*tau yy*/
		vectauyy->SetValue(elem_vertices[0],weight*deviatoricstressyy[0],ADD_VAL);
	   vectauyy->SetValue(elem_vertices[1],weight*deviatoricstressyy[1],ADD_VAL);
      vectauyy->SetValue(elem_vertices[2],weight*deviatoricstressyy[2],ADD_VAL);
      /*tau xy*/
		vectauxy->SetValue(elem_vertices[0],weight*deviatoricstressxy[0],ADD_VAL);
      vectauxy->SetValue(elem_vertices[1],weight*deviatoricstressxy[1],ADD_VAL);
      vectauxy->SetValue(elem_vertices[2],weight*deviatoricstressxy[2],ADD_VAL);
		/*total weight*/
		vectotalweight->SetValue(elem_vertices[0],weight,ADD_VAL);
      vectotalweight->SetValue(elem_vertices[1],weight,ADD_VAL);
      vectotalweight->SetValue(elem_vertices[2],weight,ADD_VAL);
   }

   /*Assemble*/
   vectauxx->Assemble();
   vectauyy->Assemble();
   vectauxy->Assemble();
   vectotalweight->Assemble();

   /*Serialize*/
   tauxx			= vectauxx->ToMPISerial();
   tauyy			= vectauyy->ToMPISerial();
   tauxy			= vectauxy->ToMPISerial();
   totalweight	= vectotalweight->ToMPISerial();

	/*Divide for the total weight*/
	for(int i=0;i<numberofvertices;i++){
		_assert_(totalweight[i]>0);
		tauxx[i] = tauxx[i]/totalweight[i];
		tauyy[i] = tauyy[i]/totalweight[i];
		tauxy[i] = tauxy[i]/totalweight[i];
	}

	/*Set output*/
	(*ptauxx) = tauxx;
	(*ptauyy) = tauyy;
	(*ptauxy) = tauxy;

   /*Cleanup*/
   delete vectauxx;
   delete vectauyy;
   delete vectauxy;
	delete vectotalweight;
   xDelete<IssmDouble>(deviatoricstressxx);
   xDelete<IssmDouble>(deviatoricstressyy);
   xDelete<IssmDouble>(deviatoricstressxy);
   xDelete<IssmDouble>(totalweight);
   xDelete<int>(elem_vertices);
}
/*}}}*/
void FemModel::ZZErrorEstimator(IssmDouble** pelementerror){/*{{{*/

	/*Compute the Zienkiewicz and Zhu (ZZ) error estimator for the deviatoric stress tensor.
	 * Ref.: Zienkiewicz and Zhu, A Simple Error Estimator and Adaptive Procedure for Practical Engineering Analysis, Int. J. Numer. Meth. Eng, 1987*/

	IssmDouble Jdet,error,ftxx,ftyy,ftxy;
	int sid;
	int numnodes							= this->GetElementsWidth();//just 2D mesh, tria elements, P1
	int numberofelements 				= this->elements->NumberOfElements();
	IssmDouble* xyz_list 				= NULL;
	IssmDouble* smoothedtauxx			= NULL;
	IssmDouble* smoothedtauyy			= NULL;
	IssmDouble* smoothedtauxy			= NULL;
	IssmDouble* tauxx						= xNew<IssmDouble>(numnodes);
   IssmDouble* tauyy						= xNew<IssmDouble>(numnodes);
   IssmDouble* tauxy						= xNew<IssmDouble>(numnodes);
	IssmDouble* basis 					= xNew<IssmDouble>(numnodes);
	int* elem_vertices 					= xNew<int>(numnodes);
   Vector<IssmDouble>* velementerror= new Vector<IssmDouble>(numberofelements);

	/*Get smoothed deviatoric stress tensor*/
	this->SmoothedDeviatoricStressTensor(&smoothedtauxx,&smoothedtauyy,&smoothedtauxy);

	/*Integrate the error over elements*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		element->GetInputListOnVertices(tauxx,DeviatoricStressxxEnum);
      element->GetInputListOnVertices(tauyy,DeviatoricStressyyEnum);
      element->GetInputListOnVertices(tauxy,DeviatoricStressxyEnum);
      element->GetVerticesSidList(elem_vertices);

		/*Integrate*/
		element->GetVerticesCoordinates(&xyz_list);
		Gauss* gauss=element->NewGauss(2);
   	error=0.;
		while(gauss->next()){
			element->JacobianDeterminant(&Jdet,xyz_list,gauss);
			element->NodalFunctions(basis,gauss);
			ftxx=0;ftyy=0;ftxy=0;
			for(int n=0;n<numnodes;n++) {
				ftxx+=(tauxx[n]-smoothedtauxx[elem_vertices[n]])*basis[n];
				ftyy+=(tauyy[n]-smoothedtauyy[elem_vertices[n]])*basis[n];
				ftxy+=(tauxy[n]-smoothedtauxy[elem_vertices[n]])*basis[n];
			}
			error+=Jdet*gauss->weight*( pow(ftxx,2)+pow(ftyy,2)+pow(ftxy,2) ); //e^2
		}
		/*Set the error in the global vector*/
      sid=element->Sid();
		error = sqrt(error);//sqrt(e^2)
		velementerror->SetValue(sid,error,INS_VAL);
		/*Cleanup intermediaries*/
		xDelete<IssmDouble>(xyz_list);
		delete gauss;
	}

	/*Assemble*/
   velementerror->Assemble();

   /*Serialize and set output*/
   (*pelementerror)=velementerror->ToMPISerial();

	/*Cleanup*/
	xDelete<IssmDouble>(smoothedtauxx);
	xDelete<IssmDouble>(smoothedtauyy);
	xDelete<IssmDouble>(smoothedtauxy);
	xDelete<IssmDouble>(tauxx);
	xDelete<IssmDouble>(tauyy);
	xDelete<IssmDouble>(tauxy);
	xDelete<IssmDouble>(basis);
	xDelete<int>(elem_vertices);
	delete velementerror;
}
/*}}}*/
void FemModel::SmoothedGradThickness(IssmDouble** pdHdx,IssmDouble** pdHdy){/*{{{*/

   int elementswidth                   = this->GetElementsWidth();//just 2D mesh, tria elements
   int numberofvertices                = this->vertices->NumberOfVertices();

   IssmDouble weight                   = 0.;
   IssmDouble* dHdx                    = NULL;
   IssmDouble* dHdy                    = NULL;
   IssmDouble* totalweight             = NULL;
   IssmDouble* xyz_list                = NULL;
   IssmDouble* H                       = xNew<IssmDouble>(elementswidth);
   IssmDouble* GradH                   = xNew<IssmDouble>(2);
   int* elem_vertices                  = xNew<int>(elementswidth);
   Vector<IssmDouble>* vecdHdx         = new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vecdHdy         = new Vector<IssmDouble>(numberofvertices);
   Vector<IssmDouble>* vectotalweight  = new Vector<IssmDouble>(numberofvertices);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      element->GetInputListOnVertices(H,ThicknessEnum);
      element->GetVerticesSidList(elem_vertices);
      element->GetVerticesCoordinates(&xyz_list);

      /*Get the gradient of thickness at the center point (in fact, GradH is constante over the element)*/
      Gauss* gauss=element->NewGauss(1);
      gauss->GaussPoint(0);
      element->ValueP1DerivativesOnGauss(GradH,H,xyz_list,gauss);

      /*weight to calculate the smoothed grad H*/
      Tria* triaelement = xDynamicCast<Tria*>(element);
      weight            = triaelement->GetArea();//the tria area is a choice for the weight

		/*dH/dx*/
      vecdHdx->SetValue(elem_vertices[0],weight*GradH[0],ADD_VAL);
      vecdHdx->SetValue(elem_vertices[1],weight*GradH[0],ADD_VAL);
      vecdHdx->SetValue(elem_vertices[2],weight*GradH[0],ADD_VAL);
      /*dH/dy*/
      vecdHdy->SetValue(elem_vertices[0],weight*GradH[1],ADD_VAL);
      vecdHdy->SetValue(elem_vertices[1],weight*GradH[1],ADD_VAL);
      vecdHdy->SetValue(elem_vertices[2],weight*GradH[1],ADD_VAL);
      /*total weight*/
      vectotalweight->SetValue(elem_vertices[0],weight,ADD_VAL);
      vectotalweight->SetValue(elem_vertices[1],weight,ADD_VAL);
      vectotalweight->SetValue(elem_vertices[2],weight,ADD_VAL);
      /*Cleanup intermediaries*/
      xDelete<IssmDouble>(xyz_list);
      delete gauss;
   }

   /*Assemble*/
   vecdHdx->Assemble();
   vecdHdy->Assemble();
   vectotalweight->Assemble();

   /*Serialize*/
   dHdx        = vecdHdx->ToMPISerial();
   dHdy        = vecdHdy->ToMPISerial();
   totalweight = vectotalweight->ToMPISerial();

   /*Divide for the total weight*/
   for(int i=0;i<numberofvertices;i++){
      _assert_(totalweight[i]>0);
      dHdx[i] = dHdx[i]/totalweight[i];
      dHdy[i] = dHdy[i]/totalweight[i];
   }

   /*Set output*/
   (*pdHdx) = dHdx;
   (*pdHdy) = dHdy;

 	/*Cleanup*/
   delete vecdHdx;
   delete vecdHdy;
   delete vectotalweight;
   xDelete<IssmDouble>(H);
   xDelete<IssmDouble>(GradH);
   xDelete<IssmDouble>(totalweight);
   xDelete<int>(elem_vertices);
}
/*}}}*/
void FemModel::ThicknessZZErrorEstimator(IssmDouble** pelementerror){/*{{{*/
   /*Compute the Zienkiewicz and Zhu (ZZ) error estimator for the thickness
    * Ref.: Zienkiewicz and Zhu, A Simple Error Estimator and Adaptive Procedure for Practical Engineering Analysis, Int. J. Numer. Meth. Eng, 1987*/

   IssmDouble Jdet,error,fdHdx,fdHdy;
   int sid;
   int numnodes                     = this->GetElementsWidth();//just 2D mesh, tria elements, P1
   int numberofelements             = this->elements->NumberOfElements();
   IssmDouble* xyz_list             = NULL;
   IssmDouble* smoothed_dHdx        = NULL;
   IssmDouble* smoothed_dHdy        = NULL;
   IssmDouble* H                    = xNew<IssmDouble>(numnodes);
   IssmDouble* GradH                = xNew<IssmDouble>(2);
   IssmDouble* basis                = xNew<IssmDouble>(numnodes);
   int* elem_vertices               = xNew<int>(numnodes);
   Vector<IssmDouble>* velementerror= new Vector<IssmDouble>(numberofelements);

   /*Get smoothed deviatoric stress tensor*/
   this->SmoothedGradThickness(&smoothed_dHdx,&smoothed_dHdy);

	/*Integrate the error over elements*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      element->GetInputListOnVertices(H,ThicknessEnum);
      element->GetVerticesSidList(elem_vertices);
      element->GetVerticesCoordinates(&xyz_list);
      /*Get the gradient of thickness*/
      Gauss* gaussH=element->NewGauss(1);
      gaussH->GaussPoint(0);
      element->ValueP1DerivativesOnGauss(GradH,H,xyz_list,gaussH);
      /*Integrate*/
      Gauss* gauss=element->NewGauss(2);
      error=0.;
		while(gauss->next()){
         element->JacobianDeterminant(&Jdet,xyz_list,gauss);
         element->NodalFunctions(basis,gauss);
         fdHdx=0;fdHdy=0;
         for(int n=0;n<numnodes;n++) {
            fdHdx+=(GradH[0]-smoothed_dHdx[elem_vertices[n]])*basis[n];
            fdHdy+=(GradH[1]-smoothed_dHdy[elem_vertices[n]])*basis[n];
         }
         error+=Jdet*gauss->weight*(pow(fdHdx,2)+pow(fdHdy,2) ); //e^2
      }
      /*Set the error in the global vector*/
      sid=element->Sid();
		error = sqrt(error); //sqrt( e^2 )
      velementerror->SetValue(sid,error,INS_VAL);
      /*Cleanup intermediaries*/
      xDelete<IssmDouble>(xyz_list);
      delete gaussH;
      delete gauss;
   }

   /*Assemble*/
   velementerror->Assemble();

   /*Serialize and set output*/
   (*pelementerror)=velementerror->ToMPISerial();

   /*Cleanup*/
   xDelete<IssmDouble>(smoothed_dHdx);
   xDelete<IssmDouble>(smoothed_dHdy);
   xDelete<IssmDouble>(H);
   xDelete<IssmDouble>(GradH);
   xDelete<IssmDouble>(basis);
   xDelete<int>(elem_vertices);
   delete velementerror;
}
/*}}}*/
void FemModel::MeanGroundedIceLevelSet(IssmDouble** pmasklevelset){/*{{{*/

   int elementswidth                   = this->GetElementsWidth();
   int numberofelements                = this->elements->NumberOfElements();
   IssmDouble* elementlevelset         = xNew<IssmDouble>(elementswidth);
   Vector<IssmDouble>* vmasklevelset   = new Vector<IssmDouble>(numberofelements);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      element->GetInputListOnVertices(elementlevelset,MaskOceanLevelsetEnum);
      int sid = element->Sid();
      vmasklevelset->SetValue(sid,(elementlevelset[0]+elementlevelset[1]+elementlevelset[2])/3.,INS_VAL);
   }

   /*Assemble*/
   vmasklevelset->Assemble();

   /*Serialize and set output*/
   (*pmasklevelset)=vmasklevelset->ToMPISerial();

   /*Cleanup*/
   xDelete<IssmDouble>(elementlevelset);
   delete vmasklevelset;
}
/*}}}*/
void FemModel::GetElementCenterCoordinates(IssmDouble** pxc,IssmDouble** pyc){/*{{{*/

	/*Intermediaries*/
   int elementswidth          = this->GetElementsWidth();
   int numberofelements       = this->elements->NumberOfElements();
   int* elem_vertices			= xNew<int>(elementswidth);
   Vector<IssmDouble>* vxc		= new Vector<IssmDouble>(numberofelements);
   Vector<IssmDouble>* vyc		= new Vector<IssmDouble>(numberofelements);
	IssmDouble* x					= NULL;
	IssmDouble* y					= NULL;
	IssmDouble* z					= NULL;
	IssmDouble* xyz_list			= NULL;
	IssmDouble x1,y1,x2,y2,x3,y3;

	/*Insert the element center coordinates*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      //element->GetVerticesSidList(elem_vertices);
      int sid = element->Sid();
		element->GetVerticesCoordinates(&xyz_list);
		x1 = xyz_list[3*0+0];y1 = xyz_list[3*0+1];
		x2 = xyz_list[3*1+0];y2 = xyz_list[3*1+1];
		x3 = xyz_list[3*2+0];y3 = xyz_list[3*2+1];
		vxc->SetValue(sid,(x1+x2+x3)/3.,INS_VAL);
      vyc->SetValue(sid,(y1+y2+y3)/3.,INS_VAL);
   }

   /*Assemble*/
   vxc->Assemble();
   vyc->Assemble();

   /*Serialize and set output*/
   (*pxc)=vxc->ToMPISerial();
   (*pyc)=vyc->ToMPISerial();

   /*Cleanup*/
	xDelete<IssmDouble>(x);
	xDelete<IssmDouble>(y);
	xDelete<IssmDouble>(z);
	xDelete<IssmDouble>(xyz_list);
   xDelete<int>(elem_vertices);
   delete vxc;
   delete vyc;
}
/*}}}*/
void FemModel::GetZeroLevelSetPoints(IssmDouble** pzerolevelset_points,int &numberofpoints,int levelset_type){/*{{{*/

	/*Here, "zero level set" means grounding line or ice front, depending on the level set type*/
	/*pzerolevelset_points are the element center points with zero level set. X and Y coords*/
	if(levelset_type!=MaskOceanLevelsetEnum && levelset_type!=MaskIceLevelsetEnum){
		_error_("level set type not implemented yet!");
	}

	/*Outputs*/
	IssmDouble* zerolevelset_points			= NULL;
	int npoints										= 0;

	/*Intermediaries*/
 	int elementswidth                   	= this->GetElementsWidth();
   int numberofelements                	= this->elements->NumberOfElements();
	int* elem_vertices         				= xNew<int>(elementswidth);
   IssmDouble* levelset      					= xNew<IssmDouble>(elementswidth);
   IssmDouble* xyz_list							= NULL;
	Vector<IssmDouble>* vx_zerolevelset		= new Vector<IssmDouble>(numberofelements);
	Vector<IssmDouble>* vy_zerolevelset		= new Vector<IssmDouble>(numberofelements);
	IssmDouble* x_zerolevelset					= NULL;
	IssmDouble* y_zerolevelset					= NULL;
	int count,sid;
	IssmDouble xc,yc,x1,y1,x2,y2,x3,y3;

	/*Use the element center coordinate if level set is zero (grounding line or ice front), otherwise set NAN*/
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      element->GetInputListOnVertices(levelset,levelset_type);
		element->GetVerticesSidList(elem_vertices);
		sid= element->Sid();
		element->GetVerticesCoordinates(&xyz_list);
		x1 = xyz_list[3*0+0];y1 = xyz_list[3*0+1];
		x2 = xyz_list[3*1+0];y2 = xyz_list[3*1+1];
		x3 = xyz_list[3*2+0];y3 = xyz_list[3*2+1];
		xc	= NAN;
		yc	= NAN;
     	Tria* tria 	= xDynamicCast<Tria*>(element);
		if(tria->IsIceInElement()){/*verify if there is ice in the element*/
			if(levelset[0]*levelset[1]<0. || levelset[0]*levelset[2]<0. ||
				abs(levelset[0]*levelset[1])<DBL_EPSILON || abs(levelset[0]*levelset[2])<DBL_EPSILON) {
				xc=(x1+x2+x3)/3.;
				yc=(y1+y2+y3)/3.;
			}
		}
		vx_zerolevelset->SetValue(sid,xc,INS_VAL);
		vy_zerolevelset->SetValue(sid,yc,INS_VAL);
		xDelete<IssmDouble>(xyz_list);
	}
   /*Assemble and serialize*/
   vx_zerolevelset->Assemble();
   vy_zerolevelset->Assemble();
   x_zerolevelset=vx_zerolevelset->ToMPISerial();
   y_zerolevelset=vy_zerolevelset->ToMPISerial();

	/*Find the number of points*/
	npoints=0;
	for(int i=0;i<numberofelements;i++) if(!xIsNan<IssmDouble>(x_zerolevelset[i])) npoints++;

	/*Keep just the element center coordinates with zero level set (compact the structure)*/
	zerolevelset_points=xNew<IssmDouble>(2*npoints);//x and y
	count=0;
	for(int i=0;i<numberofelements;i++){
		if(!xIsNan<IssmDouble>(x_zerolevelset[i])){
			zerolevelset_points[2*count]	 = x_zerolevelset[i];
			zerolevelset_points[2*count+1] = y_zerolevelset[i];
			count++;
		}
	}

	/*Assign outputs*/
	numberofpoints				= npoints;
	(*pzerolevelset_points) = zerolevelset_points;

	/*Cleanup*/
   xDelete<int>(elem_vertices);
   xDelete<IssmDouble>(levelset);
	xDelete<IssmDouble>(x_zerolevelset);
	xDelete<IssmDouble>(y_zerolevelset);
   xDelete<IssmDouble>(xyz_list);
	delete vx_zerolevelset;
	delete vy_zerolevelset;
}
/*}}}*/
#endif

#ifdef  _HAVE_DAKOTA_
void FemModel::DakotaResponsesx(double* d_responses,char** responses_descriptors,int numresponsedescriptors,int d_numresponses){/*{{{*/

	/*intermediary: */
	int    i,j;
	char   root[50];
	int    index;
	double femmodel_response;
	int    flag;
	double *vertex_response   = NULL;
	double *qmu_response      = NULL;
	double *responses_pointer = NULL;

	IssmDouble **response_partitions         = NULL;
	IssmDouble * response_partition         = NULL;
	int * response_partitions_npart         = NULL;
	int          response_partitions_num;
	int          npart;

	/*retrieve partition vectors for responses that are scaled:*/
	this->parameters->FindParam(&response_partitions,&response_partitions_num,NULL,NULL,QmuResponsePartitionsEnum);
	this->parameters->FindParam(&response_partitions_npart,NULL,NULL,QmuResponsePartitionsNpartEnum);

	/*retrieve my_rank: */
	int my_rank=IssmComm::GetRank();

	/*save the d_responses pointer: */
	responses_pointer=d_responses;

	//watch out, we have more d_numresponses than numresponsedescriptors, because the responses have been expanded if they were scaled.
	//because we don't know the d_responses descriptors (the scaled ones) we can't key off them, so we will key off the responses_descriptors: */

	for(i=0;i<numresponsedescriptors;i++){

		flag=DescriptorIndex(root,&index,responses_descriptors[i]);

		if(flag==ScaledEnum){

			/*this response was scaled. pick up the response from the inputs: */
			GetVectorFromInputsx(&vertex_response,this, StringToEnumx(root),VertexPIdEnum);

			/*recover partition vector: */
			response_partition=response_partitions[i];
			npart=response_partitions_npart[i];

			/*Now, average it onto the partition nodes: */
			AverageOntoPartitionx(&qmu_response,elements,nodes,vertices,loads,materials,parameters,vertex_response,response_partition,npart);

			/*Copy onto our dakota responses: */
			if(my_rank==0){
				/*plug response: */
				for(j=0;j<npart;j++)responses_pointer[j]=qmu_response[j];

				/*increment response_pointer :*/
				responses_pointer+=npart;
			}

			/*Free resources:*/
			xDelete<double>(vertex_response);
			xDelete<double>(qmu_response);

		}
		else if (flag==IndexedEnum){

			/*indexed response: plug index into parameters and call response module: */
			parameters->SetParam(index,IndexEnum);

			this->Responsex(&femmodel_response,root);

			if(my_rank==0){
				/*plug response: */
				responses_pointer[0]=femmodel_response;

				/*increment response_pointer :*/
				responses_pointer++;
			}
		}
		else if (flag==NodalEnum){
			_error_("nodal response functions not supported yet!");

			/*increment response_pointer :*/
			responses_pointer++;
		}
		else if (flag==RegularEnum){

			/*perfectly normal response function: */
			this->Responsex(&femmodel_response,root);

			if(my_rank==0){
				/*plug response: */
				responses_pointer[0]=femmodel_response;

				/*increment response_pointer :*/
				responses_pointer++;
			}
		}
		else _error_("flag type " << flag << " not supported yet for response analysis");
	}

	/*Synthesize echo: {{{*/
	if(my_rank==0){
		_printf_("   responses: " << d_numresponses << ": ");
		for(i=0;i<d_numresponses-1;i++)_printf_(d_responses[i] << "|");
		_printf_(d_responses[d_numresponses-1]);
		_printf_("\n");
	}
	/*}}}*/

	/*Free resources:*/
	for(i=0;i<response_partitions_num;i++){
		IssmDouble* matrix=response_partitions[i];
		xDelete<IssmDouble>(matrix);
	}
	xDelete<IssmDouble*>(response_partitions);

}
/*}}}*/
#endif
#ifdef _HAVE_ESA_
void FemModel::EsaGeodetic2D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pGravity, Vector<IssmDouble>* pX, Vector<IssmDouble>* pY, IssmDouble* xx, IssmDouble* yy){/*{{{*/

	int         ns,nsmax;

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the esa geodetic core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->EsaGeodetic2D(pUp,pNorth,pEast,pGravity,pX,pY,xx,yy);
		}
		if(i%100==0){
			pUp->Assemble();
			pNorth->Assemble();
			pEast->Assemble();
			pGravity->Assemble();
			pX->Assemble();
			pY->Assemble();
		}
	}

	/*One last time: */
	pUp->Assemble();
	pNorth->Assemble();
	pEast->Assemble();
	pGravity->Assemble();
	pX->Assemble();
	pY->Assemble();

	/*Free resources:*/
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
}
/*}}}*/
void FemModel::EsaGeodetic3D(Vector<IssmDouble>* pUp, Vector<IssmDouble>* pNorth, Vector<IssmDouble>* pEast, Vector<IssmDouble>* pGravity, IssmDouble* latitude, IssmDouble* longitude, IssmDouble* radius, IssmDouble* xx, IssmDouble* yy, IssmDouble* zz){/*{{{*/

	int         ns,nsmax;

	/*Go through elements, and add contribution from each element to the deflection vector wg:*/
	ns = elements->Size();

	/*Figure out max of ns: */
	ISSM_MPI_Reduce(&ns,&nsmax,1,ISSM_MPI_INT,ISSM_MPI_MAX,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&nsmax,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Call the esa geodetic core: */
	for(int i=0;i<nsmax;i++){
		if(i<ns){
			Element* element=xDynamicCast<Element*>(elements->GetObjectByOffset(i));
			element->EsaGeodetic3D(pUp,pNorth,pEast,pGravity,latitude,longitude,radius,xx,yy,zz);
		}
		if(i%100==0){
			pUp->Assemble();
			pNorth->Assemble();
			pEast->Assemble();
			pGravity->Assemble();
		}
	}

	/*One last time: */
	pUp->Assemble();
	pNorth->Assemble();
	pEast->Assemble();
	pGravity->Assemble();

	/*Free resources:*/
	xDelete<IssmDouble>(latitude);
	xDelete<IssmDouble>(longitude);
	xDelete<IssmDouble>(radius);
	xDelete<IssmDouble>(xx);
	xDelete<IssmDouble>(yy);
	xDelete<IssmDouble>(zz);
}
/*}}}*/
#endif
void FemModel::HydrologyEPLupdateDomainx(IssmDouble* pEplcount){ /*{{{*/

	Vector<IssmDouble> *mask = NULL;
	Vector<IssmDouble> *recurence = NULL;
	Vector<IssmDouble> *active = NULL;
	IssmDouble         *serial_mask       = NULL;
	IssmDouble         *serial_rec        = NULL;
	IssmDouble         *serial_active     = NULL;
	IssmDouble         *old_active        = NULL;
	int                *eplzigzag_counter = NULL;
	int                 eplflip_lock;

	HydrologyDCEfficientAnalysis   *effanalysis  = new HydrologyDCEfficientAnalysis();
	HydrologyDCInefficientAnalysis *inefanalysis = new HydrologyDCInefficientAnalysis();

	/*Step 1: update mask, the mask might be extended by residual and/or using downstream sediment head*/
	int numnodes = this->nodes->NumberOfNodes();
	mask=new Vector<IssmDouble>(numnodes);
	recurence=new Vector<IssmDouble>(numnodes);
	this->parameters->FindParam(&eplzigzag_counter,NULL,EplZigZagCounterEnum);
	this->parameters->FindParam(&eplflip_lock,HydrologydcEplflipLockEnum);
	GetVectoronBaseFromInputsx(&old_active,this,HydrologydcMaskEplactiveNodeEnum,NodeSIdEnum);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		effanalysis->HydrologyEPLGetMask(mask,recurence,element);
	}

	/*check for changes and increment zigzag counter, change the mask if necessary*/
	recurence->Assemble();
	serial_rec=recurence->ToMPISerial();
	for(Object* & object : this->nodes->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(serial_rec[node->Sid()]==1.)eplzigzag_counter[node->Lid()] ++;
		if(eplzigzag_counter[node->Lid()]>eplflip_lock && eplflip_lock!=0){
			mask->SetValue(node->Sid(),old_active[node->Sid()],INS_VAL);
		}
	}
	this->parameters->SetParam(eplzigzag_counter,this->nodes->Size(),EplZigZagCounterEnum);
	/*Assemble and serialize*/
	mask->Assemble();
	serial_mask=mask->ToMPISerial();

	xDelete<int>(eplzigzag_counter);
	xDelete<IssmDouble>(serial_rec);
	xDelete<IssmDouble>(old_active);
	delete mask;
	delete recurence;

	/*Update Mask*/
	InputUpdateFromVectorx(this,serial_mask,HydrologydcMaskEplactiveNodeEnum,NodeSIdEnum);
	xDelete<IssmDouble>(serial_mask);
	inefanalysis->ElementizeEplMask(this);
	/*Step 2: update node activity. If one element is connected to mask=1, all nodes are active*/
	active=new Vector<IssmDouble>(nodes->NumberOfNodes());
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		effanalysis->HydrologyEPLGetActive(active,element);
	}

	/*Assemble and serialize*/
	active->Assemble();
	serial_active=active->ToMPISerial();
	delete active;

	/*Update node activation accordingly*/
	int         counter  = 0; //this is probably not acurate but we are only interested in positivity
	Element*    basalelement=NULL;
	int         domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		switch(domaintype){
			case Domain2DhorizontalEnum:
				basalelement = element;
				break;
			case Domain3DEnum:
				if(!element->IsOnBase()) continue;
				basalelement = element->SpawnBasalElement();
				break;
			default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
		}

		int         numnodes = basalelement->GetNumberOfNodes();
		IssmDouble *base     = xNew<IssmDouble>(numnodes);
		basalelement->GetInputListOnNodes(&base[0],BaseEnum);
		for(int in=0;in<numnodes;in++){
			Node* node=basalelement->GetNode(in);
			if(serial_active[node->Sid()]==1.){
				node->Activate();
				if(!node->IsClone()) counter++;
			}
			else{
				node->Deactivate();
				node->ApplyConstraint(0,base[in]);
			}
		}
		xDelete<IssmDouble>(base);
		if(basalelement->IsSpawnedElement()){basalelement->DeleteMaterials(); delete basalelement;};
	}
	xDelete<IssmDouble>(serial_active);
	delete effanalysis;
	delete inefanalysis;
	int sum_counter;
	ISSM_MPI_Reduce(&counter,&sum_counter,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	counter=sum_counter;
	*pEplcount = counter;
	if(VerboseSolution()) {
		if(counter==0){
			_printf0_("   No nodes are active in EPL layer \n");
		}
		else {
			_printf0_("   Some active nodes in EPL layer \n");
		}
	}

	/*Update dof indexings*/
	this->UpdateConstraintsx();
}
/*}}}*/
void FemModel::HydrologyIDSupdateDomainx(IssmDouble* pIDScount){ /*{{{*/

	bool                isthermal;
	Vector<IssmDouble>* mask				= NULL;
	Vector<IssmDouble>* active				= NULL;
	IssmDouble*         serial_mask	= NULL;
	IssmDouble*         serial_active	= NULL;

	HydrologyDCInefficientAnalysis* inefanalysis =  new HydrologyDCInefficientAnalysis();
	parameters->FindParam(&isthermal,TransientIsthermalEnum);

	/*When solving a thermal model we update the thawed nodes*/
	if(isthermal){
		/*Step 1: update mask, the mask correspond to thawed nodes (that have a meltingrate)*/
		mask=new Vector<IssmDouble>(this->nodes->NumberOfNodes());

		for(Object* & object : this->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			inefanalysis->HydrologyIDSGetMask(mask,element);
		}
		/*Assemble and serialize*/
		mask->Assemble();
		serial_mask=mask->ToMPISerial();
		delete mask;
	}
	/*for other cases we just grab the mask from the initialisation value*/
	else{
		GetVectoronBaseFromInputsx(&serial_mask,this,HydrologydcMaskThawedNodeEnum,NodeSIdEnum);
	}
	/*Update Mask and elementize*/
	InputUpdateFromVectorx(this,serial_mask,HydrologydcMaskThawedNodeEnum,NodeSIdEnum);
	xDelete<IssmDouble>(serial_mask);
	inefanalysis->ElementizeIdsMask(this);

	/*get node mask coherent with element mask*/
	active=new Vector<IssmDouble>(nodes->NumberOfNodes());
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		inefanalysis->HydrologyIdsGetActive(active,element);
	}

	/*Assemble and serialize*/
	active->Assemble();
	serial_active=active->ToMPISerial();
	delete active;

	/*Update node activation accordingly*/
	int         counter  = 0; //this is probably not acurate but we are only interested in positivity
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		int         numnodes = element->GetNumberOfNodes();
		IssmDouble *base     = xNew<IssmDouble>(numnodes);

		element->GetInputListOnNodes(&base[0],BaseEnum);

		for(int in=0;in<numnodes;in++){
			Node* node=element->GetNode(in);
			if(serial_active[node->Sid()]==1.){
				node->Activate();
				if(!node->IsClone()) counter++;
			}
			else{
				node->Deactivate();
				node->ApplyConstraint(0,base[in]);
			}
		}
		xDelete<IssmDouble>(base);
	}
	xDelete<IssmDouble>(serial_active);
	delete inefanalysis;
	int sum_counter;
	ISSM_MPI_Reduce(&counter,&sum_counter,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	counter=sum_counter;
	*pIDScount = counter;
	if(VerboseSolution()) {
		if(counter==0){
			_printf0_("   No nodes are active in IDS layer \n");
		}
		else {
			_printf0_("   Some active nodes in IDS layer \n");
		}
	}
	/*Update dof indexings*/
	this->UpdateConstraintsx();

}
/*}}}*/
void FemModel::UpdateConstraintsL2ProjectionEPLx(IssmDouble* pL2count){ /*{{{*/

	Vector<IssmDouble>* active        = NULL;
	IssmDouble*         serial_active = NULL;
	HydrologyDCEfficientAnalysis* effanalysis = new HydrologyDCEfficientAnalysis();

	/*update node activity. If one element is connected to mask=1, all nodes are active*/
	this->SetCurrentConfiguration(HydrologyDCEfficientAnalysisEnum);
	active=new Vector<IssmDouble>(nodes->NumberOfNodes());
	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
		effanalysis->HydrologyEPLGetActive(active,element);
	}

	/*Assemble and serialize*/
	active->Assemble();
	serial_active=active->ToMPISerial();
	delete active;
	delete effanalysis;

	/*Update node activation accordingly*/
	int counter =0;
	this->SetCurrentConfiguration(L2ProjectionEPLAnalysisEnum);
	for(Object* & object : this->nodes->objects){
		Node* node = xDynamicCast<Node*>(object);
		if(serial_active[node->Sid()]==1.){
			node->Activate();
			if(!node->IsClone()) counter++;
		}
		else{
			node->Deactivate();
		}
	}
	xDelete<IssmDouble>(serial_active);
	int sum_counter;
	ISSM_MPI_Reduce(&counter,&sum_counter,1,ISSM_MPI_INT,ISSM_MPI_SUM,0,IssmComm::GetComm() );
	ISSM_MPI_Bcast(&sum_counter,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	counter=sum_counter;
	*pL2count = counter;
	if(VerboseSolution()) _printf0_("   Number of active nodes L2 Projection: "<< counter <<"\n");
}
/*}}}*/
void FemModel::InitTransientInputx(int* transientinput_enum,int numoutputs){ /*{{{*/

	for(int i=0;i<numoutputs;i++){
		this->inputs->DeleteInput(transientinput_enum[i]);
		this->inputs->SetTransientInput(transientinput_enum[i],NULL,0);
		/*We need to configure this input!*/
		TransientInput* transientinput = this->inputs->GetTransientInput(transientinput_enum[i]); _assert_(transientinput);
		transientinput->Configure(this->parameters);
	}
}
/*}}}*/
void FemModel::StackTransientInputx(int* input_enum,int* transientinput_enum,IssmDouble subtime,int numoutputs){ /*{{{*/

  for(int i=0;i<numoutputs;i++){
		if(input_enum[i]<0){
			_error_("Can't deal with non enum fields for result Stack");
		}
		else{
			for(Object* & object : this->elements->objects){
				/*Get the right transient input*/
				Element* element = xDynamicCast<Element*>(object);
				TransientInput* transientinput = this->inputs->GetTransientInput(transientinput_enum[i]);

				/*Get values and lid list*/
				const int   numvertices = element->GetNumberOfVertices();
				IssmDouble* values=xNew<IssmDouble>(numvertices);
				int        *vertexlids = xNew<int>(numvertices);
				element->GetInputListOnVertices(&values[0],input_enum[i]);   //this is the enum to stack

				element->GetVerticesLidList(vertexlids);

				switch(element->ObjectEnum()){
					case TriaEnum:  transientinput->AddTriaTimeInput(subtime,numvertices,vertexlids,values,P1Enum); break;
					case PentaEnum: transientinput->AddPentaTimeInput(subtime,numvertices,vertexlids,values,P1Enum); break;
					default: _error_("Not implemented yet");
				}
				xDelete<IssmDouble>(values);
				xDelete<int>(vertexlids);
			}
		}
	}
}
/*}}}*/
void FemModel::StackTransientInputonBasex(int* input_enum,int* transientinput_enum,IssmDouble subtime,int numoutputs){ /*{{{*/

	Element*   basalelement=NULL;
	int      domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);

	for(int i=0;i<numoutputs;i++){
		if(input_enum[i]<0){
			_error_("Can't deal with non enum fields for result Stack");
		}
		else{
			for(Object* & object : this->elements->objects){
				/*Get the right transient input*/
				Element* element = xDynamicCast<Element*>(object);
				/*Get basal element*/
				switch(domaintype){
					case Domain2DhorizontalEnum:
						break;
					case Domain3DEnum:
						if(!element->IsOnBase()) continue;
						break;
					default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
				}
				/*Get values and lid list*/
				TransientInput* transientinput = this->inputs->GetTransientInput(transientinput_enum[i]);
				const int   numvertices = element->GetNumberOfVertices();
				IssmDouble* values      = xNew<IssmDouble>(numvertices);
				int        *vertexlids  = xNew<int>(numvertices);
				switch(domaintype){
					case Domain2DhorizontalEnum:
						element->GetInputListOnVertices(&values[0],input_enum[i]); //this is the enum to stack
						element->GetVerticesLidList(vertexlids);
						transientinput->AddTriaTimeInput(subtime,numvertices,vertexlids,values,P1Enum);
						break;
					case Domain3DEnum:{
						element->GetInputListOnVertices(&values[0],input_enum[i]); //this is the enum to stack
						element->GetVerticesLidList(vertexlids);
						Penta* penta=xDynamicCast<Penta*>(object);
						for(;;){
							transientinput->AddPentaTimeInput(subtime,numvertices,vertexlids,values,P1Enum);
							if (penta->IsOnSurface()) break;
							penta=penta->GetUpperPenta(); _assert_(penta->Id()!=element->id);
							penta->GetVerticesLidList(vertexlids);
						}
						break;
					}
					default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
				}
				xDelete<IssmDouble>(values);
				xDelete<int>(vertexlids);
			}
		}
	}
}
/*}}}*/
void FemModel::AverageTransientInputx(int* transientinput_enum,int* averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int numoutputs, int averaging_method){ /*{{{*/

	for(int i=0;i<numoutputs;i++){
		for(Object* & object : this->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			element->CreateInputTimeAverage(transientinput_enum[i],averagedinput_enum[i],init_time,end_time,averaging_method);
		}
	}
}/*}}}*/
void FemModel::AverageTransientInputonBasex(int* transientinput_enum,int* averagedinput_enum,IssmDouble init_time,IssmDouble end_time,int numoutputs, int averaging_method){ /*{{{*/

	int      domaintype;
	this->parameters->FindParam(&domaintype,DomainTypeEnum);

	for(int i=0;i<numoutputs;i++){
		for(Object* & object : this->elements->objects){
			Element* element = xDynamicCast<Element*>(object);
			/*Get basal element*/
			switch(domaintype){
				case Domain2DhorizontalEnum:
					element->CreateInputTimeAverage(transientinput_enum[i],averagedinput_enum[i],init_time,end_time,averaging_method);;
					break;
				case Domain3DEnum:
					if(!element->IsOnBase()) continue;
					element->CreateInputTimeAverage(transientinput_enum[i],averagedinput_enum[i],init_time,end_time,averaging_method);
					break;
				default: _error_("mesh "<<EnumToStringx(domaintype)<<" not supported yet");
			}
		}
	}
}/*}}}*/
#ifdef _HAVE_JAVASCRIPT_
FemModel::FemModel(IssmDouble* buffer, int buffersize, char* toolkits, char* solution, char* modelname,ISSM_MPI_Comm incomm, bool trace){ /*{{{*/
	/*configuration: */
	int  solution_type;
	int  ierr;

	/*First things first, store the communicator, and set it as a global variable: */
	IssmComm::SetComm(incomm);

	/*Start profiler: */
	this->profiler=new Profiler();
	profiler->Start(TOTAL);

	/*From command line arguments, retrieve different filenames needed to create the FemModel: */
	solution_type=StringToEnumx(solution);

	/*Create femmodel from input files: */
	profiler->Start(MPROCESSOR);
	this->InitFromBuffers((char*)buffer,buffersize,toolkits, solution_type,trace,NULL);
	profiler->Stop(MPROCESSOR);

	/*Save communicator in the parameters dataset: */
	this->parameters->AddObject(new GenericParam<ISSM_MPI_Comm>(incomm,FemModelCommEnum));

}
/*}}}*/
void FemModel::CleanUpJs(char** poutput, size_t* psize){/*{{{*/

	/*Intermediary*/
	FILE *output_fid;
	GenericParam<char**>* outputbufferparam=NULL;
	GenericParam<size_t*>* outputbuffersizeparam=NULL;
	char** poutputbuffer;
	size_t* poutputbuffersize;

	/*Before we delete the profiler, report statistics for this run: */
	profiler->Stop(TOTAL);  //final tagging
	_printf0_("\n");
	_printf0_("   "<<setw(40)<<left<<"FemModel initialization elapsed time:"<<profiler->TotalTime(MPROCESSOR) << "\n");
	_printf0_("   "<<setw(40)<<left<<"Core solution elapsed time:"<<profiler->TotalTime(CORE) << "\n");
	_printf0_("\n");
	_printf0_("   Total elapsed time: "
				<<profiler->TotalTimeModHour(TOTAL)<<" hrs "
				<<profiler->TotalTimeModMin(TOTAL)<<" min "
				<<profiler->TotalTimeModSec(TOTAL)<<" sec"
				);
	_printf0_("\n");

	/*Before we close the output file, recover the buffer and size:*/
	outputbufferparam = xDynamicCast<GenericParam<char**>*>(this->parameters->FindParamObject(OutputBufferPointerEnum));
	poutputbuffer=outputbufferparam->GetParameterValue();
	outputbuffersizeparam = xDynamicCast<GenericParam<size_t*>*>(this->parameters->FindParamObject(OutputBufferSizePointerEnum));
	poutputbuffersize=outputbuffersizeparam->GetParameterValue();

	/*Assign output values: */
	*poutput=*poutputbuffer;
	*psize=*poutputbuffersize;
}
/*}}}*/
void FemModel::InitFromBuffers(char* buffer, int buffersize, char* toolkits, int in_solution_type, bool trace, IssmPDouble* X){/*{{{*/

	/*intermediary*/
	FILE       *IOMODEL = NULL;
	FILE       *toolkitsoptionsfid = NULL;
	FILE       *output_fid = NULL;
	int         my_rank;
	size_t      outputsize;
	char       *outputbuffer;
	const char *rootpath = "";   //needed for Dakota runs only, which we won't do here.

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*Open input file descriptor on cpu 0: */
	if(my_rank==0) IOMODEL = fmemopen((void*)buffer, buffersize, "rb");

	/*Open toolkits file descriptor: */
	toolkitsoptionsfid=fmemopen((void*)toolkits, strlen(toolkits)+1, "r");

	/*Now, go create FemModel:*/
	this->InitFromFids((char*)rootpath,IOMODEL,toolkitsoptionsfid,in_solution_type,trace,X);

	/*Close input file and toolkits file descriptors: */
	if(my_rank==0) fclose(IOMODEL);
	fclose(toolkitsoptionsfid);

	/*Open output file once for all and add output file descriptor to parameters*/
	output_fid=open_memstream(&outputbuffer,&outputsize);
	if(output_fid==NULL)_error_("could not initialize output stream");
	this->parameters->SetParam(output_fid,OutputFilePointerEnum);
	this->parameters->AddObject(new GenericParam<char**>(&outputbuffer,OutputBufferPointerEnum));
	this->parameters->AddObject(new GenericParam<size_t*>(&outputsize,OutputBufferSizePointerEnum));

}/*}}}*/
#endif

#if defined(_HAVE_BAMG_) && !defined(_HAVE_AD_)
void FemModel::ReMeshBamg(int* pnewnumberofvertices,int* pnewnumberofelements,IssmDouble** pnewx,IssmDouble** pnewy,IssmDouble** pnewz,int** pnewelementslist){/*{{{*/

	/*Output*/
	IssmDouble *newx			= NULL;
	IssmDouble *newy			= NULL;
	IssmDouble *newz			= NULL;
	IssmDouble *newxylist   = NULL;
	int *newelementslist		= NULL;
	int* newdatalist        = NULL;
	int newnumberofvertices	= -1;
	int newnumberofelements = -1;

	/*Get Rank*/
	int my_rank	= IssmComm::GetRank();

	/*Intermediaries*/
	int numberofvertices 				= this->vertices->NumberOfVertices();
	IssmDouble* vector_serial			= NULL;
	IssmDouble* hmaxvertices_serial	= NULL;
	Vector<IssmDouble> *vector			= NULL;

	/*Get vector to create metric*/
	if(this->amrbamg->fieldenum!=NoneEnum){
		GetVectorFromInputsx(&vector,this,this->amrbamg->fieldenum,VertexSIdEnum);
		vector->Assemble();
		vector_serial = vector->ToMPISerial();
	}

	/*Get hmaxVertices to create metric*/
	if(this->amrbamg->groundingline_distance>0||this->amrbamg->icefront_distance>0||
		this->amrbamg->thicknesserror_threshold>0||this->amrbamg->deviatoricerror_threshold>0){
		/*Initialize hmaxvertices with NAN*/
		hmaxvertices_serial=xNew<IssmDouble>(numberofvertices);
		for(int i=0;i<numberofvertices;i++) hmaxvertices_serial[i]=NAN;
		/*Fill hmaxvertices*/
		if(this->amrbamg->thicknesserror_threshold>0)	this->GethmaxVerticesFromEstimators(hmaxvertices_serial,ThicknessErrorEstimatorEnum);
		if(this->amrbamg->deviatoricerror_threshold>0)	this->GethmaxVerticesFromEstimators(hmaxvertices_serial,DeviatoricStressErrorEstimatorEnum);
		if(this->amrbamg->groundingline_distance>0)		this->GethmaxVerticesFromZeroLevelSetDistance(hmaxvertices_serial,MaskOceanLevelsetEnum);
		if(this->amrbamg->icefront_distance>0)				this->GethmaxVerticesFromZeroLevelSetDistance(hmaxvertices_serial,MaskIceLevelsetEnum);
	}

	if(my_rank==0){
		this->amrbamg->ExecuteRefinementBamg(vector_serial,hmaxvertices_serial,&newdatalist,&newxylist,&newelementslist);
		if(newdatalist[0]<=0 || newdatalist[1]<=0) _error_("Error in the refinement process.");
	}

   /*Send new mesh to others CPU's*/
   if(my_rank) newdatalist=xNew<int>(2);
   ISSM_MPI_Bcast(newdatalist,2,ISSM_MPI_INT,0,IssmComm::GetComm());
   newnumberofvertices=newdatalist[0];
   newnumberofelements=newdatalist[1];
   if(my_rank){
      newxylist      =xNew<IssmDouble>(newnumberofvertices*2);
      newelementslist=xNew<int>(newnumberofelements*this->GetElementsWidth());
   }
   ISSM_MPI_Bcast(newxylist,newnumberofvertices*2,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
   ISSM_MPI_Bcast(newelementslist,newnumberofelements*this->GetElementsWidth(),ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Reorganize the data*/
   newx=xNew<IssmDouble>(newnumberofvertices);
   newy=xNew<IssmDouble>(newnumberofvertices);
   newz=xNewZeroInit<IssmDouble>(newnumberofvertices);
   for(int i=0;i<newnumberofvertices;i++){
      newx[i] = newxylist[2*i];
      newy[i] = newxylist[2*i+1];
   }

   /*Assign output pointers*/
   *pnewnumberofvertices = newnumberofvertices;
   *pnewnumberofelements = newnumberofelements;
   *pnewx = newx;
   *pnewy = newy;
   *pnewz = newz;
   *pnewelementslist = newelementslist;

   /*Cleanup*/
   xDelete<int>(newdatalist);
   xDelete<IssmDouble>(newxylist);
   xDelete<IssmDouble>(vector_serial);
   xDelete<IssmDouble>(hmaxvertices_serial);
   delete vector;
}
/*}}}*/
void FemModel::InitializeAdaptiveRefinementBamg(void){/*{{{*/

	/*Define variables*/
	int numberofvertices      = this->vertices->NumberOfVertices();
	int numberofelements      = this->elements->NumberOfElements();
	IssmDouble* x             = NULL;
	IssmDouble* y             = NULL;
	int* elements             = NULL;
	IssmDouble hmin,hmax,err,gradation;

   /*Get rank*/
	int my_rank = IssmComm::GetRank();

	/*Initialize field as NULL for now*/
	this->amrbamg = NULL;

	/*Get vertices coordinates of the coarse mesh (father mesh)*/
	this->GetMesh(this->vertices,this->elements,&x,&y,&elements);

	/*Create bamg data structures for bamg*/
	this->amrbamg = new AmrBamg();

	/*Get amr parameters*/
	this->parameters->FindParam(&hmin,AmrHminEnum);
	this->parameters->FindParam(&hmax,AmrHmaxEnum);
	this->parameters->FindParam(&err,AmrErrEnum);
	this->parameters->FindParam(&gradation,AmrGradationEnum);
	this->parameters->FindParam(&this->amrbamg->fieldenum,AmrFieldEnum);
	this->parameters->FindParam(&this->amrbamg->keepmetric,AmrKeepMetricEnum);
	this->parameters->FindParam(&this->amrbamg->groundingline_resolution,AmrGroundingLineResolutionEnum);
	this->parameters->FindParam(&this->amrbamg->groundingline_distance,AmrGroundingLineDistanceEnum);
	this->parameters->FindParam(&this->amrbamg->icefront_resolution,AmrIceFrontResolutionEnum);
	this->parameters->FindParam(&this->amrbamg->icefront_distance,AmrIceFrontDistanceEnum);
	this->parameters->FindParam(&this->amrbamg->thicknesserror_resolution,AmrThicknessErrorResolutionEnum);
	this->parameters->FindParam(&this->amrbamg->thicknesserror_threshold,AmrThicknessErrorThresholdEnum);
	this->parameters->FindParam(&this->amrbamg->thicknesserror_groupthreshold,AmrThicknessErrorGroupThresholdEnum);
	this->parameters->FindParam(&this->amrbamg->thicknesserror_maximum,AmrThicknessErrorMaximumEnum);
	this->parameters->FindParam(&this->amrbamg->deviatoricerror_resolution,AmrDeviatoricErrorResolutionEnum);
	this->parameters->FindParam(&this->amrbamg->deviatoricerror_threshold,AmrDeviatoricErrorThresholdEnum);
	this->parameters->FindParam(&this->amrbamg->deviatoricerror_groupthreshold,AmrDeviatoricErrorGroupThresholdEnum);
	this->parameters->FindParam(&this->amrbamg->deviatoricerror_maximum,AmrDeviatoricErrorMaximumEnum);

	/*Set BamgOpts*/
	this->amrbamg->SetBamgOpts(hmin,hmax,err,gradation);

	/*Re-create original mesh and put it in bamg structure (only cpu 0)*/
	this->amrbamg->SetMesh(&elements,&x,&y,&numberofvertices,&numberofelements);
	if(my_rank==0){
		this->amrbamg->Initialize();
	}
}
/*}}}*/
void FemModel::GethmaxVerticesFromZeroLevelSetDistance(IssmDouble* hmaxvertices,int levelset_type){/*{{{*/

	if(!hmaxvertices) _error_("hmaxvertices is NULL!\n");

	/*Intermediaries*/
	int numberofvertices			 				= this->vertices->NumberOfVertices();
   Vector<IssmDouble>* vminvertexdistance = new Vector<IssmDouble>(numberofvertices);
	IssmDouble* pminvertexdistance 			= NULL;
	IssmDouble* levelset_points				= NULL;
	IssmDouble x,y;
	IssmDouble threshold,resolution;
	IssmDouble minvertexdistance,distance;
	int sid,numberofpoints;

	switch(levelset_type){
		case MaskOceanLevelsetEnum:
			threshold	= this->amrbamg->groundingline_distance;
			resolution	= this->amrbamg->groundingline_resolution;
			break;
		case MaskIceLevelsetEnum:
			threshold	= this->amrbamg->icefront_distance;
			resolution	= this->amrbamg->icefront_resolution;
			break;
		default: _error_("not implemented yet");
	}

	/*Get points which level set is zero (center of elements with zero level set)*/
   this->GetZeroLevelSetPoints(&levelset_points,numberofpoints,levelset_type);//levelset_points is serial (global)

	for(Object* & object : this->vertices->objects){
		Vertex* vertex = xDynamicCast<Vertex*>(object);
      /*Attention: no spherical coordinates*/
      x	= vertex->GetX();
      y 	= vertex->GetY();
      sid= vertex->Sid();
		minvertexdistance=INFINITY;

		/*Find the minimum vertex distance*/
		for(int j=0;j<numberofpoints;j++){
         distance=sqrt((x-levelset_points[2*j])*(x-levelset_points[2*j])+(y-levelset_points[2*j+1])*(y-levelset_points[2*j+1]));
         minvertexdistance=min(distance,minvertexdistance);
    	}
		/*Now, insert in the vector*/
		vminvertexdistance->SetValue(sid,minvertexdistance,INS_VAL);
   }
	/*Assemble*/
   vminvertexdistance->Assemble();

   /*Assign the pointer*/
  	pminvertexdistance=vminvertexdistance->ToMPISerial();

	/*Fill hmaxVertices*/
	for(int i=0;i<numberofvertices;i++){
		if(pminvertexdistance[i]<threshold){ //hmaxvertices is serial (global)
			if(xIsNan<IssmDouble>(hmaxvertices[i])) hmaxvertices[i]=resolution;
			else hmaxvertices[i]=min(resolution,hmaxvertices[i]);
		}
	}

	/*Cleanup*/
	xDelete<IssmDouble>(pminvertexdistance);
   xDelete<IssmDouble>(levelset_points);
   delete vminvertexdistance;
}
/*}}}*/
void FemModel::GethmaxVerticesFromEstimators(IssmDouble* hmaxvertices,int errorestimator_type){/*{{{*/

	if(!hmaxvertices) _error_("hmaxvertices is NULL!\n");

	/*Intermediaries*/
	int elementswidth							= this->GetElementsWidth();
	int numberofelements						= -1;
	int numberofvertices						= -1;
	IssmDouble	hmax							= this->amrbamg->GetBamgOpts()->hmax;
	IssmDouble* maxlength					= NULL;
	IssmDouble* error_vertices				= NULL;
	IssmDouble* error_elements				= NULL;
	IssmDouble* x								= NULL;
	IssmDouble* y								= NULL;
	int* index									= NULL;
	IssmDouble maxerror,threshold,groupthreshold,resolution,length;
	IssmDouble L1,L2,L3;
	int vid,v1,v2,v3;
	bool refine;

	/*Fill variables*/
	switch(errorestimator_type){
		case ThicknessErrorEstimatorEnum:
			threshold		= this->amrbamg->thicknesserror_threshold;
			groupthreshold	= this->amrbamg->thicknesserror_groupthreshold;
			resolution		= this->amrbamg->thicknesserror_resolution;
			maxerror			= this->amrbamg->thicknesserror_maximum;
			this->ThicknessZZErrorEstimator(&error_elements);//error is serial, but the calculation is parallel
			break;
		case DeviatoricStressErrorEstimatorEnum:
			threshold		= this->amrbamg->deviatoricerror_threshold;
			groupthreshold	= this->amrbamg->deviatoricerror_groupthreshold;
			resolution		= this->amrbamg->deviatoricerror_resolution;
			maxerror			= this->amrbamg->deviatoricerror_maximum;
			this->ZZErrorEstimator(&error_elements);//error is serial, but the calculation is parallel
			break;
		default: _error_("not implemented yet");
	}
	if(!error_elements) _error_("error_elements is NULL!\n");
	if(groupthreshold<DBL_EPSILON) _error_("group threshold is too small!");

	/*Get mesh*/
	this->GetMesh(&index,&x,&y,&numberofvertices,&numberofelements);
	if(numberofelements<0) _error_("number of elements is negative!\n");
	if(numberofvertices<0) _error_("number of vertices is negative!\n");
	maxlength		= xNew<IssmDouble>(numberofelements);
	error_vertices	= xNewZeroInit<IssmDouble>(numberofvertices);

	/*Find the max of the estimators if it was not provided*/
	if(maxerror<DBL_EPSILON){
		for(int i=0;i<numberofelements;i++) maxerror=max(maxerror,error_elements[i]);
		switch(errorestimator_type){
      	case ThicknessErrorEstimatorEnum:			this->amrbamg->thicknesserror_maximum 	= maxerror;break;
      	case DeviatoricStressErrorEstimatorEnum: 	this->amrbamg->deviatoricerror_maximum = maxerror;break;
   	}
	}

	/*Fill error_vertices (this is the sum of all elements connected to the vertex)*/
	for(int i=0;i<numberofelements;i++){
		v1=index[i*elementswidth+0]-1;//Matlab to C indexing
		v2=index[i*elementswidth+1]-1;//Matlab to C indexing
		v3=index[i*elementswidth+2]-1;//Matlab to C indexing
		L1=sqrt(pow(x[v2]-x[v1],2)+pow(y[v2]-y[v1],2));
		L2=sqrt(pow(x[v3]-x[v2],2)+pow(y[v3]-y[v2],2));
		L3=sqrt(pow(x[v1]-x[v3],2)+pow(y[v1]-y[v3],2));
		/*Fill the vectors*/
		maxlength[i]		=max(L1,max(L2,L3));
		error_vertices[v1]+=error_elements[i];
		error_vertices[v2]+=error_elements[i];
		error_vertices[v3]+=error_elements[i];
	}

	/*Fill hmaxvertices with the criteria*/
	for(int i=0;i<numberofelements;i++){
		/*Refine any element if its error > phi*maxerror*/
		if(error_elements[i]>threshold*maxerror){
			/*Now, fill the hmaxvertices if requested*/
			for(int j=0;j<elementswidth;j++){
				vid=index[i*elementswidth+j]-1;//Matlab to C indexing
				if(xIsNan<IssmDouble>(hmaxvertices[vid])) hmaxvertices[vid]=max(maxlength[i]/2.,resolution);//Try first dividing the element
				else hmaxvertices[vid]=min(max(maxlength[i]/2.,resolution),hmaxvertices[vid]);//Try first dividing the element
			}
		}
		else {
			/*Try unrefine the element*/
			if(maxlength[i] < 1.1*hmax/2.){
				for(int j=0;j<elementswidth;j++){
					vid=index[i*elementswidth+j]-1;//Matlab to C indexing
					if(error_vertices[vid]>groupthreshold*maxerror) hmaxvertices[vid]=maxlength[i]; //keep the current resolution
					else{
						if(xIsNan<IssmDouble>(hmaxvertices[vid])) hmaxvertices[vid]=min(maxlength[i]*2.,hmax);
						else hmaxvertices[vid]=min(min(maxlength[i]*2.,hmax),hmaxvertices[vid]);//Try first to duplicate the element
					}
				}
			}
		}
	}

	/*Cleanup*/
   xDelete<IssmDouble>(error_elements);
   xDelete<IssmDouble>(error_vertices);
   xDelete<IssmDouble>(maxlength);
}
/*}}}*/
void FemModel::GetVerticeDistanceToZeroLevelSet(IssmDouble** pverticedistance,int levelset_type){/*{{{*/

	//itapopo esse metodo pode ser deletado

	/*Here, "zero level set" means grounding line or ice front, depending on the level set type*/
	/*pverticedistance is the minimal vertice distance to the grounding line or ice front*/
	if(levelset_type!=MaskOceanLevelsetEnum && levelset_type!=MaskIceLevelsetEnum){
		_error_("level set type not implemented yet!");
	}

	/*Output*/
	IssmDouble* verticedistance;

	/*Intermediaries*/
   int numberofvertices       = -1;
	int numberofelements			= -1;
   IssmDouble* levelset_points= NULL;
   IssmDouble* x					= NULL;
   IssmDouble* y					= NULL;
	int* elementslist				= NULL;
	int numberofpoints;
	IssmDouble distance;

	/*Get vertices coordinates*/
	this->GetMesh(&elementslist,&x,&y,&numberofvertices,&numberofelements);
	//this->GetMeshOnPartition(this->vertices,this->elements,&x,&y,&z,&elementslist,&sidtoindex);

	/*Get points which level set is zero (center of elements with zero level set)*/
	this->GetZeroLevelSetPoints(&levelset_points,numberofpoints,levelset_type);

	/*Find the minimal vertice distance to the zero levelset (grounding line or ice front)*/
	verticedistance=xNew<IssmDouble>(numberofvertices);
	for(int i=0;i<numberofvertices;i++){
		verticedistance[i]=INFINITY;
		for(int j=0;j<numberofpoints;j++){
			distance=sqrt((x[i]-levelset_points[2*j])*(x[i]-levelset_points[2*j])+(y[i]-levelset_points[2*j+1])*(y[i]-levelset_points[2*j+1]));
			verticedistance[i]=min(distance,verticedistance[i]);
		}
	}

	/*Assign the pointer*/
	(*pverticedistance)=verticedistance;

	/*Cleanup*/
   xDelete<IssmDouble>(levelset_points);
}
/*}}}*/
#endif

#if defined(_HAVE_NEOPZ_) && !defined(_HAVE_AD_)
void FemModel::ReMeshNeopz(int* pnewnumberofvertices,int* pnewnumberofelements,IssmDouble** pnewx,IssmDouble** pnewy,IssmDouble** pnewz,int** pnewelementslist){/*{{{*/

	/*pnewelementslist keep vertices in Matlab indexing*/
   int my_rank						= IssmComm::GetRank();
	IssmDouble* gl_distance		= NULL;
	IssmDouble* if_distance		= NULL;
	IssmDouble* deviatoricerror= NULL;
	IssmDouble* thicknesserror	= NULL;
	IssmDouble* newx				= NULL;
   IssmDouble* newy				= NULL;
   IssmDouble* newz				= NULL;
	IssmDouble* newxylist		= NULL;
   int* newelementslist			= NULL;
	int* newdatalist				= NULL;
	int newnumberofvertices		= -1;
	int newnumberofelements		= -1;

	/*Get fields, if requested*/
	if(this->amr->groundingline_distance>0)		this->GetElementDistanceToZeroLevelSet(&gl_distance,MaskOceanLevelsetEnum);
   if(this->amr->icefront_distance>0)				this->GetElementDistanceToZeroLevelSet(&if_distance,MaskIceLevelsetEnum);
   if(this->amr->thicknesserror_threshold>0)		this->ThicknessZZErrorEstimator(&thicknesserror);
	if(this->amr->deviatoricerror_threshold>0)	this->ZZErrorEstimator(&deviatoricerror);

	if(my_rank==0){
		this->amr->ExecuteRefinement(gl_distance,if_distance,deviatoricerror,thicknesserror,
												&newdatalist,&newxylist,&newelementslist);
		if(newdatalist[0]<=0 || newdatalist[1]<=0) _error_("Error in the ReMeshNeopz.");
	}

	/*Send new mesh to others CPU's*/
	if(my_rank) newdatalist=xNew<int>(2);
   ISSM_MPI_Bcast(newdatalist,2,ISSM_MPI_INT,0,IssmComm::GetComm());
   newnumberofvertices=newdatalist[0];
   newnumberofelements=newdatalist[1];
   if(my_rank){
      newxylist      =xNew<IssmDouble>(newnumberofvertices*2);
      newelementslist=xNew<int>(newnumberofelements*this->GetElementsWidth());
   }
   ISSM_MPI_Bcast(newxylist,newnumberofvertices*2,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
   ISSM_MPI_Bcast(newelementslist,newnumberofelements*this->GetElementsWidth(),ISSM_MPI_INT,0,IssmComm::GetComm());

   /*Reorganize the data*/
   newx=xNew<IssmDouble>(newnumberofvertices);
   newy=xNew<IssmDouble>(newnumberofvertices);
   newz=xNewZeroInit<IssmDouble>(newnumberofvertices);
   for(int i=0;i<newnumberofvertices;i++){
      newx[i] = newxylist[2*i];
      newy[i] = newxylist[2*i+1];
   }

	/*Assign the pointers*/
   (*pnewelementslist)  = newelementslist; //Matlab indexing
   (*pnewx)             = newx;
   (*pnewy)             = newy;
   (*pnewz)             = newz;
   *pnewnumberofvertices= newnumberofvertices;
   *pnewnumberofelements= newnumberofelements;

	/*Cleanup*/
   xDelete<int>(newdatalist);
   xDelete<IssmDouble>(newxylist);
   xDelete<IssmDouble>(deviatoricerror);
   xDelete<IssmDouble>(thicknesserror);
   xDelete<IssmDouble>(gl_distance);
   xDelete<IssmDouble>(if_distance);
}
/*}}}*/
void FemModel::InitializeAdaptiveRefinementNeopz(void){/*{{{*/

	/*Define variables*/
	int my_rank										= IssmComm::GetRank();
	int numberofvertices							= this->vertices->NumberOfVertices();
	int numberofelements							= this->elements->NumberOfElements();
	IssmDouble* x									= NULL;
	IssmDouble* y									= NULL;
	int* elements									= NULL;
	int amr_restart;

	/*Initialize field as NULL for now*/
	this->amr = NULL;

	/*Get vertices coordinates of the coarse mesh (father mesh)*/
	/*elements comes in Matlab indexing*/
	this->GetMesh(this->vertices,this->elements,&x,&y,&elements);

	/*Create initial mesh (coarse mesh) in neopz data structure*/
	/*Just CPU #0 should keep AMR object*/
   /*Initialize refinement pattern*/
	this->SetRefPatterns();
	this->amr = new AmrNeopz();
	this->amr->refinement_type=1;//1 is refpattern; 0 is uniform (faster)
	/*Get amr parameters*/
	this->parameters->FindParam(&this->amr->level_max,AmrLevelMaxEnum);
	this->parameters->FindParam(&this->amr->gradation,AmrGradationEnum);
	this->parameters->FindParam(&this->amr->lag,AmrLagEnum);
	this->parameters->FindParam(&this->amr->groundingline_distance,AmrGroundingLineDistanceEnum);
	this->parameters->FindParam(&this->amr->icefront_distance,AmrIceFrontDistanceEnum);
	this->parameters->FindParam(&this->amr->thicknesserror_threshold,AmrThicknessErrorThresholdEnum);
	this->parameters->FindParam(&this->amr->thicknesserror_groupthreshold,AmrThicknessErrorGroupThresholdEnum);
	this->parameters->FindParam(&this->amr->thicknesserror_maximum,AmrThicknessErrorMaximumEnum);
	this->parameters->FindParam(&this->amr->deviatoricerror_threshold,AmrDeviatoricErrorThresholdEnum);
	this->parameters->FindParam(&this->amr->deviatoricerror_groupthreshold,AmrDeviatoricErrorGroupThresholdEnum);
	this->parameters->FindParam(&this->amr->deviatoricerror_maximum,AmrDeviatoricErrorMaximumEnum);

	/*Initialize NeoPZ data structure*/
	this->amr->SetMesh(&elements,&x,&y,&numberofvertices,&numberofelements);
	if(my_rank==0){
		this->parameters->FindParam(&amr_restart,AmrRestartEnum);
		if(amr_restart){//experimental
			this->amr->ReadMesh();
		} else {//this is the default method
			this->amr->Initialize();
		}
	}
}
/*}}}*/
void FemModel::GetElementDistanceToZeroLevelSet(IssmDouble** pelementdistance,int levelset_type){/*{{{*/

	/*Here, "zero level set" means grounding line or ice front, depending on the level set type*/
	/*pverticedistance is the minimal vertice distance to the grounding line or ice front*/
	if(levelset_type!=MaskOceanLevelsetEnum && levelset_type!=MaskIceLevelsetEnum){
		_error_("level set type not implemented yet!");
	}

	/*Output*/
	IssmDouble* elementdistance;

	/*Intermediaries*/
   int numberofelements							= this->elements->NumberOfElements();
   Vector<IssmDouble>* velementdistance	= new Vector<IssmDouble>(numberofelements);
   IssmDouble* levelset_points				= NULL;
	IssmDouble* xyz_list							= NULL;
	IssmDouble mindistance,distance;
	IssmDouble xc,yc,x1,y1,x2,y2,x3,y3;
	int numberofpoints;

	/*Get points which level set is zero (center of elements with zero level set, levelset_points is serial)*/
	this->GetZeroLevelSetPoints(&levelset_points,numberofpoints,levelset_type);

	for(Object* & object : this->elements->objects){
		Element* element = xDynamicCast<Element*>(object);
      int sid = element->Sid();
		element->GetVerticesCoordinates(&xyz_list);
		x1 = xyz_list[3*0+0];y1 = xyz_list[3*0+1];
		x2 = xyz_list[3*1+0];y2 = xyz_list[3*1+1];
		x3 = xyz_list[3*2+0];y3 = xyz_list[3*2+1];
		xc = (x1+x2+x3)/3.;
		yc = (y1+y2+y3)/3.;
		mindistance=INFINITY;
		/*Loop over each point (where level set is zero)*/
		for(int j=0;j<numberofpoints;j++){
			distance =sqrt((xc-levelset_points[2*j])*(xc-levelset_points[2*j])+(yc-levelset_points[2*j+1])*(yc-levelset_points[2*j+1]));
			mindistance=min(distance,mindistance);
		}
		velementdistance->SetValue(sid,mindistance,INS_VAL);
		xDelete<IssmDouble>(xyz_list);
	}

   /*Assemble*/
   velementdistance->Assemble();

	/*Assign the pointer*/
	(*pelementdistance)=velementdistance->ToMPISerial();

	/*Cleanup*/
   xDelete<IssmDouble>(levelset_points);
   xDelete<IssmDouble>(xyz_list);
	delete velementdistance;
}
/*}}}*/
void FemModel::SetRefPatterns(){/*{{{*/

   /*Initialize the global variable of refinement patterns*/
   gRefDBase.InitializeUniformRefPattern(ETriangle);

	/*Insert specifics patterns to ISSM core*/
   std::string filepath  = REFPATTERNDIR;
   std::string filename1 = filepath + "/2D_Triang_Rib_3.rpt";
   std::string filename2 = filepath + "/2D_Triang_Rib_4.rpt";
   std::string filename3 = filepath + "/2D_Triang_Rib_5.rpt";
   std::string filename4 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4.rpt";
   std::string filename5 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_4_permuted.rpt";
   std::string filename6 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5.rpt";
   std::string filename7 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_3_5_permuted.rpt";
   std::string filename8 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_4_5.rpt";
   std::string filename9 = filepath + "/2D_Triang_Rib_OnlyTriang_Side_4_5_permuted.rpt";

   TPZAutoPointer<TPZRefPattern> refpat1 = new TPZRefPattern(filename1);
   TPZAutoPointer<TPZRefPattern> refpat2 = new TPZRefPattern(filename2);
   TPZAutoPointer<TPZRefPattern> refpat3 = new TPZRefPattern(filename3);
   TPZAutoPointer<TPZRefPattern> refpat4 = new TPZRefPattern(filename4);
   TPZAutoPointer<TPZRefPattern> refpat5 = new TPZRefPattern(filename5);
   TPZAutoPointer<TPZRefPattern> refpat6 = new TPZRefPattern(filename6);
   TPZAutoPointer<TPZRefPattern> refpat7 = new TPZRefPattern(filename7);
   TPZAutoPointer<TPZRefPattern> refpat8 = new TPZRefPattern(filename8);
   TPZAutoPointer<TPZRefPattern> refpat9 = new TPZRefPattern(filename9);

   if(!gRefDBase.FindRefPattern(refpat1)) gRefDBase.InsertRefPattern(refpat1);
   if(!gRefDBase.FindRefPattern(refpat2)) gRefDBase.InsertRefPattern(refpat2);
   if(!gRefDBase.FindRefPattern(refpat3)) gRefDBase.InsertRefPattern(refpat3);
   if(!gRefDBase.FindRefPattern(refpat4)) gRefDBase.InsertRefPattern(refpat4);
   if(!gRefDBase.FindRefPattern(refpat5)) gRefDBase.InsertRefPattern(refpat5);
   if(!gRefDBase.FindRefPattern(refpat6)) gRefDBase.InsertRefPattern(refpat6);
   if(!gRefDBase.FindRefPattern(refpat7)) gRefDBase.InsertRefPattern(refpat7);
   if(!gRefDBase.FindRefPattern(refpat8)) gRefDBase.InsertRefPattern(refpat8);
   if(!gRefDBase.FindRefPattern(refpat9)) gRefDBase.InsertRefPattern(refpat9);
}
/*}}}*/
#endif
