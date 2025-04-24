/*! \file IoModel.cpp
 * \brief  file containing the methods that will help in processing the input data coming
 * into ISSM, from Matlab, or through a binary file opened for reading.
 */

/*CODES:
 * 1: boolean constant
 * 2: integer constant
 * 3: IssmDouble constant
 * 5: boolean vector
 * 6: int vector
 * 7: IssmDouble vector*/

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdint.h>

#include "./classes.h"
#include "../shared/io/io.h"
#include "../shared/shared.h"
#include "../classes/Inputs/TransientInput.h"
#include "../toolkits/codipack/CoDiPackGlobal.h"

/*IoConstant class and methods*/
IoConstant::IoConstant(){/*{{{*/
	this->isindependent = false;
	this->name          = NULL;
	this->constant      = NULL;
}
/*}}}*/
IoConstant::~IoConstant(){/*{{{*/
	xDelete<char>(this->name);
	delete this->constant;
}
/*}}}*/
IoConstant::IoConstant(bool value,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->constant      = new BoolParam(IoConstantEnum,value);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/
IoConstant::IoConstant(int value,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->constant      = new IntParam(IoConstantEnum, value);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/
IoConstant::IoConstant(IssmDouble value,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->constant      = new DoubleParam(IoConstantEnum, value);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/
IoConstant::IoConstant(char* value,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->constant      = new StringParam(IoConstantEnum, value);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/
IoConstant::IoConstant(char** value,int numstrings,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->constant      = new StringArrayParam(IoConstantEnum, value, numstrings);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/

/*IoData class and methods*/
IoData::IoData(){/*{{{*/
	this->isindependent = false;
	this->name          = NULL;
	this->code          = -1;
	this->layout        = -1;
	this->M             = 0;
	this->N             = 0;
	this->data          = NULL;
}
/*}}}*/
IoData::~IoData(){/*{{{*/
	xDelete<char>(this->name);
	xDelete<IssmDouble>(this->data);
}
/*}}}*/
IoData::IoData(IssmDouble* matrix,int code_in,int layout_in,int M_in,int N_in,const char* name_in){/*{{{*/
	this->isindependent = false;
	this->code          = code_in;
	this->layout        = layout_in;
	this->M             = M_in;
	this->N             = N_in;
	this->data          = matrix; /*do not copy*/
	_assert_(code_in==5 ||  code_in==6 || code_in==7);

	_assert_(name_in);
	int len=strlen(name_in);
	this->name=xNew<char>(len+1);
	memcpy(this->name,name_in,(len+1)*sizeof(char));
}
/*}}}*/

/*IoModel constructors/destructors*/
IoModel::IoModel(){/*{{{*/

	this->fid=NULL;
	this->solution_enum=-1;

	this->my_elements=NULL;
	this->my_faces=NULL;
	this->my_vfaces=NULL;
	this->my_edges=NULL;
	this->my_vedges=NULL;
	this->my_hedges=NULL;
	this->my_vertices=NULL;
	this->my_vertices_lids=NULL;
	this->epart=NULL;

	this->domaintype=-1;
	this->domaindim=-1;
	this->meshelementtype=-1;
	this->numberofvertices=-1;
	this->numberofelements=-1;
	this->numberoffaces=-1;
	this->numberofverticalfaces=-1;
	this->numberofedges=-1;
	this->numberofverticaledges=-1;
	this->numberofhorizontaledges=-1;
	this->facescols=-1;
	this->elements=NULL;
	this->faces=NULL;
	this->verticalfaces=NULL;
	this->edges=NULL;
	this->verticaledges=NULL;
	this->horizontaledges=NULL;
	this->elementtofaceconnectivity           = NULL;
	this->elementtoverticalfaceconnectivity   = NULL;
	this->elementtoedgeconnectivity           = NULL;
	this->elementtoverticaledgeconnectivity   = NULL;
	this->elementtohorizontaledgeconnectivity = NULL;
	this->singlenodetoelementconnectivity     = NULL;
	this->numbernodetoelementconnectivity     = NULL;
}/*}}}*/
IoModel::IoModel(FILE* iomodel_handle,int solution_enum_in,bool trace,IssmPDouble* X){/*{{{*/

	bool autodiff=false;
	bool iscontrol=false;

	/*First, keep track of the file handle: */
	this->fid=iomodel_handle;

	/*Check that Enums are Synchronized*/
	this->CheckFile();

	/*Keep track of solution*/
	this->solution_enum = solution_enum_in;

	/*If we are running in AD mode, we need to start the trace and declare our independent variables now,
	 *and prevent them from being erased during successive calls to iomodel->FetchConstants, iomodel->FetchData and
	 iomodel->DeleteData:*/
	this->StartTrace(trace);
	this->DeclareIndependents(trace,X);

	/*Initialize and read constants:*/
	this->FetchConstants(); /*this routine goes through the input file, and fetches bool, int, IssmDouble and string only, nothing memory intensive*/

	/*Is this an autodiff run?*/
	this->FindConstant(&autodiff,"md.autodiff.isautodiff");
	this->FindConstant(&iscontrol,"md.inversion.iscontrol");
	if(trace){
		autodiff=true;
	}
	else{
		if(autodiff && !iscontrol)
		 autodiff=true;
		else
		 autodiff=false;
	}
	this->AddConstant(new IoConstant(autodiff,"md.autodiff.isautodiff"));

	/*Initialize permanent data: */
	this->my_elements      = NULL;
	this->my_faces         = NULL;
	this->my_vfaces        = NULL;
	this->my_edges         = NULL;
	this->my_vedges        = NULL;
	this->my_hedges        = NULL;
	this->my_vertices      = NULL;
	this->my_vertices_lids = NULL;
	this->epart            = NULL;

	FindConstant(&this->domaintype,"md.mesh.domain_type");
	FindConstant(&this->meshelementtype,"md.mesh.elementtype");

	FetchData(&this->domaindim,"md.mesh.domain_dimension");
	FetchData(&this->numberofvertices,"md.mesh.numberofvertices");
	FetchData(&this->numberofelements,"md.mesh.numberofelements");
	FetchData(&this->elements,NULL,NULL,"md.mesh.elements");
	this->facescols                           = -1;
	this->faces                               = NULL;
	this->verticalfaces                       = NULL;
	this->edges                               = NULL;
	this->verticaledges                       = NULL;
	this->horizontaledges                     = NULL;
	this->elementtofaceconnectivity           = NULL;
	this->elementtoverticalfaceconnectivity   = NULL;
	this->elementtoedgeconnectivity           = NULL;
	this->elementtoverticaledgeconnectivity   = NULL;
	this->elementtohorizontaledgeconnectivity = NULL;
	this->singlenodetoelementconnectivity     = NULL;
	this->numbernodetoelementconnectivity     = NULL;
}/*}}}*/
IoModel::~IoModel(){/*{{{*/

	/*Delete constants*/
	vector<IoConstant*>::iterator iter1;
	for(iter1=constants.begin();iter1<constants.end();iter1++){
		delete *iter1;
	}
	this->constants.clear();

	/*Delete data*/
	vector<IoData*>::iterator iter2;
	for(iter2=data.begin();iter2<data.end();iter2++){
		#if defined(_ISSM_DEBUG_)
		if(!(*iter2)->isindependent){
			_printf0_("WARNING: IoData \"" << (*iter2)->name << "\" has not been freed (DeleteData has not been called)\n");
		}
		#endif
		delete *iter2;
	}
	this->data.clear();

	xDelete<bool>(this->my_elements);
	xDelete<bool>(this->my_faces);
	xDelete<bool>(this->my_vfaces);
	xDelete<bool>(this->my_edges);
	xDelete<bool>(this->my_vedges);
	xDelete<bool>(this->my_hedges);
	xDelete<bool>(this->my_vertices);
	xDelete<int>(this->my_vertices_lids);
	xDelete<int>(this->epart);

	xDelete<int>(this->elements);
	xDelete<int>(this->faces);
	xDelete<int>(this->verticalfaces);
	xDelete<int>(this->edges);
	xDelete<int>(this->verticaledges);
	xDelete<int>(this->horizontaledges);
	xDelete<int>(this->elementtofaceconnectivity);
	xDelete<int>(this->elementtoverticalfaceconnectivity);
	xDelete<int>(this->elementtoedgeconnectivity);
	xDelete<int>(this->elementtoverticaledgeconnectivity);
	xDelete<int>(this->elementtohorizontaledgeconnectivity);
	xDelete<int>(this->singlenodetoelementconnectivity);
	xDelete<int>(this->numbernodetoelementconnectivity);
}
/*}}}*/

/*IoModel methods*/
void  IoModel::AddConstant(IoConstant* in_constant){/*{{{*/

	_assert_(in_constant);

	/*Go through dataset of constant and check whether it already exists */
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		if(strcmp((*iter)->name,in_constant->name)==0){
			delete in_constant;
			return;
		}
	}

	this->constants.push_back(in_constant);
}
/*}}}*/
void  IoModel::AddConstantIndependent(IoConstant* in_constant){/*{{{*/

	_assert_(in_constant);

	/*Set constant as independent*/
	in_constant->isindependent = true;

	/*Add to constnats*/
	this->AddConstant(in_constant);
}
/*}}}*/
void  IoModel::AddData(IoData* in_data){/*{{{*/

	_assert_(in_data);

	/*Go through dataset of data and check whether it already exists */
	vector<IoData*>::iterator iter;

	for(iter=data.begin();iter<data.end();iter++){
		if(strcmp((*iter)->name,in_data->name)==0){
			delete in_data;
			return;
		}
	}

	this->data.push_back(in_data);
}
/*}}}*/
void  IoModel::AddDataIndependent(IoData* in_data){/*{{{*/

	_assert_(in_data);

	/*Set data as independent*/
	in_data->isindependent = true;

	/*Add to constnats*/
	this->AddData(in_data);
}
/*}}}*/
void  IoModel::CheckFile(void){/*{{{*/

	bool        found;
	int         record_enum,record_name_size;
   long long   record_length;
	char       *record_name = NULL;
	const char *mddot = "md.";

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Check that some fields have been allocated*/
	_assert_(this->fid || my_rank);

	/*Go find in the binary file, the position of the data we want to fetch: */
	if(my_rank==0){ //cpu 0

		/*First set FILE* position to the beginning of the file: */
		fseek(this->fid,0,SEEK_SET);

		for(;;){
			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				xDelete<char>(record_name);
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				_error_("error while looking in binary file. Found a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				/*we have reached the end of the file. break: */
				found=false;
				xDelete<char>(record_name);
				break;
			}
			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}

			/*Have we found the last string?*/
			if(strncmp(record_name,"md.EOF",6)==0){
				found = true;
				xDelete<char>(record_name);
				break;
			}

			/*Go to next Enum*/
			if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
			fseek(fid,record_length,SEEK_CUR);
			xDelete<char>(record_name);
		}
		if(!found){
			_printf0_("\n");
			_printf0_("=========================================================================\n");
			_printf0_(" Marshalled file is corrupted                                            \n");
			_printf0_("                                                                         \n");
			_printf0_("    Last record found is:                                                \n");
			_printf0_("    the corresponding model field has probably been marshalled           \n");
			_printf0_("    incorrectly                                                          \n");
			_printf0_("                                                                         \n");
			_printf0_("=========================================================================\n");
			_printf0_("\n");
			_error_("Binary file corrupted (See error message above)");
		}
	}
}
/*}}}*/
Param* IoModel::CopyConstantObject(const char* constant_name,int param_enum){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			Param* output = ioconstant->constant->copy();
			output->SetEnum(param_enum);
			return output;
		}
	}

	_error_("Constant \"" << constant_name << "\" not found in iomodel");
	return NULL;
}
/*}}}*/
IssmDouble* IoModel::Data(const char* data_name){/*{{{*/

	/*Intermediary*/
	vector<IoData*>::iterator iter;

	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,data_name)==0) return iodata->data;
	}

	return NULL;
}
/*}}}*/
void  IoModel::ConstantToInput(Inputs* inputs,Elements* elements,IssmDouble value, int vector_enum,int type){/*{{{*/

	if (type==P1Enum){
		for(Object* & object : elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			element->InputCreateP1FromConstant(inputs,this,value,vector_enum);
		}
	}
	else if (type==P0Enum){
		for(Object* & object : elements->objects){
			Element* element=xDynamicCast<Element*>(object);
			element->InputCreateP0FromConstant(inputs,this,value,vector_enum);
		}
	}
	else _error_("not supported yet!");
	return;
}
/*}}}*/
void  IoModel::DeclareIndependents(bool trace,IssmPDouble* X){/*{{{*/

	bool autodiff,iscontrol;
	int  num_independent_objects,temp;
	int  Xcount=0;

	char** names = NULL;
	int *types = NULL;

	/*Initialize array detecting whether data[i] is an independent AD mode variable: */
	this->FetchData(&autodiff,"md.autodiff.isautodiff");
	this->FetchData(&iscontrol,"md.inversion.iscontrol");

	if(trace || (autodiff && !iscontrol)){

		#ifdef _HAVE_AD_
		// FIXME codi here we should be able to execute codi version as normal
		this->FetchData(&num_independent_objects,"md.autodiff.num_independent_objects");
		if(num_independent_objects){
			this->FetchMultipleData(&names,&temp,"md.autodiff.independent_name"); _assert_(temp==num_independent_objects);

			/*create independent objects, and at the same time, fetch the corresponding independent variables,
			 *and declare them as such in ADOLC: */
			for(int i=0;i<num_independent_objects;i++){
				this->FetchIndependentData(&Xcount,X,names[i]);
			}
			for(int i=0;i<num_independent_objects;i++) xDelete<char>(names[i]);
			xDelete<char*>(names);
		}
		#else
		/*if we asked for AD computations, we have a problem!: */
		_error_("Cannot carry out AD mode computations without support of ADOLC or CoDiPack compiled in!");
		#endif
	}
}
/*}}}*/
void  IoModel::DeleteData(int num,...){/*{{{*/

	/*Intermediaries*/
	va_list     ap;
	char       *data_name = NULL;
	const char *mddot     = "md.";
	vector<IoData *>::iterator iter;

	/*Go through the entire list of data and delete the corresponding data from the iomodel-data dataset: */
	va_start(ap,num);
	for(int i=0;i<num;i++){
		data_name=va_arg(ap,char*);

		if(strncmp(data_name,mddot,3)!=0) _error_("String provided does not start with \"md.\" ("<<data_name<<")");

		for(iter=data.begin();iter<data.end();iter++){
			IoData* iodata=*iter;
			if(strcmp(iodata->name,data_name)==0 && !iodata->isindependent){
				delete *iter;
				this->data.erase(iter);
				break;
			}
		}
	}
	va_end(ap);
} /*}}}*/
void  IoModel::DeleteData(IssmDouble* vector_in,const char* data_name){/*{{{*/

	vector<IoData*>::iterator iter;

	/*do not do anything if pointer is NULL*/
	if(!vector_in) return;

	/*do not delete if this is an independent variable*/
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,data_name)==0 && iodata->isindependent){
			return;
		}
	}

	/*Go ahead and delete*/
	xDelete<IssmDouble>(vector_in);
} /*}}}*/
void  IoModel::DeleteData(char*** pstringarray, int numstrings,const char* data_name){/*{{{*/

	char** stringarray=*pstringarray;

	if(numstrings){
		for(int i=0;i<numstrings;i++){
			char* string=stringarray[i];
			xDelete<char>(string);
		}
		xDelete<char*>(stringarray);
	}
	*pstringarray=NULL;
} /*}}}*/
void  IoModel::FetchConstants(void){/*{{{*/

	/*record descriptions; */
	const char* mddot = "md.";
	char* record_name = NULL;
	int   record_name_size;
	long long record_length;
	int record_code; //1 to 7 number

	/*records: */
	int          booleanint  = 0;
	int          integer     = 0;
	IssmPDouble  pscalar     = 0;
	IssmDouble   scalar      = 0;
	char        *string      = NULL;
	char       **strings     = NULL;
	int          string_size,numstrings;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Check that some fields have been allocated*/
	_assert_(this->fid || my_rank);

	/*Go find in the binary file, the position of the data we want to fetch: */
	if(my_rank==0){ //cpu 0{{{

		/*First set FILE* position to the beginning of the file: */
		fseek(this->fid,0,SEEK_SET);

		/*Now march through file looking for the correct data identifiers (bool,int,IssmDouble or string): */
		for(;;){

			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				record_code=0; //0 means bailout
				ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*tell others cpus we are bailing: */
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				/*we are going to error out, still try and get informatoin for the user:*/
				record_name=xNew<char>(record_name_size+1);
				record_name[record_name_size]='\0';

				/*Read record_name: */
				if(fread(record_name,record_name_size*sizeof(char),1,fid)!=0){};

				_error_("error while looking in binary file. String \"" << record_name << "\" is a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				_error_("Could not read record name");
			}

			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}

			/* Read the record length and the data type code: */
			if(fread(&record_length,sizeof(long long),1,this->fid)!=1) _error_("Cound not read record_length");
			if(fread(&record_code  ,sizeof(int),1,this->fid)!=1) _error_("Cound not read record_code");

			/*Tell other cpus what we are doing: */
			ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*tell other cpus what we are going to do: */

			/*Tell other cpus the name of the data, then branch according to the data type: */
			ISSM_MPI_Bcast(&record_name_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
			ISSM_MPI_Bcast(record_name,record_name_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
			ISSM_MPI_Bcast(&record_length,1,ISSM_MPI_LONG_LONG_INT,0,IssmComm::GetComm());

			switch(record_code){
				case 1:
					/*Read the boolean and broadcast it to other cpus:*/
					if(fread(&booleanint,sizeof(int),1,this->fid)!=1) _error_("could not read boolean ");
					ISSM_MPI_Bcast(&booleanint,1,ISSM_MPI_INT,0,IssmComm::GetComm());

					/*create BoolParam: */
					this->AddConstant(new IoConstant((bool)booleanint,record_name)); //cast to boolean

					break;
				case 2:
					/*Read the integer and broadcast it to other cpus:*/
					if(fread(&integer,sizeof(int),1,this->fid)!=1) _error_("could not read integer ");

					/*Convert codes to Enums if needed*/
					if(strcmp(record_name,"md.smb.model")==0) integer = IoCodeToEnumSMB(integer);
					if(strcmp(record_name,"md.basalforcings.model")==0) integer = IoCodeToEnumBasal(integer);
					if(strcmp(record_name,"md.calving.law")==0) integer = IoCodeToEnumCalving(integer);
					if(strcmp(record_name,"md.frontalforcings.parameterization")==0) integer = IoCodeToEnumFrontalforcings(integer);
					if(strcmp(record_name,"md.hydrology.model")==0) integer = IoCodeToEnumHydrology(integer);
					if(strcmp(record_name,"md.materials.type")==0) integer = IoCodeToEnumMaterials(integer);
					if(strcmp(record_name,"md.materials.nature")==0) integer = IoCodeToEnumNature(integer);
					if(strcmp(record_name,"md.timestepping.type")==0) integer = IoCodeToEnumTimestepping(integer);
					if(strcmp(record_name,"md.amr.type")==0) integer = IoCodeToEnumAmr(integer);
					if(strcmp(record_name,"md.solidearth.settings.grdmodel")==0) integer = IoCodeToEnumGrd(integer);

					/*Broadcast to other cpus*/
					ISSM_MPI_Bcast(&integer,1,ISSM_MPI_INT,0,IssmComm::GetComm());

					/*create IntParam: */
					this->AddConstant(new IoConstant(integer,record_name));

					break;
				case 3:
					  {
						/*IssmDouble, check whether it is already there (from "declare independents")*/
						bool exists = false;
						vector<IoConstant*>::iterator iter;
						for(iter=constants.begin();iter<constants.end();iter++){
							IoConstant* ioconstant=*iter;
							if(strcmp(ioconstant->name,record_name)==0){
								exists = true;
								break;
							}
						}
						if(!exists){
							if(fread(&pscalar,sizeof(IssmPDouble),1,this->fid)!=1) _error_("could not read scalar ");
							ISSM_MPI_Bcast(&pscalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
							scalar=pscalar;

							/*create DoubleParam: */
							this->AddConstant(new IoConstant(scalar,record_name));
						}
						else{
							if(fread(&pscalar,sizeof(IssmPDouble),1,this->fid)!=1) _error_("could not read scalar ");
						}
					  }
					break;
				case 4:
					/*We have to read a string from disk. First read the dimensions of the string, then the string: */
					if(fread(&string_size,sizeof(int),1,this->fid)!=1) _error_("could not read length of string ");
					ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

					if(string_size){
						string=xNew<char>(string_size+1);
						string[string_size]='\0';

						/*Read string, then broadcast: */
						if(fread(string,string_size*sizeof(char),1,this->fid)!=1)_error_(" could not read string ");
						ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
					}
					else{
						string=xNew<char>(1);
						string[0]='\0';
					}
					/*Convert strings to enums if needed*/
					if(strcmp(record_name,"md.flowequation.fe_SSA")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.flowequation.fe_HO")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.flowequation.fe_FS")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.thermal.fe")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.levelset.fe")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.groundingline.migration")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.groundingline.friction_interpolation")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.groundingline.melt_interpolation")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.masstransport.hydrostatic_adjustment")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.materials.rheology_law")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.damage.elementinterp")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.mesh.domain_type")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else if(strcmp(record_name,"md.mesh.elementtype")==0){
						this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
					} else {
						/*Add string to parameters: */
						this->AddConstant(new IoConstant(string,record_name));
					}

					/*Free string*/
					xDelete<char>(string);
					break;
				case 5:
				case 6:
				case 7:
				case 8:
				case 10:
					/*We are not interested in this record, too memory intensive. Skip it: */
					/*skip: */
					fseek(fid,-sizeof(int),SEEK_CUR); //backtrak 1 integer
					fseek(fid,record_length,SEEK_CUR);
					break;
				case 9:
					/*String Array*/
					if(fread(&numstrings,sizeof(int),1,fid)!=1) _error_("could not read length of string array");
					ISSM_MPI_Bcast(&numstrings,1,ISSM_MPI_INT,0,IssmComm::GetComm());
					/*Now allocate string array: */
					if(numstrings){
						strings=xNew<char*>(numstrings);
						for(int i=0;i<numstrings;i++)strings[i]=NULL;

						/*Go through strings, and read: */
						for(int i=0;i<numstrings;i++){

							if(fread(&string_size,sizeof(int),1,fid)!=1) _error_("could not read length of string ");
							ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
							if(string_size){
								string=xNew<char>((string_size+1));
								string[string_size]='\0';
								if(fread(string,string_size*sizeof(char),1,fid)!=1)_error_(" could not read string ");
								ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
							}
							else{
								string=xNew<char>(1);
								string[0]='\0';
							}
							strings[i]=string;
						}
					}

					/*Add strings to parameters: */
					this->AddConstant(new IoConstant(strings,numstrings,record_name));

					/*Free string*/
					for(int i=0;i<numstrings;i++) xDelete<char>(strings[i]);
					xDelete<char*>(strings);
					break;
				default:
					_error_("unknown record type:" << record_code);
					break;
			}
			xDelete<char>(record_name);
		}
	} //}}}
	else{ //cpu ~0 {{{
		for(;;){ //wait on cpu 0
			ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*get from cpu 0 what we are going to do: */
			if(record_code==0){
				break; //we are done, break from the loop
			}
			else{
				ISSM_MPI_Bcast(&record_name_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
				_assert_(record_name_size);
				record_name=xNew<char>((record_name_size+1)); record_name[record_name_size]='\0';
				ISSM_MPI_Bcast(record_name,record_name_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
				ISSM_MPI_Bcast(&record_length,1,ISSM_MPI_LONG_LONG_INT,0,IssmComm::GetComm());
				switch(record_code){
					case 1:
						/*boolean. get it from cpu 0 */
						ISSM_MPI_Bcast(&booleanint,1,ISSM_MPI_INT,0,IssmComm::GetComm());

						/*create BoolParam: */
						this->AddConstant(new IoConstant((bool)booleanint,record_name)); //cast to a boolean
						break;

					case 2:
						/*integer. get it from cpu 0 */
						ISSM_MPI_Bcast(&integer,1,ISSM_MPI_INT,0,IssmComm::GetComm());

						/*create IntParam: */
						this->AddConstant(new IoConstant(integer,record_name));
						break;
					case 3:
						/*scalar. get it from cpu 0 */
						  {
							/*IssmDouble, check whether it is already there (from "declare independents")*/
							bool exists = false;
							vector<IoConstant*>::iterator iter;
							for(iter=constants.begin();iter<constants.end();iter++){
								IoConstant* ioconstant=*iter;
								if(strcmp(ioconstant->name,record_name)==0){
									exists = true;
									break;
								}
							}
							if(!exists){
								ISSM_MPI_Bcast(&pscalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
								scalar=pscalar;
								/*create DoubleParam: */
								this->AddConstant(new IoConstant(scalar,record_name));
							}
						  }
						break;
					case 4:
						ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
						if(string_size){
							string=xNew<char>((string_size+1));
							string[string_size]='\0';

							/*Read string from cpu 0: */
							ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
						}
						else{
							string=xNew<char>(1);
							string[0]='\0';
						}

						if(strcmp(record_name,"md.flowequation.fe_SSA")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.flowequation.fe_HO")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.flowequation.fe_FS")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.thermal.fe")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.levelset.fe")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.groundingline.migration")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.groundingline.friction_interpolation")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.groundingline.melt_interpolation")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.masstransport.hydrostatic_adjustment")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.materials.rheology_law")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.damage.elementinterp")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.mesh.domain_type")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else if(strcmp(record_name,"md.mesh.elementtype")==0){
							this->AddConstant(new IoConstant(StringToEnumx(string),record_name));
						} else {
							/*Add string to parameters: */
							this->AddConstant(new IoConstant(string,record_name));
						}

						/*Free string*/
						xDelete<char>(string);
						break;
					case 5: break; //do nothing. not interested in this type of data, which is memory intensive.
					case 6: break; //do nothing. not interested in this type of data, which is memory intensive.
					case 7: break; //do nothing. not interested in this type of data, which is memory intensive.
					case 8: break; //do nothing. not interested in this type of data, which is memory intensive.
					case 10: break; //do nothing. not interested in this type of data, which is memory intensive.
					case 9:
							  ISSM_MPI_Bcast(&numstrings,1,ISSM_MPI_INT,0,IssmComm::GetComm());
							  /*Now allocate string array: */
							  if(numstrings){
								  strings=xNew<char*>(numstrings);
								  for(int i=0;i<numstrings;i++)strings[i]=NULL;

								  /*Go through strings, and read: */
								  for(int i=0;i<numstrings;i++){

									  ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
									  if(string_size){
										  string=xNew<char>((string_size+1));
										  string[string_size]='\0';
										  ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
									  }
									  else{
										  string=xNew<char>(1);
										  string[0]='\0';
									  }
									  strings[i]=string;
								  }
							  }

							  /*Add strings to parameters: */
							  this->AddConstant(new IoConstant(strings,numstrings,record_name));

							  /*Free string*/
							  for(int i=0;i<numstrings;i++) xDelete<char>(strings[i]);
							  xDelete<char*>(strings);
							  break;
					default:
							  _error_("unknown record type:" << record_code);
							  break;
				}
				xDelete<char>(record_name);
			}
		}
	} //}}}
}/*}}}*/
void  IoModel::FetchData(bool* pboolean,const char* data_name){/*{{{*/

	/*output: */
	int   booleanint;
	int   code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=1)_error_("expecting a boolean for \"" << data_name<<"\"");

	/*We have to read a boolean from disk. */
	if(my_rank==0){
		if(fread(&booleanint,sizeof(int),1,fid)!=1) _error_("could not read boolean ");
	}
	ISSM_MPI_Bcast(&booleanint,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*cast to bool: */
	/*Assign output pointers: */
	*pboolean=(bool)booleanint;

}
/*}}}*/
void  IoModel::FetchData(int* pinteger,const char* data_name){/*{{{*/

	/*output: */
	int   integer;
	int   code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();
	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=2)_error_("expecting an integer for \"" << data_name<<"\"");

	/*We have to read a integer from disk. First read the dimensions of the integer, then the integer: */
	if(my_rank==0){
		if(fread(&integer,sizeof(int),1,fid)!=1) _error_("could not read integer ");
	}

	ISSM_MPI_Bcast(&integer,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pinteger=integer;
}
/*}}}*/
void  IoModel::FetchData(IssmDouble* pscalar,const char* data_name){/*{{{*/

	/*output: */
	IssmPDouble   scalar;
	int      code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=3)_error_("expecting a IssmDouble for \""<<data_name<<"\"");

	/*We have to read a scalar from disk. First read the dimensions of the scalar, then the scalar: */
	if(my_rank==0){
		if(fread(&scalar,sizeof(IssmPDouble),1,fid)!=1)_error_("could not read scalar ");
	}
	ISSM_MPI_Bcast(&scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

	/*Assign output pointers: */
	*pscalar=scalar;

}
/*}}}*/
void  IoModel::FetchData(IssmDouble** pscalar, const char* data_name){/*{{{*/

   /*output: */
   IssmPDouble *scalar = NULL;
   int          code   = 0;

   /*recover my_rank:*/
   int my_rank=IssmComm::GetRank();

   /*Set file pointer to beginning of the data: */
   fid=this->SetFilePointerToData(&code,NULL,data_name);
   if(code!=3)_error_("expecting a IssmDouble for \""<<data_name<<"\"");

   /*Now fetch: */

   /*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */

   /*Now allocate matrix: */
   /*Read matrix on node 0, then broadcast: */
   scalar=xNew<IssmPDouble>(1);
   if(my_rank==0) if(fread(scalar,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
   ISSM_MPI_Bcast(scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

   _printf0_("scalar: " << *scalar << "\n");
   *pscalar=xNew<IssmDouble>(1);
   *pscalar[0]=scalar[0];
   xDelete<IssmPDouble>(scalar);
}
/*}}}*/
void  IoModel::FetchData(char** pstring,const char* data_name){/*{{{*/

	/*output: */
	char* string=NULL;
	int   string_size;
	int code=0;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=4)_error_("expecting a string for \""<<data_name<<"\"");

	/*Now fetch: */

	/*We have to read a string from disk. First read the dimensions of the string, then the string: */
	if(my_rank==0){
		if(fread(&string_size,sizeof(int),1,fid)!=1) _error_("could not read length of string ");
	}

	ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate string: */
	if(string_size){
		string=xNew<char>((string_size+1));
		string[string_size]='\0';

		/*Read string on node 0, then broadcast: */
		if(my_rank==0){
			if(fread(string,string_size*sizeof(char),1,fid)!=1)_error_(" could not read string ");
		}
		ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
	}
	else{
		string=xNew<char>(1);
		string[0]='\0';
	}

	/*Assign output pointers: */
	*pstring=string;
}
/*}}}*/
void  IoModel::FetchData(char*** pstrings,int* pnumstrings,const char* data_name){/*{{{*/

	/*output: */
	char** strings = NULL;
	char*  string  = NULL;
	int    numstrings;
	int    string_size;
	int    code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=9)_error_("expecting a string array for \""<<data_name<<"\"");

	/*Now fetch: */

	if(my_rank==0){
		if(fread(&numstrings,sizeof(int),1,fid)!=1) _error_("could not read length of string array");
	}
	ISSM_MPI_Bcast(&numstrings,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate string array: */
	if(numstrings){
		strings=xNew<char*>(numstrings);
		for(int i=0;i<numstrings;i++) strings[i]=NULL;

		/*Go through strings, and read: */
		for(int i=0;i<numstrings;i++){

			if(my_rank==0){
				if(fread(&string_size,sizeof(int),1,fid)!=1) _error_("could not read length of string ");
			}
			ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
			if(string_size){
				string=xNew<char>((string_size+1));
				string[string_size]='\0';
				if(my_rank==0){
					if(fread(string,string_size*sizeof(char),1,fid)!=1)_error_(" could not read string ");
				}
				ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
			}
			else{
				string=xNew<char>(1);
				string[0]='\0';
			}
			strings[i]=string;
		}
	}

	/*Assign output pointers: */
	*pstrings = strings;
	if(pnumstrings) *pnumstrings = numstrings;
}
/*}}}*/
void  IoModel::FetchData(bool** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/

	/*output: */
	int M,N;
	IssmPDouble* matrix=NULL;
	bool*        bool_matrix=NULL;
	int code=0;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=5 && code!=6 && code!=7)_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")");

	/*Now fetch: */

	/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
	/*numberofelements: */
	if(my_rank==0){
		if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
	}

	ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(my_rank==0){
		if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
	}
	ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate matrix: */
	if(M*N){
		matrix=xNew<IssmPDouble>(M*N);

		/*Read matrix on node 0, then broadcast: */
		if(my_rank==0){
			if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
		}

		ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
	}

	/*Now cast to bool: */
	if(M*N){
		bool_matrix=xNew<bool>(M*N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				bool_matrix[i*N+j]=(bool)matrix[i*N+j];
			}
		}
	}
	else{
		bool_matrix=NULL;
	}
	/*Free resources:*/
	xDelete<IssmPDouble>(matrix);

	/*Assign output pointers: */
	*pmatrix=bool_matrix;
	if (pM)*pM=M;
	if (pN)*pN=N;

}
/*}}}*/
void  IoModel::FetchData(int** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/
	int i,j;

	/*output: */
	int M,N;
	IssmPDouble* matrix=NULL;
	int*    integer_matrix=NULL;
	int code=0;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);

	if(code!=5 && code!=6 && code!=7)_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")");

	/*Now fetch: */

	/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
	/*numberofelements: */
	if(my_rank==0){
		if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
	}

	ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(my_rank==0){
		if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
	}
	ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate matrix: */
	if(M*N){
		matrix=xNew<IssmPDouble>(M*N);

		/*Read matrix on node 0, then broadcast: */
		if(my_rank==0){
			if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
		}

		ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
	}

	/*Now cast to integer: */
	if(M*N){
		integer_matrix=xNew<int>(M*N);
		for (i=0;i<M;i++){
			for (j=0;j<N;j++){
				integer_matrix[i*N+j]=(int)matrix[i*N+j];
			}
		}
	}
	else{
		integer_matrix=NULL;
	}
	/*Free resources:*/
	xDelete<IssmPDouble>(matrix);

	/*Assign output pointers: */
	*pmatrix=integer_matrix;
	if (pM)*pM=M;
	if (pN)*pN=N;

}
/*}}}*/
void  IoModel::FetchData(IssmDouble** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/

	/*First, look if has already been loaded (might be an independent variable)*/
	vector<IoData*>::iterator iter;
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,data_name)==0){
			*pmatrix=iodata->data;
			if(pM) *pM=iodata->M;
			if(pN) *pN=iodata->N;
			return;
		}
	}

	/*output: */
	int          M,N;
	IssmPDouble *matrix = NULL;
	int          code   = 0;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);
	if(code!=5 && code!=6 && code!=7 && code!=10)_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")");

	/*Now fetch: */

	/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
	/*numberofelements: */
	if(my_rank==0){
		if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
	}
	ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(my_rank==0){
		if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
	}
	ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate matrix: */
	if(M*N){
		if(code==10){
			/*Special case for Compressed mat*/
			IssmPDouble offset,range;
			if(my_rank==0) if(fread(&offset,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read offset");
			ISSM_MPI_Bcast(&offset,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			if(my_rank==0) if(fread(&range,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read range");
			ISSM_MPI_Bcast(&range,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			*pmatrix=xNew<IssmDouble>(M*N);

			/*Read matrix*/
			uint8_t* rawmatrix=xNew<uint8_t>((M-1)*N);
			if(my_rank==0) if(fread(rawmatrix,(M-1)*N*sizeof(char),1,fid)!=1) _error_("could not read matrix ");
			ISSM_MPI_Bcast(rawmatrix,(M-1)*N,ISSM_MPI_CHAR,0,IssmComm::GetComm());

			for(int i=0;i<(M-1)*N;++i) (*pmatrix)[i]=offset+range*reCast<IssmDouble>(rawmatrix[i])/255.;
			xDelete<uint8_t>(rawmatrix);

			/*read time now*/
			IssmPDouble* timematrix=xNew<IssmPDouble>(N);
			if(my_rank==0) if(fread(timematrix,N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read time in compressed matrix");
			ISSM_MPI_Bcast(timematrix,N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			for(int i=0;i<N;++i) (*pmatrix)[(M-1)*N+i]=timematrix[i];
			xDelete<IssmPDouble>(timematrix);

		}
		else{
			/*Read matrix on node 0, then broadcast: */
			matrix=xNew<IssmPDouble>(M*N);
			if(my_rank==0) if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
			ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			*pmatrix=xNew<IssmDouble>(M*N);
			for(int i=0;i<M*N;++i) (*pmatrix)[i]=matrix[i];
			xDelete<IssmPDouble>(matrix);
		}
	}
	else{
		*pmatrix=NULL;
	}
	/*Assign output pointers: */
	if(pM) *pM=M;
	if(pN) *pN=N;
}
/*}}}*/
#ifdef _HAVE_AD_
void  IoModel::FetchData(IssmPDouble** pmatrix,int* pM,int* pN,const char* data_name){/*{{{*/

	/*First, look if has already been loaded (might be an independent variable)*/
	vector<IoData*>::iterator iter;
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,data_name)==0){
			_error_(data_name <<" has already been loaded as in IssmDouble and cannot be converted to IssmPDouble");
		}
	}

	/*output: */
	int          M,N;
	IssmPDouble *matrix = NULL;
	int          code   = 0;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);
	if(code!=5 && code!=6 && code!=7 && code!=10)_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\""<<" (Code is "<<code<<")");

	/*Now fetch: */

	/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
	/*numberofelements: */
	if(my_rank==0){
		if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
	}
	ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(my_rank==0){
		if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
	}
	ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate matrix: */
	if(M*N){
		if(code==10){
			/*Special case for Compressed mat*/
			IssmPDouble offset,range;
			if(my_rank==0) if(fread(&offset,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read offset");
			ISSM_MPI_Bcast(&offset,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			if(my_rank==0) if(fread(&range,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read range");
			ISSM_MPI_Bcast(&range,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			*pmatrix=xNew<IssmPDouble>(M*N);

			/*Read matrix*/
			uint8_t* rawmatrix=xNew<uint8_t>((M-1)*N);
			if(my_rank==0) if(fread(rawmatrix,(M-1)*N*sizeof(char),1,fid)!=1) _error_("could not read matrix ");
			ISSM_MPI_Bcast(rawmatrix,(M-1)*N,ISSM_MPI_CHAR,0,IssmComm::GetComm());

			for(int i=0;i<(M-1)*N;++i) (*pmatrix)[i]=offset+range*reCast<IssmPDouble>(rawmatrix[i])/255.;
			xDelete<uint8_t>(rawmatrix);

			/*read time now*/
			IssmPDouble* timematrix=xNew<IssmPDouble>(N);
			if(my_rank==0) if(fread(timematrix,N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read time in compressed matrix");
			ISSM_MPI_Bcast(timematrix,N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			for(int i=0;i<N;++i) (*pmatrix)[(M-1)*N+i]=timematrix[i];
			xDelete<IssmPDouble>(timematrix);

		}
		else{
			/*Read matrix on node 0, then broadcast: */
			matrix=xNew<IssmPDouble>(M*N);
			if(my_rank==0) if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
			ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
			*pmatrix=matrix;
		}
	}
	else{
		*pmatrix=NULL;
	}
	/*Assign output pointers: */
	if(pM) *pM=M;
	if(pN) *pN=N;
}
/*}}}*/
void  IoModel::FetchData(IssmPDouble** pscalar,const char* data_name){/*{{{*/

   /*output: */
   IssmPDouble   *scalar = NULL;
   int      code;

   /*recover my_rank:*/
   int my_rank=IssmComm::GetRank();

   /*Set file pointer to beginning of the data: */
   fid=this->SetFilePointerToData(&code,NULL,data_name);

   if(code!=3)_error_("expecting a IssmDouble for \""<<data_name<<"\"");

   /*We have to read a scalar from disk. First read the dimensions of the scalar, then the scalar: */
   scalar=xNew<IssmPDouble>(1);
   if(my_rank==0){
      if(fread(scalar,sizeof(IssmPDouble),1,fid)!=1)_error_("could not read scalar ");
   }
   ISSM_MPI_Bcast(scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

   /*Assign output pointers: */
   *pscalar=scalar;
}
/*}}}*/
#endif
void  IoModel::FetchData(IssmDouble*** pmatrices,int** pmdims,int** pndims, int* pnumrecords,const char* data_name){/*{{{*/

	int i;
	/*output: */
	IssmDouble** matrices=NULL;
	int*     mdims=NULL;
	int*     ndims=NULL;
	int      numrecords=0;

	/*intermediary: */
	int     M, N;
	IssmPDouble *matrix = NULL;
	int     code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);
	if(code!=8)_error_("expecting a IssmDouble mat array for \""<<data_name<<"\"");

	/*Now fetch: */
	if(my_rank==0){
		if(fread(&numrecords,sizeof(int),1,fid)!=1) _error_("could not read number of records in matrix array ");
	}
	ISSM_MPI_Bcast(&numrecords,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(numrecords){

		/*Allocate matrices :*/
		matrices=xNew<IssmDouble*>(numrecords);
		mdims=xNew<int>(numrecords);
		ndims=xNew<int>(numrecords);

		for(i=0;i<numrecords;i++){
			matrices[i]=NULL;
			mdims[i]=0;
			ndims[i]=0;
		}

		/*Loop through records and fetch matrix: */
		for(i=0;i<numrecords;i++){

			if(my_rank==0){
				if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows in " << i << "th matrix of matrix array");
			}
			ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			if(my_rank==0){
				if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns in " << i << "th matrix of matrix array");
			}
			ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			/*Now allocate matrix: */
			if(M*N){
				matrix=xNew<IssmPDouble>(M*N);

				/*Read matrix on node 0, then broadcast: */
				if(my_rank==0){
					if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
				}

				ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
				matrices[i]=xNew<IssmDouble>(M*N);
				for (int j=0;j<M*N;++j) {matrices[i][j]=matrix[j];}
				xDelete<IssmPDouble>(matrix);
			}
			else
			  matrices[i]=NULL;
			/*Assign: */
			mdims[i]=M;
			ndims[i]=N;
		}
	}

	/*Assign output pointers: */
	*pmatrices=matrices;
	*pmdims=mdims;
	*pndims=ndims;
	*pnumrecords=numrecords;
}
/*}}}*/
void  IoModel::FetchData(IssmDouble** pmatrix,int* pM,int* pN, int layer_number,const char* data_name){/*{{{*/
	/*Same function as above, but here we want to fetch only one specific "layer" of the dataset to avoid memory issues*/

	/*output: */
	int          M=0,N=0;
	IssmPDouble *matrix = NULL;

	/*intermediary: */
	int  numrecords=0;
	int  code,i;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,data_name);
	if(code!=8)_error_("expecting a IssmDouble mat array for \""<<data_name<<"\"");

	if(my_rank==0){
		if(fread(&numrecords,sizeof(int),1,fid)!=1) _error_("could not read number of records in matrix array ");
	}
	ISSM_MPI_Bcast(&numrecords,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(numrecords){

		/*Check consistency of layer number*/
		_assert_(layer_number>=0);
		if(layer_number>=numrecords) _error_("layer number of "<<data_name<<" cannot exceed "<< numrecords-1);

		/*Skip until we get to the right layer*/
		if(my_rank==0){
			for(i=0;i<layer_number;i++){
				if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows in " << i << "th matrix of matrix array");
				if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns in " << i << "th matrix of matrix array");
				fseek(fid,M*N*sizeof(double),SEEK_CUR);
			}
		}

		/*fetch this slice!*/
		if(my_rank==0){
			if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows in " << i << "th matrix of matrix array");
		}
		ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

		if(my_rank==0){
			if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns in " << i << "th matrix of matrix array");
		}
		ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

		/*Now allocate matrix: */
		if(M*N){
			matrix=xNew<IssmPDouble>(M*N);

			/*Read matrix on node 0, then broadcast: */
			if(my_rank==0) if(fread(matrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
			ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			*pmatrix=xNew<IssmDouble>(M*N);
			for(int i=0;i<M*N;++i) (*pmatrix)[i]=matrix[i];
			xDelete<IssmPDouble>(matrix);
		}
		else{
			*pmatrix=NULL;
		}
	}

	/*Assign output pointers: */
	if(pM) *pM = M;
	if(pN) *pN = N;
}
/*}}}*/
void  IoModel::FetchData(Options* options,const char* lastnonoption){/*{{{*/

	/*record descriptions; */
	const char* mddot = "md.";
	char* record_name = NULL;
	int   record_name_size;
	long long record_length;
	int   record_code;

	/*records: */
	IssmDouble   scalar = 0;
	char        *string = NULL;
	int          string_size;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Go find in the binary file, the position of the data we want to fetch: */
	if(my_rank==0){
		fseek(fid,0,SEEK_SET);
		for(;;){
			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0) _error_("could not read record_name");
			if(record_name_size<3 || record_name_size>80) _error_("error while looking in binary file. Found a string of size "<<record_name_size);

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0)_error_("Could not find field "<<lastnonoption);
			if(strncmp(record_name,mddot,3)!=0) _error_("error while reading binary file: record does not start with \"md.\": "<<record_name);

			/*Is this the record sought for? : */
			if(strcmp(record_name,lastnonoption)==0){
				if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
				fseek(fid,record_length,SEEK_CUR);
				xDelete<char>(record_name);
				break;
			}
			else{
				if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
				fseek(fid,record_length,SEEK_CUR);
				xDelete<char>(record_name);
			}
		}
	}

	/*Go find in the binary file, the position of the data we want to fetch: */
	if(my_rank==0){ //cpu 0{{{

		/*Now march through file looking for the correct data identifiers (bool,int,IssmDouble or string): */
		for(;;){

			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				record_code=0; //0 means bailout
				ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*tell others cpus we are bailing: */
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				_error_("error while looking in binary file. Found a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				_error_("Could not read record name");
			}
			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}
			if(strcmp(record_name,"md.EOF")==0){
				xDelete<char>(record_name);
				record_code=0; //0 means bailout
				ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*tell others cpus we are bailing: */
				break;
			}

			/* Read the record length and the data type code: */
			if(fread(&record_length,sizeof(long long),1,this->fid)!=1) _error_("Cound not read record_length");
			if(fread(&record_code  ,sizeof(int),1,this->fid)!=1) _error_("Cound not read record_code");

			/*Tell other cpus what we are doing: */
			ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*tell other cpus what we are going to do: */

			/*Tell other cpus the name of the data, then branch according to the data type: */
			ISSM_MPI_Bcast(&record_name_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
			ISSM_MPI_Bcast(record_name,record_name_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
			ISSM_MPI_Bcast(&record_length,1,ISSM_MPI_LONG_LONG_INT,0,IssmComm::GetComm());

			switch(record_code){
				case 3:
					  {
						if(fread(&scalar,sizeof(IssmPDouble),1,this->fid)!=1) _error_("could not read scalar ");
						ISSM_MPI_Bcast(&scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
						GenericOption<IssmDouble>* option = new GenericOption<IssmDouble>();
						char* optionname=xNew<char>(strlen(record_name)-3+1);
						xMemCpy(optionname,&record_name[3],strlen(record_name)-3+1);
						option->value = scalar;
						option->name  = optionname;
						option->size[0] = 1;
						option->size[1] = 1;
						options->AddOption(option);
					  }
					break;
				case 4:
					  {
					/*We have to read a string from disk. First read the dimensions of the string, then the string: */
					if(fread(&string_size,sizeof(int),1,this->fid)!=1) _error_("could not read length of string ");
					ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

					if(string_size){
						string=xNew<char>(string_size+1);
						string[string_size]='\0';

						/*Read string, then broadcast: */
						if(fread(string,string_size*sizeof(char),1,this->fid)!=1)_error_(" could not read string ");
						ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
					}
					else{
						string=xNew<char>(1);
						string[0]='\0';
					}

					/*Add string to parameters: */
					GenericOption<char*>* option = new GenericOption<char*>();
					char* optionname=xNew<char>(strlen(record_name)-3+1);
					xMemCpy(optionname,&record_name[3],strlen(record_name)-3+1);
					option->value = string;
					option->name  = optionname;
					option->size[0] = 1;
					option->size[1] = 1;
					options->AddOption(option);

					  }
					break;
				default:
					_error_("record type not supported:" << record_code);
					break;
			}
			xDelete<char>(record_name);
		}
	} //}}}
	else{ //cpu ~0 {{{
		for(;;){ //wait on cpu 0
			ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());  /*get from cpu 0 what we are going to do: */
			if(record_code==0){
				break; //we are done, break from the loop
			}
			else{
				ISSM_MPI_Bcast(&record_name_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());
				_assert_(record_name_size);
				record_name=xNew<char>((record_name_size+1)); record_name[record_name_size]='\0';
				ISSM_MPI_Bcast(record_name,record_name_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
				ISSM_MPI_Bcast(&record_length,1,ISSM_MPI_LONG_LONG_INT,0,IssmComm::GetComm());
				switch(record_code){
					case 3:
						  {
							if(fread(&scalar,sizeof(IssmPDouble),1,this->fid)!=1) _error_("could not read scalar ");
							ISSM_MPI_Bcast(&scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());
							char* optionname=xNew<char>(strlen(record_name)-3+1);
							xMemCpy(optionname,&record_name[3],strlen(record_name)-3+1);
							GenericOption<IssmDouble>* option = new GenericOption<IssmDouble>();
							option->value = scalar;
							option->name  = optionname;
							option->size[0] = 1;
							option->size[1] = 1;
							options->AddOption(option);
						  }
						break;
					case 4:
						  {
						/*We have to read a string from disk. First read the dimensions of the string, then the string: */
						if(fread(&string_size,sizeof(int),1,this->fid)!=1) _error_("could not read length of string ");
						ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

						if(string_size){
							string=xNew<char>(string_size+1);
							string[string_size]='\0';

							/*Read string, then broadcast: */
							if(fread(string,string_size*sizeof(char),1,this->fid)!=1)_error_(" could not read string ");
							ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
						}
						else{
							string=xNew<char>(1);
							string[0]='\0';
						}

						/*Add string to parameters: */
						char* optionname=xNew<char>(strlen(record_name)-3+1);
						xMemCpy(optionname,&record_name[3],strlen(record_name)-3+1);
						GenericOption<char*>* option = new GenericOption<char*>();
						option->value = string;
						option->name  = optionname;
						option->size[0] = 1;
						option->size[1] = 1;
						options->AddOption(option);
						  }
						break;
					default:
						_error_("record type not supported:" << record_code);
						break;
				}

			}
		}
	} //}}}

}
/*}}}*/
void  IoModel::FetchData(int num,...){/*{{{*/

	va_list     ap;
	int         code,layout;
	IssmDouble *matrix   = NULL;
	char*       data_name;
	int         M,N;
	bool        exists;
	const char *mddot     = "md.";
	vector<IoData*>::iterator iter;

	/*Go through the entire list of names and fetch the corresponding data. Add it to the iomodel->data dataset. Everything
	 *we fetch is a IssmDouble* : */

	va_start(ap,num);
	for(int i=0; i<num; i++){

		data_name=va_arg(ap,char*);
		if(strncmp(data_name,mddot,3)!=0) _error_("String provided does not start with \"md.\" ("<<data_name<<")");

		exists = false;

		for(iter=data.begin();iter<data.end();iter++){
			IoData* iodata=*iter;
			if(strcmp(iodata->name,data_name)==0){
				/*Already there, no need to fetch it*/
				_assert_(iodata->isindependent);
				exists = true;
				break;
			}
		}

		if(exists){
			/*this data has already been checked out! Continue: */
			continue;
		}
		else{
			/*Add to this->data: */
			this->SetFilePointerToData(&code,&layout,data_name);
			this->FetchData(&matrix,&M,&N,data_name);
			this->AddData(new IoData(matrix,code,layout,M,N,data_name));
		}
	}
	va_end(ap);
}
/*}}}*/
void  IoModel::FetchDataToInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum,IssmDouble default_value){/*{{{*/

	/*First, look whether it is not already loaded in this->data*/
	vector<IoData*>::iterator iter;
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,vector_name)==0){
			_assert_(iodata->code==7);
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->InputCreate(iodata->data,inputs,this,iodata->M,iodata->N,iodata->layout,input_enum,iodata->code);//we need i to index into elements.
			}
			return;
		}
	}

	/*intermediary: */
	int         code,vector_layout;
	IssmDouble *doublearray = NULL;
	int         M,N;

	/*First of, find the record for the name, and get code  of data type: */
	this->SetFilePointerToData(&code, &vector_layout,vector_name);

	/*Defaulting only supported for double arrays*/
	if(code!=7 && code!=10) _error_(vector_name<<" is not a double array");

	this->FetchData(&doublearray,&M,&N,vector_name);

	for(Object* & object : elements->objects){
		Element* element=xDynamicCast<Element*>(object);
		if(!doublearray){
			element->SetElementInput(inputs,input_enum,default_value);
		}
		else{
			element->InputCreate(doublearray,inputs,this,M,N,vector_layout,input_enum,code);//we need i to index into elements.
		}
	}

	/*Free resources:*/
	xDelete<IssmDouble>(doublearray);
}
/*}}}*/
void  IoModel::FetchDataToInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum){/*{{{*/

	/*First, look whether it is not already loaded in this->data*/
	vector<IoData*>::iterator iter;
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,vector_name)==0){
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->InputCreate(iodata->data,inputs,this,iodata->M,iodata->N,iodata->layout,input_enum,iodata->code);//we need i to index into elements.
			}
			return;
		}
	}

	/*intermediary: */
	int     i;
	int     code,vector_layout;

	/*variables being fetched: */
	bool        boolean;
	int         integer;
	IssmDouble  scalar;
	IssmDouble *doublearray = NULL;
	int         M,N;

	/*First of, find the record for the name, and get code  of data type: */
	this->SetFilePointerToData(&code, &vector_layout,vector_name);

	switch(code){
		case 1: //boolean constant
			this->FetchData(&boolean,vector_name);
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->SetBoolInput(inputs,input_enum,boolean);
			}
			break;
		case 2: //integer constant
			this->FetchData(&integer,vector_name);
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->SetIntInput(inputs,input_enum,integer);
			}
			break;
		case 3: //IssmDouble constant
			this->FetchData(&scalar,vector_name);
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->SetElementInput(inputs,input_enum,scalar);
			}
			break;
		case 5: //boolean vector
			this->FetchData(&doublearray,&M,&N,vector_name); //we still have a doublearray, because it might include times in transient mode
			if(!doublearray) _error_("\""<<vector_name<<"\" not found in binary file");
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->InputCreate(doublearray,inputs,this,M,N,vector_layout,input_enum,code);//we need i to index into elements.
			}
			break;
		case 6: //int vector
			this->FetchData(&doublearray,&M,&N,vector_name); //we still have a doublearray, because it might include times in transient mode
			if(!doublearray) _error_("\""<<vector_name<<"\" not found in binary file");
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->InputCreate(doublearray,inputs,this,M,N,vector_layout,input_enum,code);//we need i to index into elements.
			}
			break;
		case 8: { //MatArray {{{

			/*variables:*/
			int numarray;
			IssmDouble** array=NULL; 
			IssmDouble*  matrix=NULL;
			IssmDouble*  times = NULL;
			int* pM = NULL;
			int* pN = NULL;
			int M,N;

			/*fetch array of matrices:*/
			this->FetchData(&array,&pM,&pN,&numarray,vector_name);

			for (int i=0;i<numarray;i++){

				M=pM[i];
				N=pN[i];
				matrix=array[i];

				//initialize times:
				if(M==this->numberofvertices || M==(this->numberofvertices+1)){
					times=xNew<IssmDouble>(N);
					if(M==this->numberofvertices) times[0] = matrix[M-1];
					if(M==this->numberofvertices+1) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];
				}
				else if(M==this->numberofelements || M==(this->numberofelements+1)){
					times=xNew<IssmDouble>(N);
					if(M==this->numberofelements) times[0] = matrix[M-1];
					if(M==this->numberofelements+1) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];
				}
				else if(M==2 || M==1){
					times=xNew<IssmDouble>(N);
					if(M==1) times[0] = 0;
					if(M==2) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];
				}
				else _error_("FetchDataToInput error message: row size of MatArray elements should be either numberofelements (+1) or numberofvertices (+1)");


				//initialize transient input dataset:
				TransientInput* transientinput=inputs->SetDatasetTransientInput(input_enum,i, times,N);
				for(Object* & object : elements->objects){

					/*Get the right transient input*/
					Element* element=xDynamicCast<Element*>(object);

					/*Get values and lid list*/
					const int   numvertices = element->GetNumberOfVertices();
					int        *vertexlids = xNew<int>(numvertices);
					int        *vertexsids = xNew<int>(numvertices);

					/*Recover vertices ids needed to initialize inputs*/
					_assert_(this->elements);
					for(int k=0;k<numvertices;k++){
						vertexsids[k] =reCast<int>(this->elements[numvertices*element->Sid()+k]-1); //ids for vertices are in the elements array from Matlab
						vertexlids[k]=this->my_vertices_lids[vertexsids[k]];
					}

					if(M==this->numberofvertices || M==(this->numberofvertices+1)){

						//recover time vector: 
						times=xNew<IssmDouble>(N);
						if(M==this->numberofvertices) times[0] = matrix[M-1];
						if(M==this->numberofvertices+1) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];

						IssmDouble* values=xNew<IssmDouble>(numvertices);

						for(int t=0;t<N;t++){
							for (int k=0;k<numvertices;k++)values[k]=matrix[N*vertexsids[k]+t];

							switch(element->ObjectEnum()){
								case TriaEnum:  transientinput->AddTriaTimeInput( t,numvertices,vertexlids,values,P1Enum); break;
								case PentaEnum: transientinput->AddPentaTimeInput(t,numvertices,vertexlids,values,P1Enum); break;
								default: _error_("Not implemented yet");
							}
						}
						xDelete<IssmDouble>(values);
					}
					else if(M==this->numberofelements || M==(this->numberofelements+1)){

						IssmDouble value;

						//recover time vector: 
						times=xNew<IssmDouble>(N);
						if(M==this->numberofelements) times[0] = matrix[M-1];
						if(M==this->numberofelements+1) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];

						for(int t=0;t<N;t++){ 

							value=matrix[N*element->Sid()+t];
							switch(element->ObjectEnum()){
								case TriaEnum:  transientinput->AddTriaTimeInput( t,1,&element->lid,&value,P0Enum); break;
								case PentaEnum:  transientinput->AddPentaTimeInput( t,1,&element->lid,&value,P0Enum); break;
								default: _error_("Not implemented yet");
							}
						}
					}
					else if(M==2 || M==1){
						IssmDouble value;

						//recover time vector: 
						times=xNew<IssmDouble>(N);
						if(M==1) times[0] = 0;
						if(M==2) for(int t=0;t<N;t++) times[t] = matrix[(M-1)*N+t];

						for(int t=0;t<N;t++){ 

							value=matrix[t];
							switch(element->ObjectEnum()){
								case TriaEnum:  transientinput->AddTriaTimeInput( t,1,&element->lid,&value,P0Enum); break;
								case PentaEnum:  transientinput->AddPentaTimeInput( t,1,&element->lid,&value,P0Enum); break;
								default: _error_("Not implemented yet");
							}
						}

					}
					else _error_("FetchDataToInput error message: row size of MatArray elements should be either numberofelements (+1) or numberofvertices (+1)");

					xDelete<int>(vertexlids);
					xDelete<int>(vertexsids);
				}

				xDelete<IssmDouble>(times);
			}

			/*Delete data:*/
			for(int i=0;i<numarray;i++){
				IssmDouble* matrix=array[i];
				xDelete<IssmDouble>(matrix);
			}
			xDelete<IssmDouble*>(array);
			xDelete<int>(pM);
			xDelete<int>(pN); 
			} //}}}
			break;
		case 7: //IssmDouble vector
		case 10:
			this->FetchData(&doublearray,&M,&N,vector_name);
			if(!doublearray) _error_("\""<<vector_name<<"\" not found in binary file");
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->InputCreate(doublearray,inputs,this,M,N,vector_layout,input_enum,code);//we need i to index into elements.
			}
			break;
		default:
			_error_("data code " << code << " not supported yet (detected while processing \""<<vector_name<<"\")");
			break;
	}
	/*Free resources:*/
	xDelete<IssmDouble>(doublearray);
}
/*}}}*/
void  IoModel::FetchDataToDatasetInput(Inputs* inputs,Elements* elements,const char* vector_name,int input_enum){/*{{{*/

	/*First, look whether it is not already loaded in this->data*/
	vector<IoData*>::iterator iter;
	for(iter=data.begin();iter<data.end();iter++){
		IoData* iodata=*iter;
		if(strcmp(iodata->name,vector_name)==0){
			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				_error_("to be implemented...");
				//element->InputCreate(iodata->data,inputs,this,iodata->M,iodata->N,iodata->layout,input_enum,iodata->code);//we need i to index into elements.
			}
			return;
		}
	}

	/*intermediary: */
	int         code,vector_layout;
	IssmDouble *doublearray = NULL;
	int         M,N;

	/*First of, find the record for the name, and get code  of data type: */
	this->SetFilePointerToData(&code,&vector_layout,vector_name);

	switch(code){
		case 1: //boolean constant
			_error_("not implemented yet");
			break;
		case 2: //integer constant
			_error_("not implemented yet");
			break;
		case 3: //IssmDouble constant
			_error_("not implemented yet");
			break;
		case 5: //boolean vector
			_error_("not implemented yet");
			break;
		case 6: //int vector
			_error_("not implemented yet");
			break;
		case 7: //IssmDouble vector
			  {
			this->FetchData(&doublearray,&M,&N,vector_name);
			if(!doublearray) _error_("\""<<vector_name<<"\" not found in binary file");

			int* ids = xNew<int>(N);
			for(int i=0;i<N;i++) ids[i] = i;

			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->DatasetInputCreate(doublearray,M,N,ids,N,inputs,this,input_enum);
			}
			xDelete<int>(ids);
			  }
			break;
		case 10: //Compressed matrix
			  {
			this->FetchData(&doublearray,&M,&N,vector_name);
			if(!doublearray) _error_("\""<<vector_name<<"\" not found in binary file");

			int* ids = xNew<int>(N);
			for(int i=0;i<N;i++) ids[i] = i;

			for(Object* & object : elements->objects){
				Element* element=xDynamicCast<Element*>(object);
				element->DatasetInputCreate(doublearray,M,N,ids,N,inputs,this,input_enum);
			}
			xDelete<int>(ids);
			  }
			break;

		default:
			_error_("data code " << code << " not supported yet (detected while processing \""<<vector_name<<"\")");
			break;
	}
	/*Free resources:*/
	xDelete<IssmDouble>(doublearray);
}
/*}}}*/
void  IoModel::FetchIndependentConstant(int* pXcount,IssmPDouble* X,const char* constant_name){/*{{{*/

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*recover Xcount if X is not NULL:*/
	int Xcount = 0;
	if(X) Xcount=*pXcount;

	#ifdef _HAVE_AD_ //cannot come here unless you are running AD mode, from DeclaredIndependents:

	/*output: */
	IssmPDouble  pscalar;
	IssmDouble   scalar; //same as pscalar, except it's an ADOLC independent variable
	int          code;

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,NULL,constant_name);
	if(code!=3) _error_("expecting a IssmDouble for \"" << constant_name<<"\"");

	/*We have to read a scalar from disk. First read the dimensions of the scalar, then the scalar: */
	if(my_rank==0){
		if(fread(&pscalar,sizeof(IssmPDouble),1,fid)!=1)_error_("could not read scalar ");

		/*Now, before we even broadcast this to other nodes, declare the scalar  as an independent variable!. If we
		 *have been supplied an X vector, use it instead of what we just read: */
		#if defined(_HAVE_CODIPACK_)
			if(X){
				scalar=X[Xcount];
			} else {
				scalar=pscalar;
			}
			codi_global.registerInput(scalar);
		#else
			if(X){
				scalar<<=X[Xcount];
			}
			else{
				scalar<<=pscalar;
			}
		#endif
	}

	ISSM_MPI_Bcast(&scalar,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
	this->AddConstantIndependent(new IoConstant(scalar,constant_name));

	/*increment offset into X vector, now that we have read 1 value:*/
	Xcount++; *pXcount=Xcount;
	#endif
}
/*}}}*/
void  IoModel::FetchIndependentData(int* pXcount,IssmPDouble* X,const char* data_name){/*{{{*/

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*recover Xcount if X is not NULL:*/
	int Xcount = 0;
	if(X) Xcount=*pXcount;

	#ifdef _HAVE_AD_ //cannot come here unless you are running AD mode, from DeclaredIndependents:

	/*Intermediaries*/
	int M,N;
	IssmPDouble* buffer=NULL; //a buffer to read the data from disk
	IssmDouble* matrix=NULL; //our independent variable
	int code,layout;

	/*Set file pointer to beginning of the data: */
	fid=this->SetFilePointerToData(&code,&layout,data_name);
	if((code!=5) && (code!=6) && (code!=7))_error_("expecting a IssmDouble, integer or boolean matrix for \"" << data_name<<"\"");

	/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
	/*numberofelements: */
	if(my_rank==0){
		if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
	}
	ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	if(my_rank==0){
		if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
	}
	ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Now allocate matrix: */
	if(M*N){
		buffer=xNew<IssmPDouble>(M*N);
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
		matrix=xNew<IssmDouble>(M*N,"t");
#else
		matrix=xNew<IssmDouble>(M*N);
#endif

		/*Read matrix on node 0, then broadcast: */
		if(my_rank==0){
			if(fread(buffer,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");

			/*Now, before we even broadcast this to other nodes, declare the whole matrix as a independent variable!
			  If we have been supplied an X vector, use it instead of what we just read: */
			#if defined(_HAVE_CODIPACK_)
				if(X){
					for (int i=0;i<M*N;i++) {
						matrix[i]=X[Xcount+i];
						codi_global.registerInput(matrix[i]);
					}
				}
				else{
					for (int i=0;i<M*N;i++) {
						matrix[i]=buffer[i];
						codi_global.registerInput(matrix[i]);
					}
				}
			#else /*ADOLC*/
				if(X){
					for(int i=0;i<M*N;i++) matrix[i]<<=X[Xcount+i];  /*<<= ADOLC overloaded operator to declare independent*/
				}
				else{
					for(int i=0;i<M*N;i++) matrix[i]<<=buffer[i];
				}
			#endif
		}
		ISSM_MPI_Bcast(matrix,M*N,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

		xDelete<IssmPDouble>(buffer);
	}
	else _error_("cannot declare the independent variable \"" << data_name <<  "\" if it's empty!");

	// FIXME codi is that at all relevant to CoDiPack or can we simply assume the same?

	/*Add to data as independent*/
	this->AddDataIndependent(new IoData(matrix,code,layout,M,N,data_name));

	/*increment offset into X vector, now that we have read M*N values:*/
	Xcount+=M*N; *pXcount=Xcount;
	#endif
}
/*}}}*/
void  IoModel::FetchMultipleData(char*** pstrings,int* pnumstrings,const char* data_name){/*{{{*/

	int  num_instances;

	/*output: */
	int    numstrings = 0;
	char **strings    = NULL;

	/*intermediary: */
	char   *string         = NULL;
	int     string_size;
	int    *codes          = NULL;
	int    *code           = NULL;
	fpos_t *file_positions = NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Get file pointers to beginning of the data (multiple instances of it): */
	file_positions=this->SetFilePointersToData(&codes,NULL,&num_instances,data_name);

	if(num_instances){
		strings=xNew<char*>(num_instances);

		for(int i=0;i<num_instances;i++){

			if(my_rank==0){
				/*check we are indeed finding a string, not something else: */
				if(codes[i]!=4)_error_("expecting a string for \""<<data_name<<"\" but code is "<<codes[i]<<" not 4");

				/*We have to read a string from disk. First read the dimensions of the string, then the string: */
				fsetpos(fid,file_positions+i);
				if(fread(&string_size,sizeof(int),1,fid)!=1) _error_("could not read length of string ");
			}

			ISSM_MPI_Bcast(&string_size,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			/*Now allocate string: */
			if(string_size){
				string=xNew<char>((string_size+1));
				string[string_size]='\0';

				/*Read string on node 0, then broadcast: */
				if(my_rank==0){
					if(fread(string,string_size*sizeof(char),1,fid)!=1)_error_(" could not read string ");
				}
				ISSM_MPI_Bcast(string,string_size,ISSM_MPI_CHAR,0,IssmComm::GetComm());
			}
			else{
				string=xNew<char>(1);
				string[0]='\0';
			}
			strings[i]=string;
		}
	}
	/*Free resources:*/
	xDelete<int>(codes);
	xDelete<fpos_t>(file_positions);

	/*Assign output pointers: */
	*pstrings=strings;
	if(pnumstrings) *pnumstrings=num_instances;
}
/*}}}*/
void  IoModel::FetchMultipleData(int** pvector, int* pnum_instances,const char* data_name){/*{{{*/

	int     num_instances;
	fpos_t* file_positions=NULL;

	/*output: */
	int* vector=NULL;

	/*intermediary: */
	int  integer;
	int *codes   = NULL;
	int  code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Get file pointers to beginning of the data (multiple instances of it): */
	file_positions=this->SetFilePointersToData(&codes,NULL,&num_instances,data_name);

	if(num_instances){

		/*Allocate vector :*/
		vector=xNew<int>(num_instances);

		for(int i=0;i<num_instances;i++){

			if(my_rank==0){
				code=codes[i];

				if(code!=2)_error_("expecting an integer for \""<<data_name<<"\"");

				/*We have to read a integer from disk. First read the dimensions of the integer, then the integer: */
				fsetpos(fid,file_positions+i);
				if(my_rank==0){
					if(fread(&integer,sizeof(int),1,fid)!=1) _error_("could not read integer ");
				}
			}
			ISSM_MPI_Bcast(&integer,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			/*Assign: */
			vector[i]=integer;
		}
	}

	/*Free resources:*/
	xDelete<fpos_t>(file_positions);
	xDelete<int>(codes);

	/*Assign output pointers: */
	*pvector=vector;
	if(pnum_instances) *pnum_instances=num_instances;
}
/*}}}*/
void  IoModel::FetchMultipleData(IssmDouble** pvector, int* pnum_instances,const char* data_name){/*{{{*/

	int     num_instances;
	fpos_t* file_positions=NULL;

	/*output: */
	IssmDouble* vector=NULL;

	/*intermediary: */
	IssmPDouble          scalar;
	int         *codes   = NULL;
	int          code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Get file pointers to beginning of the data (multiple instances of it): */
	file_positions=this->SetFilePointersToData(&codes,NULL,&num_instances,data_name);

	if(num_instances){

		/*Allocate vector :*/
		vector=xNew<IssmDouble>(num_instances);

		for(int i=0;i<num_instances;i++){

			if(my_rank==0){
				code=codes[i];

				if(code!=3)_error_("expecting a double for \""<<data_name<<"\"");

				/*We have to read a double from disk: */
				fsetpos(fid,file_positions+i);
				if(my_rank==0){
					if(fread(&scalar,sizeof(IssmPDouble),1,fid)!=1) _error_("could not read scalar ");
				}
			}
			ISSM_MPI_Bcast(&scalar,1,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

			/*Assign: */
			vector[i]=scalar;
		}
	}

	/*Free resources:*/
	xDelete<fpos_t>(file_positions);
	xDelete<int>(codes);

	/*Assign output pointers: */
	*pvector=vector;
	*pnum_instances=num_instances;
}
/*}}}*/
void  IoModel::FetchMultipleData(IssmDouble*** pmatrices,int** pmdims,int** pndims, int* pnumrecords,const char* data_name){/*{{{*/

	int     num_instances;
	fpos_t* file_positions=NULL;

	/*output: */
	IssmDouble **matrices = NULL;
	int         *mdims    = NULL;
	int         *ndims    = NULL;

	/*intermediary: */
	int          M, N;
	IssmPDouble *pmatrix = NULL;
	IssmDouble  *matrix  = NULL;
	int         *codes   = NULL;
	int          code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Get file pointers to beginning of the data (multiple instances of it): */
	file_positions=this->SetFilePointersToData(&codes,NULL,&num_instances,data_name);

	if(num_instances){

		/*Allocate matrices :*/
		matrices=xNew<IssmDouble*>(num_instances);
		mdims=xNew<int>(num_instances);
		ndims=xNew<int>(num_instances);

		for(int i=0;i<num_instances;i++){

			if(my_rank==0){
				code=codes[i];

				if((code!=5) && (code!=6) && (code!=7))_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\"");

				/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
				/*numberofelements: */
				fsetpos(fid,file_positions+i);
				if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
			}
			ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			if(my_rank==0){
				if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
			}
			ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			/*Now allocate matrix: */
			if(M*N){
				pmatrix=xNew<IssmPDouble>(M*N);

				/*Read matrix on node 0, then broadcast: */
				if(my_rank==0){
					if(fread(pmatrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
				}
				ISSM_MPI_Bcast(pmatrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

				//if(this->independents[data_enum]){ FIXME
				//	/*this data has already been checked out! So cancel all that we've done here, and return
				//	 * the data[data_enum] directly: */
				//	matrix=this->data[data_enum];
				//}
				//else{
					matrix=xNew<IssmDouble>(M*N);
					for (int i=0;i<M*N;++i) matrix[i]=pmatrix[i];
				//}
				xDelete<IssmPDouble>(pmatrix);
			}
			else
				matrix=NULL;

			/*Assign: */
			mdims[i]=M;
			matrices[i]=matrix;
			ndims[i]=N;
		}
	}

	/*Free resources:*/
	xDelete<fpos_t>(file_positions);
	xDelete<int>(codes);

	/*Assign output pointers: */
	*pmatrices=matrices;
	if(pmdims){
		*pmdims=mdims;
	}
	else{
		xDelete<int>(mdims);
	}
	if(pndims){
		*pndims=ndims;
	}
	else{
		xDelete<int>(ndims);
	}
	if(pnumrecords){
		*pnumrecords=num_instances;
	}
}
/*}}}*/
void  IoModel::FetchMultipleData(int*** pmatrices,int** pmdims,int** pndims, int* pnumrecords,const char* data_name){/*{{{*/

	int     num_instances;
	fpos_t* file_positions=NULL;

	/*output: */
	int        **matrices = NULL;
	int         *mdims    = NULL;
	int         *ndims    = NULL;

	/*intermediary: */
	int          M, N;
	IssmPDouble *pmatrix = NULL;
	IssmDouble  *matrix  = NULL;
	int         *integer_matrix=NULL;
	int         *codes   = NULL;
	int          code;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();

	/*Get file pointers to beginning of the data (multiple instances of it): */
	file_positions=this->SetFilePointersToData(&codes,NULL,&num_instances,data_name);

	if(num_instances){

		/*Allocate matrices :*/
		matrices=xNew<int*>(num_instances);
		mdims=xNew<int>(num_instances);
		ndims=xNew<int>(num_instances);

		for(int i=0;i<num_instances;i++){

			if(my_rank==0){
				code=codes[i];

				if((code!=5) && (code!=6) && (code!=7))_error_("expecting a IssmDouble, integer or boolean matrix for \""<<data_name<<"\"");

				/*We have to read a matrix from disk. First read the dimensions of the matrix, then the whole matrix: */
				/*numberofelements: */
				fsetpos(fid,file_positions+i);
				if(fread(&M,sizeof(int),1,fid)!=1) _error_("could not read number of rows for matrix ");
			}
			ISSM_MPI_Bcast(&M,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			if(my_rank==0){
				if(fread(&N,sizeof(int),1,fid)!=1) _error_("could not read number of columns for matrix ");
			}
			ISSM_MPI_Bcast(&N,1,ISSM_MPI_INT,0,IssmComm::GetComm());

			/*Now allocate matrix: */
			if(M*N){
				pmatrix=xNew<IssmPDouble>(M*N);
				integer_matrix=xNew<int>(M*N);

				/*Read matrix on node 0, then broadcast: */
				if(my_rank==0){
					if(fread(pmatrix,M*N*sizeof(IssmPDouble),1,fid)!=1) _error_("could not read matrix ");
				}
				ISSM_MPI_Bcast(pmatrix,M*N,ISSM_MPI_PDOUBLE,0,IssmComm::GetComm());

				//if(this->independents[data_enum]){ FIXME
				//	/*this data has already been checked out! So cancel all that we've done here, and return
				//	 * the data[data_enum] directly: */
				//	matrix=this->data[data_enum];
				//	for (int i=0;i<M*N;++i) integer_matrix[i]=reCast<int>(matrix[i]);
				//}
				//else{
					for (int i=0;i<M*N;++i) integer_matrix[i]=pmatrix[i];
				//}
				xDelete<IssmPDouble>(pmatrix);
			}
			else
				integer_matrix=NULL;

			/*Assign: */
			mdims[i]=M;
			matrices[i]=integer_matrix;
			ndims[i]=N;
		}
	}

	/*Free resources:*/
	xDelete<fpos_t>(file_positions);
	xDelete<int>(codes);

	/*Assign output pointers: */
	*pmatrices=matrices;
	if(pmdims){
		*pmdims=mdims;
	}
	else{
		xDelete<int>(mdims);
	}
	if(pndims){
		*pndims=ndims;
	}
	else{
		xDelete<int>(ndims);
	}
	*pnumrecords=num_instances;
}
/*}}}*/
void  IoModel::FillIndependents(IssmDouble* xp){/*{{{*/

	_assert_(xp);

	/*Initialize local num ind*/
	int local_num_ind = 0;

	/*Process constants*/
	for(vector<IoConstant*>::iterator iter=constants.begin();iter<constants.end();iter++){
		if((*iter)->isindependent){
			(*iter)->constant->GetParameterValue(&xp[local_num_ind]);
			local_num_ind += 1;
		}
	}

	/*Process data*/
	for(vector<IoData*>::iterator iter=data.begin();iter<data.end();iter++){
		if((*iter)->isindependent){
			for(int i=0;i<(*iter)->M*(*iter)->N;i++){
				xp[local_num_ind+i] = (*iter)->data[i];
			}
			local_num_ind += (*iter)->M*(*iter)->N;
		}
	}

	_assert_(local_num_ind == this->NumIndependents());
}
/*}}}*/
void  IoModel::FindConstant(bool* pvalue,const char* constant_name){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			if(ioconstant->constant->ObjectEnum()!=BoolParamEnum){
				_printf0_("=========================================================================\n");
				_printf0_(" Marshalled file is not consistent with compiled code                    \n");
				_printf0_("                                                                         \n");
				_printf0_("    This problem typically happens when two different versions of ISSM   \n");
				_printf0_("    are being used. Make sure that you are running the same version:     \n");
				_printf0_("    - to marshall the model (i.e., MATLAB/python interface)              \n");
				_printf0_("    - to run ISSM (i.e., the compiled code issm.exe)                     \n");
				_printf0_("                                                                         \n");
				_printf0_("=========================================================================\n\n");
				_error_("\""<< constant_name <<"\" cannot return a bool, it is a " << EnumToStringx(ioconstant->constant->ObjectEnum()));
			}
			ioconstant->constant->GetParameterValue(pvalue);
			return;
		}
	}

	for(vector<IoConstant*>::iterator iter=constants.begin();iter<constants.end();iter++) (*iter)->constant->Echo();
	_error_("Could not find constant \""<<constant_name<<"\"");
}
/*}}}*/
void  IoModel::FindConstant(int* pvalue,const char* constant_name){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			if(ioconstant->constant->ObjectEnum()!=IntParamEnum){
				_printf0_("=========================================================================\n");
				_printf0_(" Marshalled file is not consistent with compiled code                    \n");
				_printf0_("                                                                         \n");
				_printf0_("    This problem typically happens when two different versions of ISSM   \n");
				_printf0_("    are being used. Make sure that you are running the same version:     \n");
				_printf0_("    - to marshall the model (i.e., MATLAB/python interface)              \n");
				_printf0_("    - to run ISSM (i.e., the compiled code issm.exe)                     \n");
				_printf0_("                                                                         \n");
				_printf0_("=========================================================================\n\n");
				_error_("\""<< constant_name <<"\" cannot return an int, it is a " << EnumToStringx(ioconstant->constant->ObjectEnum()));
			}
			ioconstant->constant->GetParameterValue(pvalue);
			return;
		}
	}

	_error_("Could not find constant \""<<constant_name <<"\"");
}
/*}}}*/
void  IoModel::FindConstant(IssmDouble* pvalue,const char* constant_name){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			if(ioconstant->constant->ObjectEnum()!=DoubleParamEnum){
				_printf0_("=========================================================================\n");
				_printf0_(" Marshalled file is not consistent with compiled code                    \n");
				_printf0_("                                                                         \n");
				_printf0_("    This problem typically happens when two different versions of ISSM   \n");
				_printf0_("    are being used. Make sure that you are running the same version:     \n");
				_printf0_("    - to marshall the model (i.e., MATLAB/python interface)              \n");
				_printf0_("    - to run ISSM (i.e., the compiled code issm.exe)                     \n");
				_printf0_("                                                                         \n");
				_printf0_("=========================================================================\n\n");
				_error_("\""<< constant_name <<"\" cannot return a double, it is a " << EnumToStringx(ioconstant->constant->ObjectEnum()));
			}
			ioconstant->constant->GetParameterValue(pvalue);
			return;
		}
	}

	_error_("Could not find constant \""<<constant_name <<"\"");
}
/*}}}*/
void  IoModel::FindConstant(char** pvalue,const char* constant_name){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			if(ioconstant->constant->ObjectEnum()!=StringParamEnum){
				_printf0_("=========================================================================\n");
				_printf0_(" Marshalled file is not consistent with compiled code                    \n");
				_printf0_("                                                                         \n");
				_printf0_("    This problem typically happens when two different versions of ISSM   \n");
				_printf0_("    are being used. Make sure that you are running the same version:     \n");
				_printf0_("    - to marshall the model (i.e., MATLAB/python interface)              \n");
				_printf0_("    - to run ISSM (i.e., the compiled code issm.exe)                     \n");
				_printf0_("                                                                         \n");
				_printf0_("=========================================================================\n\n");
				_error_("\""<< constant_name <<"\" cannot return a string, it is a " << EnumToStringx(ioconstant->constant->ObjectEnum()));
			}
			ioconstant->constant->GetParameterValue(pvalue);
			return;
		}
	}

	_error_("Could not find constant \""<<constant_name <<"\"");
}
/*}}}*/
void  IoModel::FindConstant(char*** pvalue,int* psize,const char* constant_name){/*{{{*/

	/*Intermediary*/
	vector<IoConstant*>::iterator iter;

	for(iter=constants.begin();iter<constants.end();iter++){
		IoConstant* ioconstant=*iter;

		if(strcmp(ioconstant->name,constant_name)==0){
			if(ioconstant->constant->ObjectEnum()!=StringArrayParamEnum){
				_printf0_("=========================================================================\n");
				_printf0_(" Marshalled file is not consistent with compiled code                    \n");
				_printf0_("                                                                         \n");
				_printf0_("    This problem typically happens when two different versions of ISSM   \n");
				_printf0_("    are being used. Make sure that you are running the same version:     \n");
				_printf0_("    - to marshall the model (i.e., MATLAB/python interface)              \n");
				_printf0_("    - to run ISSM (i.e., the compiled code issm.exe)                     \n");
				_printf0_("                                                                         \n");
				_printf0_("=========================================================================\n\n");
				_error_("\""<< constant_name <<"\" cannot return a string array, it is a " << EnumToStringx(ioconstant->constant->ObjectEnum()));
			}
			ioconstant->constant->GetParameterValue(pvalue,psize);
			return;
		}
	}

	_error_("Could not find constant \""<<constant_name <<"\"");
}
/*}}}*/
int   IoModel::NumIndependents(void){/*{{{*/

	/*Initialize output*/
	int num_independents = 0;

	/*Process constants*/
	for(vector<IoConstant*>::iterator iter=constants.begin();iter<constants.end();iter++){
		if((*iter)->isindependent){
			num_independents+= 1;
		}
	}

	/*Process data*/
	for(vector<IoData*>::iterator iter=data.begin();iter<data.end();iter++){
		if((*iter)->isindependent){
			num_independents+= (*iter)->M*(*iter)->N;
		}
	}

	/*return*/
	return num_independents;
}
/*}}}*/
fpos_t* IoModel::SetFilePointersToData(int** pcodes,int** pvector_types, int* pnum_instances,const char* data_name){/*{{{*/

	int     found          = 0;
	const char* mddot = "md.";
	char* record_name = NULL;
	int   record_name_size;
	long long  record_length;
	int     record_code;           //1 to 7 number
	int     vector_type;           //1 to 7 number
	int    *vector_types   = NULL;
	int    *codes          = NULL;
	int     num_instances  = 0;
	int     counter;
	fpos_t *file_positions = NULL;

	/*recover my_rank:*/
	int my_rank=IssmComm::GetRank();
	_assert_(strncmp(data_name,mddot,3)==0);

	/*Go find in the binary file, the data we want to fetch and count the number of
	 * instances it appears: */
	if(my_rank==0){

		/*First set FILE* position to the beginning of the file: */
		fseek(fid,0,SEEK_SET);

		/*Now march through file looking for the correct data identifier: */
		for(;;){
			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				xDelete<char>(record_name);
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				_error_("error while looking in binary file. Found a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				break;
			}
			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}

			/*Is this the record sought for? : */
			if(strcmp(record_name,data_name)==0) num_instances++;

			/*Read the record length, and use it to skip the record: */
			if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
			fseek(fid,record_length,SEEK_CUR);
			xDelete<char>(record_name);
		}

		/*Ok, initialize the number of file handles we are going to return: */
		if(num_instances){
			file_positions = xNew<fpos_t>(num_instances);
			codes          = xNew<int>(num_instances);
			vector_types   = xNew<int>(num_instances);
		}

		/*Reset FILE* position to the beginning of the file, and start again, this time saving the data information
		 * as we find it: */
		counter=0;
		fseek(fid,0,SEEK_SET);

		for(;;){
			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				_error_("error while looking in binary file. Found a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				/*we have reached the end of the file. break: */
				xDelete<char>(record_name);
				break;
			}
			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}

			/*Is this the record sought for? : */
			if(strcmp(record_name,data_name)==0){
				/*Ok, we have found the correct string. Pass the record length, and read data type code: */
				fseek(fid,sizeof(long long),SEEK_CUR);
				if(fread(&record_code,sizeof(int),1,fid)!=1) _error_("Could not read record_code");

				/*if record_code points to a vector, get its type (nodal or elementary): */
				if(5<=record_code && record_code<=7){
					if(fread(&vector_type,sizeof(int),1,fid)!=1) _error_("Could not read vector_type");
				}
				codes[counter]        = record_code;
				vector_types[counter] = vector_type;
				fgetpos(fid,file_positions+counter);

				/*backup and skip over the record, as we have more work to do: */
				if(5<=record_code && record_code<=7) fseek(fid,-sizeof(int),SEEK_CUR); /*rewind for nodal or elementary type*/
				fseek(fid,-sizeof(int),SEEK_CUR);/*rewind for data code*/
				fseek(fid,-sizeof(long long),SEEK_CUR);/*rewind for record length*/

				/*increment counter: */
				counter++;
			}

			/*Read the record length, and use it to skip this record, as it has already been processed: */
			if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
			/*skip: */
			fseek(fid,record_length,SEEK_CUR);
			xDelete<char>(record_name);
		}
	}

	/*Broadcast data: */
	ISSM_MPI_Bcast(&num_instances,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Assign output pointers:*/
	*pcodes         = codes;
	*pnum_instances = num_instances;
	if(pvector_types){
		*pvector_types=vector_types;
	}
	else{
		xDelete<int>(vector_types);
	}
	return file_positions;
}
/*}}}*/
FILE* IoModel::SetFilePointerToData(int* pcode,int* pvector_type,const char* data_name){/*{{{*/

	int my_rank;

	int found  = 0;
	const char* mddot = "md.";
	char* record_name = NULL;
	int   record_name_size;
	long long record_length;
	int record_code;       //1 to 7 number
	int vector_type   = 0; //nodal or elementary

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();
	if(strncmp(data_name,mddot,3)!=0){
		_error_("Cannot fetch \""<<data_name<<"\" does not start with \""<<mddot<<"\"");
	}

	/*Go find in the binary file, the position of the data we want to fetch: */
	if(my_rank==0){

		/*First set FILE* position to the beginning of the file: */
		_assert_(fid);
		fseek(fid,0,SEEK_SET);

		/*Now march through file looking for the correct data identifier: */
		for(;;){
			/*Read size of first string name: */
			if(fread(&record_name_size,sizeof(int),1,fid)==0){
				/*we have reached the end of the file. break: */
				xDelete<char>(record_name);
				break;
			}
			if(record_name_size<3 || record_name_size>80){
				_error_("error while looking in binary file. Found a string of size "<<record_name_size);
			}

			/*Allocate string of correct size: */
			record_name=xNew<char>(record_name_size+1);
			record_name[record_name_size]='\0';

			/*Read record_name: */
			if(fread(record_name,record_name_size*sizeof(char),1,fid)==0){
				/*we have reached the end of the file. break: */
				found=0;
				xDelete<char>(record_name);
				break;
			}
			if(strncmp(record_name,mddot,3)!=0){
				_error_("error while reading binary file: record does not start with \"md.\": "<<record_name);
			}

			/*Is this the record sought for? : */
			if(strcmp(record_name,data_name)==0){
				/*Ok, we have found the correct string. Pass the record length, and read data type code: */
				fseek(fid,sizeof(long long),SEEK_CUR);
				if(fread(&record_code,sizeof(int),1,fid)!=1) _error_("Could not read record_code");
				/*if record_code points to a vector, get its type (nodal or elementary): */
				if((5<=record_code && record_code<=7) || record_code==10){
					if(fread(&vector_type,sizeof(int),1,fid)!=1) _error_("Could not read vector_type");
				}
				found=1;
				xDelete<char>(record_name);
				break;
			}
			else{
				/*This is not the correct string, read the record length, and use it to skip this record: */
				if(fread(&record_length,sizeof(long long),1,fid)!=1) _error_("Could not read record_length");
				/*skip: */
				fseek(fid,record_length,SEEK_CUR);
				xDelete<char>(record_name);
			}
		}
	}
	ISSM_MPI_Bcast(&found,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	if(!found) _error_("could not find data with name \"" << data_name << "\" in binary file");

	/*Broadcast code and vector type: */
	ISSM_MPI_Bcast(&record_code,1,ISSM_MPI_INT,0,IssmComm::GetComm());
	ISSM_MPI_Bcast(&vector_type,1,ISSM_MPI_INT,0,IssmComm::GetComm());

	/*Assign output pointers:*/
	*pcode=record_code;
	if(pvector_type)*pvector_type=vector_type;

	return fid;
}
/*}}}*/
void  IoModel::StartTrace(bool trace){/*{{{*/

	bool autodiff = false;
	bool iscontrol = false;
	bool keep=false;
	IssmDouble gcTriggerRatio;
	IssmDouble gcTriggerMaxSize;
	IssmDouble obufsize;
	IssmDouble lbufsize;
	IssmDouble cbufsize;
	IssmDouble tbufsize;

	int my_rank=IssmComm::GetRank();

	this->FetchData(&autodiff,"md.autodiff.isautodiff");
	this->FetchData(&iscontrol,"md.inversion.iscontrol");

	if(trace || (autodiff && !iscontrol)){

		#if defined(_HAVE_ADOLC_)
		/*Retrieve parameters: */
		this->FetchData(&keep,"md.autodiff.keep");
		int keepTaylors=keep?1:0;
		this->FetchData(&gcTriggerRatio,"md.autodiff.gcTriggerRatio");
		this->FetchData(&gcTriggerMaxSize,"md.autodiff.gcTriggerMaxSize");
		this->FetchData(&obufsize,"md.autodiff.obufsize");
		this->FetchData(&lbufsize,"md.autodiff.lbufsize");
		this->FetchData(&cbufsize,"md.autodiff.cbufsize");
		this->FetchData(&tbufsize,"md.autodiff.tbufsize");

		/*Set garbage collection parameters: */
		setStoreManagerControl(reCast<IssmPDouble>(gcTriggerRatio),reCast<size_t>(gcTriggerMaxSize));

		/*Start trace: */
		int skipFileDeletion=1;
		trace_on(my_rank,keepTaylors,reCast<size_t>(obufsize),reCast<size_t>(lbufsize),reCast<size_t>(cbufsize),reCast<size_t>(tbufsize),skipFileDeletion);

		#elif defined(_HAVE_CODIPACK_)
		//fprintf(stderr, "*** Codipack IoModel::StartTrace\n");
		codi_global.start();
		#endif
	}

}
/*}}}*/
