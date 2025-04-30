/*
 * \file Parameters.cpp
 * \brief: Implementation of the Parameters class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <vector>
#include <functional>
#include <algorithm>
#include <iostream>

#include "./Parameters.h"
#include "./Param.h"

#include "./BoolParam.h"
#include "./ControlParam.h"
#include "./DoubleMatParam.h"
#include "./DataSetParam.h"
#include "./DoubleParam.h"
#include "./DoubleVecParam.h"
#include "./IntParam.h"
#include "./IntVecParam.h"
#include "./IntMatParam.h"
#include "./FileParam.h"
#include "./MatrixParam.h"
#include "./VectorParam.h"
#include "./StringArrayParam.h"
#include "./StringParam.h"
#include "./DoubleMatArrayParam.h"
#include "./TransientParam.h"
#include "./TransientArrayParam.h"
#include "./TransientGriddedFieldParam.h"

#include "../../shared/shared.h"
#include "../../toolkits/toolkits.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Parameters::Parameters(){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++) this->params[i] = NULL;
	return;
}
/*}}}*/
Parameters::~Parameters(){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++){
		if(this->params[i]) delete this->params[i];
	}
	return;
}
/*}}}*/
int Parameters::EnumToIndex(int enum_in){/*{{{*/

	/*Make sure this parameter is at the right place*/
	#ifdef _ISSM_DEBUG_
	if(enum_in<=ParametersSTARTEnum) _error_("Enum "<<EnumToStringx(enum_in)<<" should appear after ParametersSTARTEnum");
	if(enum_in>=ParametersENDEnum)   _error_("Enum "<<EnumToStringx(enum_in)<<" should appear before ParametersENDEnum");
	#endif
	return enum_in - ParametersSTARTEnum -1;
}/*}}}*/

void Parameters::AddObject(Param* newparam){/*{{{*/

	/*Get Enum from Param*/
	_assert_(newparam);
	int param_enum = newparam->InstanceEnum();

	/*Get index in array*/
	int index = EnumToIndex(param_enum);

	/*Delete param if it already exists*/
	if(this->params[index]){
		delete this->params[index];
		this->params[index] = NULL;
	}

	/*Add param to array*/
	this->params[index] = newparam;
}
/*}}}*/
Parameters* Parameters::Copy(void){/*{{{*/

	Parameters* output = new Parameters();

	for(int i=0;i<NUMPARAMS;i++){
		if(this->params[i]){
			output->params[i]=this->params[i]->copy();
		}
	}

	return output;
}
/*}}}*/
void Parameters::DeepEcho(void){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++) {
		if(this->params[i]) this->params[i]->DeepEcho();
	}
	return;
}
/*}}}*/
void Parameters::Echo(void){/*{{{*/
	for(int i=0;i<NUMPARAMS;i++) {
		if(this->params[i]) this->params[i]->Echo();
	}
	return;
}
/*}}}*/
void Parameters::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int num_params=0;
	int obj_enum= ParametersEnum;
	marshallhandle->call(obj_enum);

	if(marshallhandle->OperationNumber()!=MARSHALLING_LOAD){

		/*Marshall num_params first*/
		for(int i=0;i<NUMPARAMS;i++){
			if(this->params[i]) num_params++;
		}
		marshallhandle->call(num_params);

		/*Marshall Parameters one by one now*/
		for(int i=0;i<NUMPARAMS;i++){
			if(this->params[i]){
				obj_enum = this->params[i]->ObjectEnum();
				marshallhandle->call(obj_enum);
				this->params[i]->Marshall(marshallhandle);
			}
		}

	}
	else{

		/*Get number of params marshalled*/
		marshallhandle->call(num_params);

		/*Recover parameters one by one*/
		for(int i=0;i<num_params;i++){

			/*Recover enum of object first: */
			marshallhandle->call(obj_enum);

			if(obj_enum==DoubleParamEnum){
				DoubleParam* doubleparam=NULL;
				doubleparam=new DoubleParam();
				doubleparam->Marshall(marshallhandle);
				this->AddObject(doubleparam);
			}
			else if(obj_enum==IntParamEnum){
				IntParam* intparam=NULL;
				intparam=new IntParam();
				intparam->Marshall(marshallhandle);
				this->AddObject(intparam);
			}
			else if(obj_enum==IntMatParamEnum){
				IntMatParam* intmparam=NULL;
				intmparam=new IntMatParam();
				intmparam->Marshall(marshallhandle);
				this->AddObject(intmparam);
			}
			else if(obj_enum==IntVecParamEnum){
				IntVecParam* intvparam=NULL;
				intvparam=new IntVecParam();
				intvparam->Marshall(marshallhandle);
				this->AddObject(intvparam);
			}
			else if(obj_enum==BoolParamEnum){
				BoolParam* boolparam=NULL;
				boolparam=new BoolParam();
				boolparam->Marshall(marshallhandle);
				this->AddObject(boolparam);
			}
			else if(obj_enum==DataSetParamEnum){
				DataSetParam* dsparam=NULL;
				dsparam=new DataSetParam();
				dsparam->Marshall(marshallhandle);
				this->AddObject(dsparam);
			}
			else if(obj_enum==DoubleMatArrayParamEnum){
				DoubleMatArrayParam* dmaparam=NULL;
				dmaparam=new DoubleMatArrayParam();
				dmaparam->Marshall(marshallhandle);
				this->AddObject(dmaparam);
			}
			else if(obj_enum==DoubleMatParamEnum){
				DoubleMatParam* dmparam=NULL;
				dmparam=new DoubleMatParam();
				dmparam->Marshall(marshallhandle);
				this->AddObject(dmparam);
			}
			else if(obj_enum==DoubleVecParamEnum){
				DoubleVecParam* dvparam=NULL;
				dvparam=new DoubleVecParam();
				dvparam->Marshall(marshallhandle);
				this->AddObject(dvparam);
			}
			else if(obj_enum==FileParamEnum){
				FileParam* fileparam=NULL;
				fileparam=new FileParam();
				fileparam->Marshall(marshallhandle);
				delete fileparam;
				/* FIXME: No need to add this object, the pointer is not valid
					The FemModel should reset all FileParams in the restart function */
			}
			else if(obj_enum==StringParamEnum){
				StringParam* sparam=NULL;
				sparam=new StringParam();
				sparam->Marshall(marshallhandle);
				this->AddObject(sparam);
			}
			else if(obj_enum==StringArrayParamEnum){
				StringArrayParam* saparam=NULL;
				saparam=new StringArrayParam();
				saparam->Marshall(marshallhandle);
				this->AddObject(saparam);
			}
			else if(obj_enum==TransientParamEnum){
				TransientParam* transparam=NULL;
				transparam=new TransientParam();
				transparam->Marshall(marshallhandle);
				this->AddObject(transparam);
			}
			else if(obj_enum==TransientArrayParamEnum){
				TransientArrayParam* transarrayparam=NULL;
				transarrayparam=new TransientArrayParam();
				transarrayparam->Marshall(marshallhandle);
				this->AddObject(transarrayparam);
			}
			else if(obj_enum==TransientGriddedFieldParamEnum){
				TransientGriddedFieldParam* transgridparam=NULL;
				transgridparam=new TransientGriddedFieldParam();
				transgridparam->Marshall(marshallhandle);
				this->AddObject(transgridparam);
			}
			else if(obj_enum==ControlParamEnum){
				ControlParam* controlparam=NULL;
				controlparam=new ControlParam();
				controlparam->Marshall(marshallhandle);
				this->AddObject(controlparam);
			}
			else if(obj_enum==GenericParamEnum){
				/*Skip for now (we don't want to Marhsall Comms)*/
			}
		}
	}
}
/*}}}*/

/*Object management*/
void Parameters::Delete(int param_enum){/*{{{*/

	int index = EnumToIndex(param_enum);
	if(this->params[index]){
		delete this->params[index];
		this->params[index] = NULL;
	}

	return;
}
/*}}}*/
bool Parameters::Exist(int param_enum){/*{{{*/

	int index = EnumToIndex(param_enum);
	if(this->params[index]) return true;

	return false;
}
/*}}}*/
void Parameters::FindParam(bool* pbool,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pbool);
}
/*}}}*/
void Parameters::FindParam(int* pinteger,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pinteger);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int param_enum){ _assert_(this);/*{{{*/
	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int param_enum,IssmDouble time){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,time);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int param_enum,IssmDouble time,int timestepping,IssmDouble dt){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,time,timestepping,dt);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int row,IssmDouble time,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,row,time);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int row,int column,IssmDouble time,int param_enum){ _assert_(this);/*{{{*/
	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,row,column,time);
}
/*}}}*/
void Parameters::FindParam(IssmDouble* pscalar,int row,IssmDouble time,int timestepping,IssmDouble dt, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pscalar,row,time,timestepping,dt);
}
/*}}}*/
void Parameters::FindParam(char** pstring,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pstring);

}
/*}}}*/
void Parameters::FindParam(char*** pstringarray,int* pM,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pstringarray,pM);
}
/*}}}*/
void Parameters::FindParam(int** pintarray,int* pM, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pintarray,pM);

}
/*}}}*/
void Parameters::FindParam(int** pintarray,int* pM,int *pN,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pintarray,pM,pN);

}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM);

}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int* pN, IssmDouble time, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM,pN,time);

}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int* pN, IssmDouble starttime, IssmDouble endtime, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM,pN,starttime,endtime);

}
/*}}}*/
void Parameters::FindParam(IssmDouble** pIssmDoublearray,int* pM, int* pN,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pIssmDoublearray,pM,pN);
}
/*}}}*/
void Parameters::FindParam(IssmDouble*** parray,int* pM,int** pmdims_array,int** pndims_array,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(parray,pM,pmdims_array,pndims_array);
}
/*}}}*/
void Parameters::FindParam(Vector<IssmDouble>** pvec,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pvec);
}
/*}}}*/
void Parameters::FindParam(Matrix<IssmDouble>** pmat,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pmat);
}
/*}}}*/
void Parameters::FindParam(FILE** pfid,int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pfid);
}
/*}}}*/
void Parameters::FindParam(DataSet** pdataset,int param_enum){ /*{{{*/
	_assert_(this);

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pdataset);
}
/*}}}*/
void Parameters::FindParamAndMakePassive(IssmPDouble* pscalar,int param_enum){ _assert_(this);/*{{{*/

	/*Get "active" parameter*/
	IssmDouble intermediary;
	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(&intermediary);

	/*cast to "passive"*/
	#ifdef _HAVE_AD_
	*pscalar=reCast<IssmPDouble>(intermediary);
	#else
	*pscalar=intermediary;
	#endif
}
/*}}}*/
void Parameters::FindParamAndMakePassive(IssmPDouble** pvec,int* pM, int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	/*Output*/
	int         n;
	IssmDouble* vector = NULL;

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(&vector,&n);

	/*Make output passive*/
	#ifdef _HAVE_AD_
	IssmPDouble* output = xNew<IssmPDouble>(n);
	for(int i=0;i<n;i++) output[i] = reCast<IssmPDouble>(vector[i]);
	xDelete<IssmDouble>(vector);
	if(pvec) *pvec = output;
	#else
	if(pvec) *pvec = vector;
	#endif

	/*assign output pointers*/
	if(pM)   *pM   = n;
}/*}}}*/
void Parameters::FindControlParam(IssmDouble** pvec,int* pM, int param_enum, const char* data){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	/*Output*/
	int         n;
	IssmDouble* vector = NULL;

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(pvec,pM,data);

}/*}}}*/
void Parameters::FindControlParamAndMakePassive(IssmPDouble** pvec,int* pM, int param_enum, const char* data){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);

	/*Output*/
	int         n;
	IssmDouble* vector = NULL;

	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");
	this->params[index]->GetParameterValue(&vector,&n,data);

	/*Make output passive*/
	#ifdef _HAVE_AD_
	IssmPDouble* output = xNew<IssmPDouble>(n);
	for(int i=0;i<n;i++) output[i] = reCast<IssmPDouble>(vector[i]);
	xDelete<IssmDouble>(vector);
	if(pvec) *pvec = output;
	#else
	if(pvec) *pvec = vector;
	#endif

	/*assign output pointers*/
	if(pM)   *pM   = n;
}/*}}}*/
IssmDouble Parameters::FindParam(int param_enum){ _assert_(this);/*{{{*/

	int index = EnumToIndex(param_enum);
	if(!this->params[index]) _error_("Parameter " << EnumToStringx(param_enum) <<" not set");

	IssmDouble value;
	this->params[index]->GetParameterValue(&value);
	return value;
}
/*}}}*/

void   Parameters::SetParam(bool boolean,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(boolean); //already exists, just set it.
	else this->AddObject(new BoolParam(enum_type,boolean)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int integer,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(integer); //already exists, just set it.
	else this->AddObject(new IntParam(enum_type,integer)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble scalar,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));
	if(param) param->SetValue(scalar); //already exists, just set it.
	else this->AddObject(new DoubleParam(enum_type,scalar)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(char* string,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(string); //already exists, just set it.
	else this->AddObject(new StringParam(enum_type,string)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(char** stringarray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(stringarray,M); //already exists, just set it.
	else this->AddObject(new StringArrayParam(enum_type,stringarray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble* IssmDoublearray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(IssmDoublearray,M); //already exists, just set it.
	else this->AddObject(new DoubleVecParam(enum_type,IssmDoublearray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble* IssmDoublearray,int M, int N, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(IssmDoublearray,M,N); //already exists, just set it.
	else this->AddObject(new DoubleMatParam(enum_type,IssmDoublearray,M,N)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(IssmDouble* IssmDoublearray, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));
	if(param) param->SetValue(IssmDoublearray); //already exists, just set it.
	else _error_("Param "<< EnumToStringx(enum_type) << " cannot setValue");

	 //this->AddObject(new ControlParam(enum_type,IssmDoublearray,M,N)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int* intarray,int M, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(intarray,M); //already exists, just set it.
	else this->AddObject(new IntVecParam(enum_type,intarray,M)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(int* intarray,int M, int N, int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(intarray,M,N); //already exists, just set it.
	else this->AddObject(new IntMatParam(enum_type,intarray,M,N)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(Vector<IssmDouble>* vector,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(vector); //already exists, just set it.
	else this->AddObject(new VectorParam(enum_type,vector)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(Matrix<IssmDouble>* matrix,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(matrix); //already exists, just set it.
	else this->AddObject(new MatrixParam(enum_type,matrix)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(FILE* fid,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(fid); //already exists, just set it.
	else this->AddObject(new FileParam(enum_type,fid)); //just add the new parameter.
}
/*}}}*/
void   Parameters::SetParam(DataSet* dataset,int enum_type){/*{{{*/

	Param* param=NULL;

	/*first, figure out if the param has already been created: */
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param){
		param->SetValue(dataset); //already exists, just set it.
	}
	else{
		this->AddObject(new DataSetParam(enum_type,dataset)); //just add the new parameter.
	}
}
/*}}}*/
void   Parameters::SetControlFromVector(IssmDouble* vector, int enum_type, int M, int N, int offset){/*{{{*/

	/*first, figure out if the param has already been created: */
	Param* param=NULL;
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetValue(&vector[offset], M, N);
	else _error_("Param "<< EnumToStringx(enum_type) << " cannot setValue");
}
/*}}}*/
void   Parameters::SetGradientFromVector(IssmDouble* vector, int enum_type, int M, int N, int offset){/*{{{*/

	/*first, figure out if the param has already been created: */
	Param* param=NULL;
	param=xDynamicCast<Param*>(this->FindParamObject(enum_type));

	if(param) param->SetGradient(&vector[offset], M, N);
	else _error_("Param "<< EnumToStringx(enum_type) << " cannot setValue");
}
/*}}}*/

void  Parameters::GetVectorFromControl(Vector<IssmDouble>* vector,int control_enum,int control_index,int N,const char* data,int offset){/*{{{*/

	/*first, figure out if the param has already been created: */
	Param* param=xDynamicCast<Param*>(this->FindParamObject(control_enum));
	if(!param) _error_("Parameter not found");

	param->GetVectorFromControl(vector, control_index, N, data, offset);
}/*}}}*/

Param* Parameters::FindParamObject(int param_enum){/*{{{*/

	return this->params[EnumToIndex(param_enum)];
}
/*}}}*/

/*Methods relating to parameters: */
char* OptionsFromAnalysis(char** pouttoolkit,Parameters* parameters,int analysis_type){ /*{{{*/

	/* figure out ISSM options for current analysis, return a string. */

	/*output: */
	char *outstring  = NULL;
	char *outtoolkit = NULL;

	/*intermediary: */
	int          dummy;
	int         *analyses    = NULL;
	char       **strings     = NULL;
	char        *string      = NULL;
	char       **toolkits    = NULL;
	char        *toolkit     = NULL;
	int          numanalyses;
	int          found       = -1;
	int          i;

	parameters->FindParam(&strings,&numanalyses,ToolkitsOptionsStringsEnum);
	parameters->FindParam(&toolkits,&dummy,ToolkitsTypesEnum); _assert_(dummy==numanalyses);
	parameters->FindParam(&analyses,&dummy,ToolkitsOptionsAnalysesEnum); _assert_(dummy==numanalyses);

	if(numanalyses==0)return NULL; //we did not find petsc options, don't bother.

	/*ok, go through analyses and figure out if it corresponds to our analysis_type: */
	for(i=0;i<numanalyses;i++){
		if(analyses[i]==analysis_type){
			found=i;
			break;
		}
	}
	if(found==-1){
		/*still haven't found a list of petsc options, go find the default one, for analysis type DefaultAnalysisEnum: */
		for(i=0;i<numanalyses;i++){
			if(analyses[i]==DefaultAnalysisEnum){
				found=i;
				break;
			}
		}
	}
	if(found==-1){
		_error_("could find neither a default analysis nor analysis " << EnumToStringx(analysis_type));
	}

	/*1. Grab the option toolkit: */
	outtoolkit=xNew<char>(strlen(toolkits[found])+1);
	strcpy(outtoolkit,toolkits[found]);
	*pouttoolkit = outtoolkit;

	/*2. Grab the option string: */
	outstring=xNew<char>(strlen(strings[found])+1);
	strcpy(outstring,strings[found]);

	/*Free resources:*/
	for(i=0;i<numanalyses;i++){
		xDelete<char>(toolkits[i]);
		xDelete<char>(strings[i]);
	}
	xDelete<char*>(toolkits);
	xDelete<char*>(strings);
	xDelete<int>(analyses);
	return outstring;
}
/*}}}*/
void ToolkitsOptionsFromAnalysis(Parameters* parameters,int analysis_type){ /*{{{*/

	/*!\file:  ToolkitsOptionsFromAnalysis.cpp
	 * \brief: for each analysis, setup the issmoptions string.
	 * This is mainly for the case where we run our toolkits using petsc. In this case, we need to
	 * plug our toolkits options directly into the petsc options database. This is the case for each analysis type
	 * and parameters
	 */

	char* options = NULL;
	char* toolkit = NULL;

	/*Recover first the options string for this analysis: */
	options=OptionsFromAnalysis(&toolkit,parameters,analysis_type);

	/*Initialize our Toolkit Options: */
	ToolkitOptions::Init(toolkit,options);

	#ifdef _HAVE_PETSC_
		/* In case we are using PETSC, we do not rely on issmoptions. Instead, 
		   we dump issmoptions into the Petsc options database */

		#if (_PETSC_MINOR_>=7 && _PETSC_MINOR_<14)
		PetscOptionsSetFromOptions(NULL);
		#endif
		#if (_PETSC_MINOR_>=7)
		PetscOptionsClear(NULL);
		int ierr = PetscOptionsInsertString(NULL,options);
		//int ierr = PetscOptionsInsertString(NULL,"-mat_type mpiaij -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 120 -mat_mumps_icntl_28 2 -mat_mumps_icntl_29 2");
		//int ierr = PetscOptionsInsertString(NULL,"-mat_type mpiaij -ksp_type preonly -pc_type lu -pc_factor_mat_solver_package mumps -mat_mumps_icntl_14 120");
		#else
		PetscOptionsSetFromOptions();
		PetscOptionsClear();
		int ierr = PetscOptionsInsertString(options);
		#endif

		if(ierr) _error_("Could not enter PETSc options");

	#endif

	xDelete<char>(options);
	xDelete<char>(toolkit);
}
/*}}}*/
