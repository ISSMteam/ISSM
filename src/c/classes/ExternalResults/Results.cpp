/*
 * \file Results.cpp
 * \brief: Implementation of the Results class, derived from DataSet class.
 */

/*Headers: {{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Results.h"
#include "./ExternalResult.h"
#include "../../shared/shared.h"
#include "../Params/Parameters.h"

using namespace std;
/*}}}*/

/*Object constructors and destructor*/
Results::Results(){/*{{{*/
	enum_type=ResultsEnum;
	return;
}
/*}}}*/
Results::~Results(){/*{{{*/
	return;
}
/*}}}*/

/*Object management*/
int Results::AddResult(ExternalResult* in_result){/*{{{*/

	/*First, go through dataset of inputs and check whether any input
	 * with the same name is already in. If so, erase the corresponding
	 * object before adding this new one: */

	/*In debugging mode, check that the input is not a NULL pointer*/
	_assert_(in_result);

	for(Object* &object : this->objects){
		ExternalResult* result=xDynamicCast<ExternalResult*>(object);

		if(result->GetStep()==in_result->GetStep()){
			char*    result_name =    result->GetResultName();
			char* in_result_name = in_result->GetResultName();
			if(strcmp(in_result_name,result_name)==0){

				this->DeleteObject(result);
				xDelete<char>(result_name);
				xDelete<char>(in_result_name);
				break;
			}
			xDelete<char>(result_name);
			xDelete<char>(in_result_name);
		}
	}
	this->AddObject(in_result);

	return 1;
}
/*}}}*/
int Results::DeleteResult(int result_enum,int result_step){/*{{{*/

	for(Object* &object : this->objects){
		ExternalResult* result=xDynamicCast<ExternalResult*>(object);
		if(result->GetStep()==result_step){
			if(strcmp(result->GetResultName(),EnumToStringx(result_enum))==0){
				this->DeleteObject(result);
				break;
			}
		}
	}

	return 1;
}
/*}}}*/
ExternalResult* Results::FindResult(int result_enum){/*{{{*/

	for(Object* &object : this->objects){
		ExternalResult* result=xDynamicCast<ExternalResult*>(object);
		if(result->GetResultEnum()==result_enum){
			return result;
		}
	}
	return NULL;
}
/*}}}*/
void Results::Write(Parameters* parameters){/*{{{*/

	FILE       *fid  = NULL;
	bool        io_gather;

	/*Recover file descriptor: */
	parameters->FindParam(&fid,OutputFilePointerEnum);
	parameters->FindParam(&io_gather,SettingsIoGatherEnum);

	for(Object* &object : this->objects){
		ExternalResult* result=xDynamicCast<ExternalResult*>(object);
		result->WriteData(fid,io_gather);
	}

}
/*}}}*/
