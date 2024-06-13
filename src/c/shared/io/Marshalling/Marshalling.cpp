/*!\file:  Marshalling.cpp
 * \brief implement marshall
 */ 

#ifdef HAVE_CONFIG_H
#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./Marshalling.h"
#include "../../Numerics/isnan.h"

WriteCheckpointFunctor::WriteCheckpointFunctor(char** pmarshalled_data_in) : MarshallHandle(MARSHALLING_WRITE){/*{{{*/
	this->pmarshalled_data = pmarshalled_data_in;
}/*}}}*/
void WriteCheckpointFunctor::Echo(void){/*{{{*/
	printf("WriteCheckpointFunctor Echo:\n");
	printf("   pmarshalled_data: %p\n",pmarshalled_data);
}/*}}}*/
void WriteCheckpointFunctor::call(char* & value){/*{{{*/
	int size = strlen(value)+1;
	this->call(size);
	this->call(value,size);
};/*}}}*/

LoadCheckpointFunctor::LoadCheckpointFunctor(char** pmarshalled_data_in) : MarshallHandle(MARSHALLING_LOAD){/*{{{*/
	this->pmarshalled_data = pmarshalled_data_in;
}/*}}}*/
void LoadCheckpointFunctor::Echo(void){/*{{{*/
	printf("LoadCheckpointFunctor Echo:\n");
	printf("   pmarshalled_data: %p\n",pmarshalled_data);
}/*}}}*/
void LoadCheckpointFunctor::call(char* & value){/*{{{*/
	int size;
	this->call(size);
	this->call(value,size);
};/*}}}*/

SizeCheckpointFunctor::SizeCheckpointFunctor(void) : MarshallHandle(MARSHALLING_SIZE){/*{{{*/
	this->marshalled_data_size = 0;
}/*}}}*/
int  SizeCheckpointFunctor::MarshalledSize(void){/*{{{*/
	return this->marshalled_data_size;
};/*}}}*/
void SizeCheckpointFunctor::Echo(void){/*{{{*/
	printf("SizeCheckpointFunctor Echo:\n");
	printf("   marshalled_data_size: %i\n",marshalled_data_size);
}/*}}}*/
void SizeCheckpointFunctor::call(char* & value){/*{{{*/
	int size = strlen(value)+1;
	this->call(size);
	this->call(value,size);
};/*}}}*/

#if defined(_HAVE_CODIPACK_) && !defined(_WRAPPERS_)
CountDoublesFunctor::CountDoublesFunctor(void) : MarshallHandle(AD_COUNTDOUBLES){/*{{{*/
	this->double_count= 0;
}/*}}}*/
int  CountDoublesFunctor::DoubleCount(void){/*{{{*/
	return this->double_count;
};/*}}}*/
void CountDoublesFunctor::Echo(void){/*{{{*/
	printf("CountDoublesFunctor Echo:\n");
	printf("   double_count: %i\n",double_count);
}/*}}}*/
void CountDoublesFunctor::call(IssmDouble & value){/*{{{*/
	this->double_count++;
}/*}}}*/
void CountDoublesFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value) this->double_count+= size;
}/*}}}*/

RegisterInputFunctor::RegisterInputFunctor(CoDi_global *data) : MarshallHandle(AD_REGISTERINPUT){/*{{{*/
	this->data = data;

}/*}}}*/
void RegisterInputFunctor::Echo(void){/*{{{*/
	printf("RegisterInputFunctor Echo:\n");
	printf("   double_count: %i\n",(int)data->input_indices.size());
}/*}}}*/
void RegisterInputFunctor::call(IssmDouble & value){/*{{{*/
	/*Comment out this assert, some parameters are NaN (e.g. abstol) by default*/
	//_assert_(!xIsNan<IssmDouble>(value));

	this->data->registerInput(value);
}/*}}}*/
void RegisterInputFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value){
		for(int i=0;i<size;i++){
			_assert_(!xIsNan<IssmDouble>(value[i]));
			this->data->registerInput(value[i]);
		}
	}
}/*}}}*/

RegisterOutputFunctor::RegisterOutputFunctor(CoDi_global *data) : MarshallHandle(AD_REGISTEROUTPUT){/*{{{*/
	this->data = data;
}/*}}}*/
void RegisterOutputFunctor::Echo(void){/*{{{*/
	printf("RegisterOutputFunctor Echo:\n");
	printf("   double_count: %i\n", (int)data->output_indices.size());
}/*}}}*/
void RegisterOutputFunctor::call(IssmDouble & value){/*{{{*/
	//_assert_(!xIsNan<IssmDouble>(value));
	this->data->registerOutput(value);
}/*}}}*/
void RegisterOutputFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value){
		for(int i=0;i<size;i++){
			_assert_(!xIsNan<IssmDouble>(value[i]));
			this->data->registerOutput(value[i]);
		}
	}
}/*}}}*/
#endif
