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

RegisterInputFunctor::RegisterInputFunctor(int* identifiers_in,int size_max_in) : MarshallHandle(AD_REGISTERINPUT){/*{{{*/
	this->double_count = 0;
	this->identifiers  = identifiers_in;
	this->size_max         = size_max_in;
	#if _CODIPACK_MAJOR_==2
	this->tape_codi    = &(IssmDouble::getTape());
	#elif _CODIPACK_MAJOR_==1
	this->tape_codi    = &(IssmDouble::getGlobalTape());
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif

}/*}}}*/
void RegisterInputFunctor::Echo(void){/*{{{*/
	printf("RegisterInputFunctor Echo:\n");
	printf("   double_count: %i\n",double_count);
}/*}}}*/
void RegisterInputFunctor::call(IssmDouble & value){/*{{{*/
	_assert_(this->double_count<size_max);

	/*Comment out this assert, some parameters are NaN (e.g. abstol) by default*/
	//_assert_(!xIsNan<IssmDouble>(value));

	this->tape_codi->registerInput(value);
	#if _CODIPACK_MAJOR_==2
	this->identifiers[this->double_count] = value.getIdentifier();
	#elif _CODIPACK_MAJOR_==1
	this->identifiers[this->double_count] = value.getGradientData();
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif

	this->double_count++;
}/*}}}*/
void RegisterInputFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value){
		for(int i=0;i<size;i++){
			_assert_(this->double_count<size_max);
			_assert_(!xIsNan<IssmDouble>(value[i]));
			this->tape_codi->registerInput(value[i]);
			#if _CODIPACK_MAJOR_==2
			this->identifiers[this->double_count] = value[i].getIdentifier();
			#elif _CODIPACK_MAJOR_==1
			this->identifiers[this->double_count] = value[i].getGradientData();
			#else
			#error "_CODIPACK_MAJOR_ not supported"
			#endif

			this->double_count++;
		}
	}
}/*}}}*/

RegisterOutputFunctor::RegisterOutputFunctor(void) : MarshallHandle(AD_REGISTEROUTPUT){/*{{{*/
	this->double_count = 0;
	#if _CODIPACK_MAJOR_==2
	this->tape_codi    = &(IssmDouble::getTape());
	#elif _CODIPACK_MAJOR_==1
	this->tape_codi    = &(IssmDouble::getGlobalTape());
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif
}/*}}}*/
void RegisterOutputFunctor::Echo(void){/*{{{*/
	printf("RegisterOutputFunctor Echo:\n");
	printf("   double_count: %i\n",double_count);
}/*}}}*/
void RegisterOutputFunctor::call(IssmDouble & value){/*{{{*/
	//_assert_(!xIsNan<IssmDouble>(value));
	this->tape_codi->registerOutput(value);
	this->double_count++;
}/*}}}*/
void RegisterOutputFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value){
		for(int i=0;i<size;i++){
			_assert_(!xIsNan<IssmDouble>(value[i]));
			this->tape_codi->registerOutput(value[i]);
			this->double_count++;
		}
	}
}/*}}}*/

SetAdjointFunctor::SetAdjointFunctor(double* adjoint_in,int size_max_in) : MarshallHandle(AD_SETADJOINT){/*{{{*/
	this->double_count = 0;
	this->adjoint      = adjoint_in;
	this->size_max     = size_max_in;
	#if _CODIPACK_MAJOR_==2
	this->tape_codi    = &(IssmDouble::getTape());
	#elif _CODIPACK_MAJOR_==1
	this->tape_codi    = &(IssmDouble::getGlobalTape());
	#else
	#error "_CODIPACK_MAJOR_ not supported"
	#endif
}/*}}}*/
void SetAdjointFunctor::Echo(void){/*{{{*/
	printf("SetAdjointFunctor Echo:\n");
	printf("   double_count: %i\n",double_count);
}/*}}}*/
void SetAdjointFunctor::call(IssmDouble & value){/*{{{*/
	_assert_(this->double_count<size_max);
	_assert_(!xIsNan<IssmDouble>(this->adjoint[this->double_count]));
	value.gradient() = this->adjoint[this->double_count];
	this->double_count++;
}/*}}}*/
void SetAdjointFunctor::call(IssmDouble* & value,int size){/*{{{*/
	if(value){
		for(int i=0;i<size;i++){
			_assert_(this->double_count<size_max);
			_assert_(!xIsNan<IssmDouble>(this->adjoint[this->double_count]));
			value[i].gradient() = this->adjoint[this->double_count];
			this->double_count++;
		}
	}
}/*}}}*/
#endif
