/*! \file GenericExternalResult.h 
 *  \brief: header file for generic external result object
 */

#ifndef _GENERIC_EXTERNAL_RESULT_
#define _GENERIC_EXTERNAL_RESULT_

/*Headers:{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cstring>
#include "./ExternalResult.h"
#include "../../shared/shared.h"
/*}}}*/

template <class ResultType> 
class GenericExternalResult: public ExternalResult {

	private: 
		int        id;
		char*      result_name;
		ResultType value;
		int        M;
		int        N;
		int        step;
		IssmDouble time;

	public:
		/*Diverse: must be in front, as it is used in what follows*/
		void GenericEcho(void){/*{{{*/
			_printf_("   id          : " << this->id << "\n");
			_printf_("   result_name : " << this->result_name<< "\n");
			_printf_("   step        : " << this->step << "\n");
			_printf_("   time        : " << this->time << "\n");
		}
		/*}}}*/
		void GenericWriteData(FILE* fid){/*{{{*/ 

			IssmPDouble  passiveDouble;

			/*First write name: */
			int length=(strlen(this->result_name)+1)*sizeof(char);
			fwrite(&length,sizeof(int),1,fid);
			fwrite(this->result_name,length,1,fid);

			/*Now write time and step: */
			passiveDouble=reCast<IssmPDouble>(time);
			fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
			fwrite(&step,sizeof(int),1,fid);
		} /*}}}*/
		void GenericMarshall(MarshallHandle* marshallhandle){/*{{{*/

			int object_enum = this->ObjectEnum();
			marshallhandle->call(object_enum);
			marshallhandle->call(this->id);
			marshallhandle->call(this->result_name);
			marshallhandle->call(this->M);
			marshallhandle->call(this->N);
			marshallhandle->call(this->step);
			marshallhandle->call(this->time);
			marshallhandle->call(this->value);
		}  /*}}}*/

		/*GenericExternalResult constructors and  destructors*/
		GenericExternalResult(){ /*{{{*/
			id          = 0;
			result_name = NULL;
			M           = 0;
			N           = 0;
			step        = 0;
			time        = 0;
			value       = 0;
		} /*}}}*/
		GenericExternalResult(int in_id, int in_enum_type,ResultType in_values, int in_M,int in_N,int in_step=UNDEF,IssmDouble in_time=UNDEF){/*{{{*/
			id        = in_id;
			step      = in_step;
			time      = in_time;
			M         = in_M;
			N         = in_N;

			/*Copy result in values*/
			if(M*N){
				value=xNew<IssmDouble>(M*N);
				xMemCpy<IssmDouble>(value,in_values,M*N);
			}
			else value=NULL;

			/*Convert enum to name*/
			EnumToStringx(&this->result_name,in_enum_type);
		}
/*}}}*/
		GenericExternalResult(int in_id,const char* name_in,ResultType in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/
			_error_("template GenericExternalResult(int in_id, int in_enum_type,double* in_values, int in_M,int in_N,int in_step,IssmDouble in_time) not implemented for this ResultType\n");
		}
		/*}}}*/
		GenericExternalResult(int in_id, int in_enum_type,ResultType in_value,int in_step=UNDEF, IssmDouble in_time=UNDEF){ /*{{{*/
			id        = in_id;
			value     = in_value;
			step      = in_step;
			time      = in_time;
			M         = 1;
			N         = 1;

			/*Convert enum to name*/
			EnumToStringx(&this->result_name,in_enum_type);
		}
		/*}}}*/
		GenericExternalResult(int in_id,const char* in_result_name,ResultType in_value,int in_step=UNDEF, IssmDouble in_time=UNDEF){ /*{{{*/
			id    = in_id;
			value = in_value;
			step  = in_step;
			time  = in_time;
			M     = 1;
			N     = 1;

			/*Copy name*/
			this->result_name = xNew<char>(strlen(in_result_name)+1);
			xMemCpy<char>(this->result_name,in_result_name,strlen(in_result_name)+1);
		}
		/*}}}*/
		~GenericExternalResult(){ /*{{{*/
			xDelete<char>(result_name);
		} /*}}}*/

		/*Object virtual functions definitions:*/
		Object* copy(void) { /*{{{*/
			return new GenericExternalResult<ResultType>(this->id,this->result_name,this->value,this->step,this->time);
		} /*}}}*/
		void Echo(void){ /*{{{*/
			this->DeepEcho();
		}
		/*}}}*/
		void DeepEcho(void){ /*{{{*/
			_error_("template DeepEcho not implemented for this ResultType\n");
		}
		/*}}}*/
		int Id(void){ /*{{{*/ 
			return -1; 
		} /*}}}*/
		int ObjectEnum(void){ /*{{{*/
			return GenericExternalResultEnum;
		} /*}}}*/
		void Marshall(MarshallHandle* marshallhandle){/*{{{*/
			_printf_("   WARNING: result "<<this->result_name<<" is a GenericExternalResult and cannot be marshalled (need overload)\n");
			/*Nothing for now*/
		} 
		/*}}}*/

		/*GenericExternalResult management: */
void  WriteData(FILE* fid,bool io_gather){ /*{{{*/

	int     my_rank;
	int     type;
	int     size;
	IssmPDouble  passiveDouble;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*return if now on cpu 0: */
	if(my_rank)return;

	/*use generic part, same for all ResultTypes: */
	this->GenericWriteData(fid);

	/*writing a IssmPDouble for Matlab or Python to post-process, type is 1, size is 1: */
	type=1;
	size=1;
	fwrite(&type,sizeof(int),1,fid);
	fwrite(&size,sizeof(int),1,fid);

	/*cast to a IssmPDouble: */
	passiveDouble=reCast<IssmPDouble>(value);
	fwrite(&passiveDouble,size*sizeof(IssmPDouble),1,fid);

} /*}}}*/
void  Transpose(void){ /*{{{*/
	_error_("not implemented yet");
} /*}}}*/
char* GetResultName(void){ /*{{{*/
		char* name = xNew<char>(strlen(this->result_name)+1);
		xMemCpy<char>(name,this->result_name,strlen(this->result_name)+1);
		return name;
} /*}}}*/
int GetResultEnum(void){ /*{{{*/
		return StringToEnumx(this->result_name,false);
} /*}}}*/
int   GetStep(void){ /*{{{*/
	return this->step;
} /*}}}*/
double GetValue(void){ /*{{{*/
	/*Only supported by IssmPDouble result, error out by default*/
	_error_("not supported for this type of result");
} /*}}}*/
double* GetValues(void){ /*{{{*/
	/*Only supported by IssmPDouble* result, error out by default*/
	_error_("not supported for this type of result");
} /*}}}*/
};

/*Specific instantiations for bool: */
template <> inline void GenericExternalResult<bool>::DeepEcho(void){ /*{{{*/

	_printf_("GenericExternalResult<bool>:\n");
	this->GenericEcho();
	_printf_("   value: " <<(this->value?"true":"false") << "\n");

} /*}}}*/
template <> inline int GenericExternalResult<bool>::ObjectEnum(void){ /*{{{*/
	return BoolExternalResultEnum;
} /*}}}*/
template <> inline void GenericExternalResult<bool>::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	this->GenericMarshall(marshallhandle);

}  /*}}}*/

/*Specific instantiations for int: */
template <> inline void GenericExternalResult<int>::DeepEcho(void){ /*{{{*/

	_printf_("GenericExternalResult<int>:\n");
	this->GenericEcho();
	_printf_("   value: " << this->value << "\n");

} /*}}}*/
template <> inline int GenericExternalResult<int>::ObjectEnum(void){ /*{{{*/
	return IntExternalResultEnum;
} /*}}}*/
template <> inline void GenericExternalResult<int>::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	this->GenericMarshall(marshallhandle);
}  /*}}}*/

/*Specific instantiations for double: */
template <> inline void GenericExternalResult<double>::DeepEcho(void){ /*{{{*/

	_printf_("GenericExternalResult<double>:\n");
	this->GenericEcho();
	_printf_("   value: " << this->value << "\n");

} /*}}}*/
template <> inline int GenericExternalResult<double>::ObjectEnum(void){ /*{{{*/
	return DoubleExternalResultEnum;
} /*}}}*/
template <> inline double GenericExternalResult<double>::GetValue(void){ /*{{{*/
	return value;
} /*}}}*/
template <> inline void GenericExternalResult<double>::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	this->GenericMarshall(marshallhandle);
}  /*}}}*/

/*Specific instantiations for char*: */
template <> inline GenericExternalResult<char*>::GenericExternalResult(int in_id, int in_enum_type,char* in_value,int in_step, IssmDouble in_time){ /*{{{*/

	id   = in_id;
	step = in_step;
	time = in_time;
	M    = 1;
	N    = 1;

	value = xNew<char>(strlen(in_value)+1);
	xMemCpy<char>(value,in_value,(strlen(in_value)+1));

	/*Convert enum to name*/
	EnumToStringx(&this->result_name,in_enum_type);

} /*}}}*/
template <> inline GenericExternalResult<char*>::~GenericExternalResult(){ /*{{{*/
	xDelete<char>(result_name);
	xDelete<char>(value);
} /*}}}*/
template <> inline void GenericExternalResult<char*>::DeepEcho(void){ /*{{{*/

	_printf_("GenericExternalResult<char*>:\n");
	this->GenericEcho();
	_printf_("   value: " << this->value << "\n");

} /*}}}*/
template <> inline void GenericExternalResult<char*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

	int     my_rank;
	int     type;
	int     length;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	/*return if now on cpu 0: */
	if(my_rank)return;

	/*use generic part, same for all ResultTypes: */
	this->GenericWriteData(fid);

	/*writing a string, type is 2: */
	type=2;
	fwrite(&type,sizeof(int),1,fid);

	length=(strlen(this->value)+1)*sizeof(char);
	fwrite(&length,sizeof(int),1,fid);
	fwrite(this->value,length,1,fid);
}
/*}}}*/
template <> inline int GenericExternalResult<char*>::ObjectEnum(void){ /*{{{*/
	return StringExternalResultEnum;
} /*}}}*/
template <> inline void GenericExternalResult<char*>::Marshall(MarshallHandle* marshallhandle){/*{{{*/
	marshallhandle->call(this->id);
	marshallhandle->call(this->result_name);
	marshallhandle->call(this->value);
	marshallhandle->call(this->step);
	marshallhandle->call(this->time);

}  /*}}}*/

/*Specific instantiations for int*: */
template <> inline GenericExternalResult<int*>::GenericExternalResult(int in_id, int in_enum_type,int* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/

	id = in_id;
	M  = in_M;
	N  = in_N;

	EnumToStringx(&this->result_name,in_enum_type);

	step = in_step;
	time = in_time;

	/*Copy result in values*/
	if(M*N){
		value=xNew<int>(M*N);
		xMemCpy<int>(value,in_values,M*N);
	}
	else value=NULL;
}
/*}}}*/
template <> inline GenericExternalResult<int*>::GenericExternalResult(int in_id,const char* in_result_name,int* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/

	id = in_id;
	M  = in_M;
	N  = in_N;

	/*Copy name*/
	this->result_name = xNew<char>(strlen(in_result_name)+1);
	xMemCpy<char>(this->result_name,in_result_name,strlen(in_result_name)+1);

	step = in_step;
	time = in_time;

	/*Copy result in values*/
	if(M*N){
		value=xNew<int>(M*N);
		xMemCpy<int>(value,in_values,M*N);
	}
	else value=NULL;
}
/*}}}*/
template <> inline GenericExternalResult<int*>::GenericExternalResult(int in_id, int in_enum_type,int* in_value,int in_step, IssmDouble in_time){ /*{{{*/
	_error_("you cannot initialize a GenericExternalResult<int*> without providing the dimensions of the matrix! Please use a more appropriate constructor!");
} /*}}}*/
template <> inline GenericExternalResult<int*>::~GenericExternalResult(){ /*{{{*/
	xDelete<char>(result_name);
	xDelete<int>(value);
} /*}}}*/
template <> inline void GenericExternalResult<int*>::Echo(void){ /*{{{*/

	_printf_("GenericExternalResult<int*>:\n");
	this->GenericEcho();
	_printf_("   matrix size: " << this->M << "-" << this->N << "\n");

} /*}}}*/
template <> inline void GenericExternalResult<int*>::DeepEcho(void){ /*{{{*/

	int i,j;

	_printf_("GenericExternalResult<int*>:\n");
	this->GenericEcho();

	_printf_("   matrix size: " << this->M << "-" << this->N << "\n");
	for (i=0;i<this->M;i++){  
		_printf_("   [ ");
		for (j=0;j<this->N;j++){
			_printf_( " " << setw(11) << this->value[i*this->N+j]);
		}  
		_printf_(" ]\n");
	}  

} /*}}}*/
template <> inline Object* GenericExternalResult<int*>::copy(void){ /*{{{*/
	return new GenericExternalResult<int*>(this->id,StringToEnumx(this->result_name),this->value,this->M,this->N,this->step,this->time);
} /*}}}*/
template <> inline void GenericExternalResult<int*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

	int     my_rank;
	int     type;
	int     rows,cols;
	char   *name    = NULL;
	IssmPDouble passiveDouble;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if(io_gather){
		/*we are gathering the data on cpu 0, don't write on other cpus: */
		if(my_rank) return;
	}

	/*First write enum: */
	int length=(strlen(this->result_name)+1)*sizeof(char);
	fwrite(&length,sizeof(int),1,fid);
	fwrite(this->result_name,length,1,fid);

	/*Now write time and step: */
	passiveDouble=reCast<IssmPDouble>(time);
	fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
	fwrite(&step,sizeof(int),1,fid);

	/*writing an int array, type is 4 (see parseresultsfromdisk.m):*/
	type=4;
	fwrite(&type,sizeof(int),1,fid);
	rows=this->M;
	fwrite(&rows,sizeof(int),1,fid);
	cols=this->N;
	fwrite(&cols,sizeof(int),1,fid);
	fwrite(value,cols*rows*sizeof(int),1,fid);

}
/*}}}*/
template <> inline int GenericExternalResult<int*>::ObjectEnum(void){ /*{{{*/
	return IntMatExternalResultEnum;
} /*}}}*/
template <> inline void GenericExternalResult<int*>::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int object_enum = this->ObjectEnum();
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->result_name);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->value,M*N);
	marshallhandle->call(this->step);
	marshallhandle->call(this->time);

}  /*}}}*/

/*Specific instantiations for IssmPDouble*: */
template <> inline GenericExternalResult<IssmPDouble*>::GenericExternalResult(int in_id, int in_enum_type,IssmPDouble* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/

	id = in_id;
	M  = in_M;
	N  = in_N;

	EnumToStringx(&this->result_name,in_enum_type);

	step = in_step;
	time = in_time;

	/*Copy result in values*/
	if(M*N){
		value=xNew<IssmPDouble>(M*N);
		xMemCpy<IssmPDouble>(value,in_values,M*N);
	}
	else value=NULL;
}
/*}}}*/
template <> inline GenericExternalResult<IssmPDouble*>::GenericExternalResult(int in_id,const char* in_result_name,IssmPDouble* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/

	id = in_id;
	M  = in_M;
	N  = in_N;

	/*Copy name*/
	this->result_name = xNew<char>(strlen(in_result_name)+1);
	xMemCpy<char>(this->result_name,in_result_name,strlen(in_result_name)+1);

	step = in_step;
	time = in_time;

	/*Copy result in values*/
	if(M*N){
		value=xNew<IssmPDouble>(M*N);
		xMemCpy<IssmPDouble>(value,in_values,M*N);
	}
	else value=NULL;
}
/*}}}*/
template <> inline GenericExternalResult<IssmPDouble*>::GenericExternalResult(int in_id, int in_enum_type,IssmPDouble* in_value,int in_step, IssmDouble in_time){ /*{{{*/
	_error_("you cannot initialize a GenericExternalResult<IssmPDouble*> without providing the dimensions of the matrix! Please use a more appropriate constructor!");
} /*}}}*/
template <> inline GenericExternalResult<IssmPDouble*>::~GenericExternalResult(){ /*{{{*/
	xDelete<char>(result_name);
	xDelete<IssmPDouble>(value);
} /*}}}*/
template <> inline void GenericExternalResult<IssmPDouble*>::Echo(void){ /*{{{*/

	_printf_("GenericExternalResult<IssmPDouble*>:\n");
	this->GenericEcho();
	_printf_("   matrix size: " << this->M << "-" << this->N << "\n");

} /*}}}*/
template <> inline void GenericExternalResult<IssmPDouble*>::DeepEcho(void){ /*{{{*/

	int i,j;

	_printf_("GenericExternalResult<IssmPDouble*>:\n");
	this->GenericEcho();

	_printf_("   matrix size: " << this->M << "-" << this->N << "\n");
	for (i=0;i<this->M;i++){  
		_printf_("   [ ");
		for (j=0;j<this->N;j++){
			_printf_( " " << setw(11) << setprecision (5) << this->value[i*this->N+j]);
		}  
		_printf_(" ]\n");
	}  

} /*}}}*/
template <> inline Object* GenericExternalResult<IssmPDouble*>::copy(void){ /*{{{*/
	return new GenericExternalResult<IssmPDouble*>(this->id,StringToEnumx(this->result_name),this->value,this->M,this->N,this->step,this->time);
} /*}}}*/
template <> inline void GenericExternalResult<IssmPDouble*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

	int     my_rank;
	int     type;
	int     rows,cols;
	char   *name    = NULL;
	IssmPDouble passiveDouble;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if(io_gather){
		/*we are gathering the data on cpu 0, don't write on other cpus: */
		if(my_rank) return;
	}

	/*First write enum: */
	int length=(strlen(this->result_name)+1)*sizeof(char);
	fwrite(&length,sizeof(int),1,fid);
	fwrite(this->result_name,length,1,fid);

	/*Now write time and step: */
	passiveDouble=reCast<IssmPDouble>(time);
	fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
	fwrite(&step,sizeof(int),1,fid);

	/*writing a IssmDouble array, type is 3:*/
	type=3;
	fwrite(&type,sizeof(int),1,fid);
	rows=this->M;
	fwrite(&rows,sizeof(int),1,fid);
	cols=this->N;
	fwrite(&cols,sizeof(int),1,fid);
	fwrite(value,cols*rows*sizeof(IssmPDouble),1,fid);

}
/*}}}*/
template <> inline int GenericExternalResult<IssmPDouble*>::ObjectEnum(void){ /*{{{*/
	return DoubleMatExternalResultEnum;
} /*}}}*/
template <> inline double* GenericExternalResult<IssmPDouble*>::GetValues(void){ /*{{{*/
	return value;
} /*}}}*/
template <> inline void GenericExternalResult<IssmPDouble*>::Marshall(MarshallHandle* marshallhandle){/*{{{*/

	int object_enum = this->ObjectEnum();
	marshallhandle->call(object_enum);

	marshallhandle->call(this->id);
	marshallhandle->call(this->result_name);
	marshallhandle->call(this->M);
	marshallhandle->call(this->N);
	marshallhandle->call(this->value,M*N);
	marshallhandle->call(this->step);
	marshallhandle->call(this->time);

}  /*}}}*/
template <> inline void GenericExternalResult<IssmPDouble*>::Transpose(void){/*{{{*/

	/*Perform transpose only if we have a matrix*/
	if(M>1 && N>1){
		IssmPDouble* temp=xNew<IssmPDouble>(M*N);
		for(int i=0;i<M;i++){
			for(int j=0;j<N;j++){
				temp[j*M+i] = value[i*N+j];
			}
		}
		xDelete<IssmPDouble>(this->value);
		this->value = temp;
	}

	/*Switch dimensions*/
	int temp2 = this->N;
	this->N = this->M;
	this->M = temp2;

} /*}}}*/

	/*Specific instantiations for IssmDouble*: */
#if defined(_HAVE_AD_) && !defined(_WRAPPERS_)  //We hook off this specific specialization when not running ADOLC, otherwise we get a redeclaration with the next specialization. 
template <> inline GenericExternalResult<IssmDouble*>::~GenericExternalResult(){ /*{{{*/
	xDelete<char>(result_name);
	xDelete<IssmDouble>(value);
} /*}}}*/
	template <> inline void GenericExternalResult<IssmDouble*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

		int     i;
		int     my_rank;
		int     type;
		int     rows,cols;
		char   *name    = NULL;
		IssmPDouble passiveDouble;
		IssmPDouble* passiveDoubles;

		/*recover my_rank:*/
		my_rank=IssmComm::GetRank();

		if(io_gather){
			/*we are gathering the data on cpu 0, don't write on other cpus: */
			if(my_rank) return;
		}

		/*First write enum: */
		int length=(strlen(this->result_name)+1)*sizeof(char);
		fwrite(&length,sizeof(int),1,fid);
		fwrite(this->result_name,length,1,fid);

		/*Now write time and step: */
		passiveDouble=reCast<IssmPDouble>(time);
		fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
		fwrite(&step,sizeof(int),1,fid);

		/*writing a IssmDouble array, type is 3:*/
		type=3;
		fwrite(&type,sizeof(int),1,fid);
		rows=this->M;
		fwrite(&rows,sizeof(int),1,fid);
		cols=this->N;
		fwrite(&cols,sizeof(int),1,fid);

		passiveDoubles=xNew<IssmPDouble>(this->M*this->N);
		for (i=0;i<this->M*this->N;i++)passiveDoubles[i]=reCast<IssmPDouble>(value[i]);
		fwrite(passiveDoubles,cols*rows*sizeof(IssmPDouble),1,fid);
		xDelete<IssmPDouble>(passiveDoubles);

	}
	/*}}}*/
#endif

	/*Specific instantiations for IssmComplex*: */
template <> inline GenericExternalResult<IssmComplex*>::GenericExternalResult(int in_id, int in_enum_type,IssmComplex* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/

	id = in_id;
	M  = in_M;
	N  = in_N;

	EnumToStringx(&this->result_name,in_enum_type);

	step = in_step;
	time = in_time;

	/*Copy result in values*/
	if(M*N){
		value=xNew<IssmComplex>(M*N);
		xMemCpy<IssmComplex>(value,in_values,M*N);
	}
	else value=NULL;
}
/*}}}*/
template <> inline void GenericExternalResult<IssmComplex*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

	int     my_rank;
	int     type;
	int     rows,cols;
	char   *name    = NULL;
	IssmPDouble passiveDouble;
	IssmDouble* reals=NULL;
	IssmDouble* imags=NULL;

	/*recover my_rank:*/
	my_rank=IssmComm::GetRank();

	if(io_gather){
		/*we are gathering the data on cpu 0, don't write on other cpus: */
		if(my_rank) return;
	}

	/*First write enum: */
	int length=(strlen(this->result_name)+1)*sizeof(char);
	fwrite(&length,sizeof(int),1,fid);
	fwrite(this->result_name,length,1,fid);

	/*Now write time and step: */
	passiveDouble=reCast<IssmPDouble>(time);
	fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
	fwrite(&step,sizeof(int),1,fid);

	/*writing a IssmComplex array, type is 3:*/
	type=5;
	fwrite(&type,sizeof(int),1,fid);
	rows=this->M;
	fwrite(&rows,sizeof(int),1,fid);
	cols=this->N;
	fwrite(&cols,sizeof(int),1,fid);

	/*write complex array into two real arrays:*/
	reals=xNew<IssmDouble>(cols*rows);
	imags=xNew<IssmDouble>(cols*rows);
	for(int i=0;i<cols*rows;i++){
		reals[i]=value[i].real();
		imags[i]=value[i].imag();
	}
	fwrite(reals,cols*rows*sizeof(IssmComplex),1,fid);
	fwrite(imags,cols*rows*sizeof(IssmComplex),1,fid);

}
/*}}}*/

	/*Specifics instantiations for Vector*/
	template <> inline GenericExternalResult<Vector<IssmPDouble>*>::GenericExternalResult(int in_id, int in_enum_type,Vector<IssmPDouble>* in_values, int in_M,int in_N,int in_step,IssmDouble in_time){/*{{{*/
		_error_("instanciation not correct");
	}
	/*}}}*/
	template <> inline GenericExternalResult<Vector<IssmPDouble>*>::GenericExternalResult(int in_id, int in_enum_type,Vector<IssmPDouble>* in_value,int in_step, IssmDouble in_time){ /*{{{*/
		id = in_id;
		M  = 0;
		N  = 0;

		/*Convert enum to name*/
		EnumToStringx(&this->result_name,in_enum_type);

		step = in_step;
		time = in_time;

		value = in_value;
	} /*}}}*/
	template <> inline GenericExternalResult<Vector<IssmPDouble>*>::~GenericExternalResult(){ /*{{{*/
		xDelete<char>(this->result_name);
		delete value;
	} /*}}}*/
	template <> inline void GenericExternalResult<Vector<IssmPDouble>*>::Echo(void){ /*{{{*/

		_printf_("GenericExternalResult<Vector<IssmPDouble>*>:\n");
		this->GenericEcho();
		this->value->Echo();

	} /*}}}*/
	template <> inline void GenericExternalResult<Vector<IssmPDouble>*>::DeepEcho(void){ /*{{{*/

		this->Echo();

	} /*}}}*/
	template <> inline Object* GenericExternalResult<Vector<IssmPDouble>*>::copy(void){ /*{{{*/
		return new GenericExternalResult<Vector<IssmPDouble>*>(this->id,StringToEnumx(this->result_name),this->value,this->step,this->time);
	} /*}}}*/
#if defined(_HAVE_AD_) && !defined(_WRAPPERS_)  //We hook off this specific specialization when not running ADOLC, otherwise we get a redeclaration with the next specialization.
	template <> inline void GenericExternalResult<Vector<IssmPDouble>*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

		char *name   = NULL;
		int   length,rows,cols=1;

		if(!io_gather){
			_error_("not supported yet");
		}

		/*Serialize vector on cpu0*/
		IssmPDouble* serialvalues = this->value->ToMPISerial0();

		if(IssmComm::GetRank()==0){
			this->value->GetSize(&rows);

			/*First write name: */
			length=(strlen(this->result_name)+1)*sizeof(char);
			fwrite(&length,sizeof(int),1,fid);
			fwrite(this->result_name,length,1,fid);

			/*Now write time and step: */
			IssmPDouble passiveDouble=reCast<IssmPDouble>(time);
			fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
			fwrite(&step,sizeof(int),1,fid);

			/*writing a IssmDouble array, type is 3:*/
			int type=3;
			fwrite(&type,sizeof(int),1,fid);
			fwrite(&rows,sizeof(int),1,fid);
			fwrite(&cols,sizeof(int),1,fid);
			fwrite(serialvalues,cols*rows*sizeof(IssmPDouble),1,fid);
		}

		/*Clean up*/
		xDelete<IssmPDouble>(serialvalues);

	}
	/*}}}*/
template <> inline GenericExternalResult<Vector<IssmDouble>*>::~GenericExternalResult(){ /*{{{*/
	xDelete<char>(this->result_name);
	delete value;
} /*}}}*/
#endif
	template <> inline int GenericExternalResult<Vector<IssmPDouble>*>::ObjectEnum(void){ /*{{{*/
		return NoneEnum;
		/*???? FIXME*/
	} /*}}}*/

	/*Specifics instantiations for Vector<IssmDouble>*/
	template <> inline void GenericExternalResult<Vector<IssmDouble>*>::WriteData(FILE* fid,bool io_gather){ /*{{{*/

		int i;
		char *name   = NULL;
		int   length,rows,cols=1;
		IssmDouble*  serialvalues = NULL;
		IssmPDouble* pserialvalues = NULL;

		if(!io_gather){
			_error_("not supported yet");
		}

		/*Serialize vector only on cpu0*/
		serialvalues = this->value->ToMPISerial0();

		if(IssmComm::GetRank()==0){

			/*Make it passive*/
			this->value->GetSize(&rows);
			pserialvalues=xNew<IssmPDouble>(rows);
			for(i=0;i<rows;i++)pserialvalues[i]=reCast<IssmPDouble>(serialvalues[i]);

			/*First write name: */
			length=(strlen(this->result_name)+1)*sizeof(char);
			fwrite(&length,sizeof(int),1,fid);
			fwrite(this->result_name,length,1,fid);

			/*Now write time and step: */
			IssmPDouble passiveDouble=reCast<IssmPDouble>(time);
			fwrite(&passiveDouble,sizeof(IssmPDouble),1,fid);
			fwrite(&step,sizeof(int),1,fid);

			/*writing a IssmDouble array, type is 3:*/
			int type=3;
			fwrite(&type,sizeof(int),1,fid);
			fwrite(&rows,sizeof(int),1,fid);
			fwrite(&cols,sizeof(int),1,fid);
			fwrite(pserialvalues,cols*rows*sizeof(IssmPDouble),1,fid);

			/*Clean up*/
			xDelete<IssmPDouble>(pserialvalues);
		}

		/*Clean up*/
		xDelete<IssmDouble>(serialvalues);
	}
	/*}}}*/
	template <> inline void GenericExternalResult<Vector<IssmDouble>*>::Marshall(MarshallHandle* marshallhandle){/*{{{*/
		_error_("GenericExternalResult instantiated for type Vector<IssmDouble>* called " << result_name << " not implemented yet");
	}  /*}}}*/

#endif  /* _EXTERNAL_RESULTOBJECT_H */
