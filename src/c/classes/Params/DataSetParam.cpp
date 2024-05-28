/*!\file DataSetParam.c
 * \brief: implementation of the DataSetParam object
 */

/*header files: */
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../classes.h"
#include "shared/shared.h"
/*}}}*/

/*DataSetParam constructors and destructor*/
DataSetParam::DataSetParam(){/*{{{*/
	value=NULL;
	return;
}
/*}}}*/
DataSetParam::DataSetParam(int in_enum_type,DataSet* in_value){/*{{{*/

	enum_type=in_enum_type;
	value=in_value->Copy();
}
/*}}}*/
DataSetParam::~DataSetParam(){/*{{{*/
	delete value;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* DataSetParam::copy() {/*{{{*/

	return new DataSetParam(this->enum_type,this->value);

}
/*}}}*/
void DataSetParam::DeepEcho(void){/*{{{*/

	_printf_(setw(22)<<"   DataSetParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" ----- begin\n");
	this->value->Echo();
	_printf_(setw(22)<<"   DataSetParam "<<setw(35)<<left<<EnumToStringx(this->enum_type)<<" ----- end\n");
}
/*}}}*/
void DataSetParam::Echo(void){/*{{{*/
	this->DeepEcho();
}
/*}}}*/
int    DataSetParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DataSetParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	if(marshallhandle->OperationNumber()==MARSHALLING_LOAD)value=new DataSet();

	int object_enum=DataSetParamEnum;
	marshallhandle->call(object_enum);
	marshallhandle->call(this->enum_type);
	value->Marshall(marshallhandle);

}
/*}}}*/
int DataSetParam::ObjectEnum(void){/*{{{*/

	return DataSetParamEnum;

}
/*}}}*/

/*DataSetParam virtual functions definitions: */
void DataSetParam::GetParameterValue(DataSet** pdataset){/*{{{*/
	*pdataset=value->Copy();
}
/*}}}*/
void DataSetParam::SetValue(DataSet* dataset){/*{{{*/
	this->value=dataset;
}
/*}}}*/
