/*!\file DoubleMatArrayParam.c
 * \brief: implementation of the DoubleMatArrayParam object
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

/*DoubleMatArrayParam constructors and destructor*/
DoubleMatArrayParam::DoubleMatArrayParam(){/*{{{*/
	return;
}
/*}}}*/
DoubleMatArrayParam::DoubleMatArrayParam(int in_enum_type,IssmDouble** in_array, int in_M, int* in_mdim_array, int* in_ndim_array){/*{{{*/

	int i;
	IssmDouble* matrix=NULL;
	int     m,n;

	enum_type=in_enum_type;
	M=in_M;
	if(M){
		array=xNew<IssmDouble*>(M);
		mdim_array=xNew<int>(M);
		ndim_array=xNew<int>(M);

		for(i=0;i<M;i++){
			m=in_mdim_array[i]; 
			n=in_ndim_array[i];

			mdim_array[i]=m;
			ndim_array[i]=n;

			if(m*n){
				matrix=xNew<IssmDouble>(m*n);
				xMemCpy<IssmDouble>(matrix,in_array[i],m*n);
			}
			else{
				matrix=NULL;
			}
			array[i]=matrix;
		}
	}
	else{
		array=NULL;
		mdim_array=NULL;
		ndim_array=NULL;
	}
}
/*}}}*/
DoubleMatArrayParam::~DoubleMatArrayParam(){/*{{{*/

	int i;
	IssmDouble* matrix=NULL;

	xDelete<int>(mdim_array);
	xDelete<int>(ndim_array);

	for(i=0;i<M;i++){
		matrix=array[i];
		xDelete<IssmDouble>(matrix);
	}

	xDelete<IssmDouble*>(array);
	return;
}
/*}}}*/

/*Object virtual functions definitions:*/
Param* DoubleMatArrayParam::copy() {/*{{{*/

	return new DoubleMatArrayParam(this->enum_type,this->array, this->M, this->mdim_array,this->ndim_array);

}
/*}}}*/
void DoubleMatArrayParam::DeepEcho(void){/*{{{*/

	int i,j,k;
	int m,n;
	IssmDouble* matrix=NULL;

	_printf_("DoubleMatArrayParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   array size: " << this->M << "\n");
	for(i=0;i<M;i++){
		_printf_("   array " << i << " (" << mdim_array[i] << "x" << ndim_array[i] << "):\n");
		matrix=array[i];
		m=mdim_array[i];
		n=ndim_array[i];

		for(j=0;j<m;j++){
			_printf_("   ");
			for(k=0;k<n;k++)_printf_(*(matrix+n*j+k) << " ");
			_printf_("\n");
		}
	}
}
/*}}}*/
void DoubleMatArrayParam::Echo(void){/*{{{*/

	_printf_("DoubleMatArrayParam:\n");
	_printf_("   enum: " << this->enum_type << " (" << EnumToStringx(this->enum_type) << ")\n");
	_printf_("   array size: " << this->M << "\n");
	_printf_("   array pointer: " << this->array << "\n");

}
/*}}}*/
int    DoubleMatArrayParam::Id(void){ return -1; }/*{{{*/
/*}}}*/
void DoubleMatArrayParam::Marshall(MarshallHandle* marshallhandle){ /*{{{*/

	int object_enum = DoubleMatArrayParamEnum;
	marshallhandle->call(object_enum);

	marshallhandle->call(this->enum_type);
	marshallhandle->call(this->M);
	if(M){
		marshallhandle->call(this->mdim_array,M);
		marshallhandle->call(this->ndim_array,M);
		if(marshallhandle->OperationNumber()==MARSHALLING_LOAD && M) array=xNew<IssmDouble*>(M);
		for(int i=0;i<M;i++){
			marshallhandle->call(this->array[i],mdim_array[i]*ndim_array[i]);
		}
	}
	else{
		array=NULL;
		mdim_array=NULL;
		ndim_array=NULL;
	}
}
/*}}}*/
int DoubleMatArrayParam::ObjectEnum(void){/*{{{*/

	return DoubleMatArrayParamEnum;

}
/*}}}*/

/*DoubleMatArrayParam virtual functions definitions: */
void  DoubleMatArrayParam::GetParameterValue(IssmDouble*** pout_array, int* pout_M,int** pout_mdim_array, int** pout_ndim_array){/*{{{*/

	int i,m,n;
	IssmDouble* matrix=NULL;
	IssmDouble* out_matrix=NULL;

	/*output: */
	IssmDouble** out_array=NULL;
	int      out_M;
	int*     out_mdim_array=NULL;
	int*     out_ndim_array=NULL;

	out_M=this->M;
	if(out_M){
		out_array=xNew<IssmDouble*>(M);
		out_mdim_array=xNew<int>(M);
		out_ndim_array=xNew<int>(M);

		xMemCpy<int>(out_mdim_array,this->mdim_array,M);
		xMemCpy<int>(out_ndim_array,this->ndim_array,M);

		for(i=0;i<this->M;i++){
			matrix=this->array[i];
			m=this->mdim_array[i];
			n=this->ndim_array[i];

			if(m*n){
				out_matrix=xNew<IssmDouble>(m*n);
				xMemCpy<IssmDouble>(out_matrix,matrix,m*n);
			}
			else{
				out_matrix=NULL;
			}
			out_array[i]=out_matrix;
		}
	}
	else{
		out_array=NULL;
		out_matrix=NULL;
		out_ndim_array=NULL;
	}

	/*Assign output pointers:*/
	if(pout_M) *pout_M=out_M;
	if(pout_mdim_array) *pout_mdim_array=out_mdim_array;
	if(pout_ndim_array) *pout_ndim_array=out_ndim_array;
	*pout_array=out_array;

}
/*}}}*/
void  DoubleMatArrayParam::SetValue(IssmDouble** in_array, int in_M, int* in_mdim_array, int* in_ndim_array){/*{{{*/

	int i,m,n;
	IssmDouble* in_matrix=NULL;
	IssmDouble* matrix=NULL;

	/*avoid leak: */
	xDelete<int>(mdim_array);
	xDelete<int>(ndim_array);
	for(i=0;i<M;i++){
		matrix=array[i];
		xDelete<IssmDouble>(matrix);
	}
	xDelete<IssmDouble*>(array);

	/*copy data: */
	this->M=in_M;
	this->array=xNew<IssmDouble*>(M);
	this->mdim_array=xNew<int>(M);
	this->ndim_array=xNew<int>(M);

	xMemCpy<int>(this->mdim_array,in_mdim_array,M);
	xMemCpy<int>(this->ndim_array,in_ndim_array,M);

	for(i=0;i<M;i++){
		in_matrix=in_array[i];
		m=in_mdim_array[i];
		n=in_ndim_array[i];

		matrix=xNew<IssmDouble>(m*n);
		xMemCpy<IssmDouble>(matrix,in_matrix,m*n);

		this->array[i]=matrix;
	}

}
/*}}}*/
