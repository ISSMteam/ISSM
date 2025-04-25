/*!\file PetscVec.cpp
 * \brief: implementation of the PetscVec object
 */

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <stdio.h>
#include <string.h>
#include "../petscincludes.h"
#include "../../../shared/shared.h"

#ifdef _HAVE_CODIPACK_
#include "../../codipack/CoDiPackDebug.h"
#endif

#include "PetscVec.h"
/*}}}*/

/*PetscVec constructors and destructor*/
template<typename doubletype>
PetscVec<doubletype>::PetscVec(){/*{{{*/
	this->vector=NULL;
}
/*}}}*/
template<typename doubletype>
PetscVec<doubletype>::PetscVec(int M,bool fromlocalsize){/*{{{*/

	this->vector=NewVec<PVec>(M,IssmComm::GetComm(),fromlocalsize);
}
/*}}}*/
template<typename doubletype>
PetscVec<doubletype>::PetscVec(int m,int M){/*{{{*/

	VecCreate(IssmComm::GetComm(),&this->vector);
	VecSetSizes(this->vector,m,M);
	VecSetFromOptions(this->vector);
	VecSetOption(this->vector,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
	VecSet(this->vector,0.0);
}
/*}}}*/
template<typename doubletype>
PetscVec<doubletype>::PetscVec(PVec petsc_vec){/*{{{*/

	if(petsc_vec==NULL){
		this->vector=NewVec<PVec>(0,IssmComm::GetComm());
	}
	else{
		/*copy vector*/
		VecDuplicate(petsc_vec,&this->vector);
		VecCopy(petsc_vec,this->vector);
	}

}
/*}}}*/
template<typename doubletype>
PetscVec<doubletype>::PetscVec(doubletype* serial_vec,int M){/*{{{*/

	int* idxm=NULL;
	if(M)idxm=xNew<int>(M);
	for(int i=0;i<M;i++) idxm[i]=i;

	this->vector=NewVec<PVec>(M,IssmComm::GetComm());
	VecSetValues(this->vector,M,idxm,serial_vec,INSERT_VALUES);
	VecAssemblyBegin(this->vector);
	VecAssemblyEnd(this->vector);

	xDelete<int>(idxm);
}
/*}}}*/
template<typename doubletype>
PetscVec<doubletype>::~PetscVec(){/*{{{*/
    VecFree(&this->vector);
}
/*}}}*/

/*PetscVec specific routines: */
template<typename doubletype>
void PetscVec<doubletype>::Echo(void){/*{{{*/

	_assert_(this->vector);
	VecView(this->vector,PETSC_VIEWER_STDOUT_WORLD);
}
/*}}}*/

template<typename doubletype>
void PetscVec<doubletype>::EchoDebug(std::string message){/*{{{*/
#if defined(_HAVE_CODIPACK_) & defined(_HAVE_ADJOINTPETSC_)
	if (std::is_same<doubletype, IssmDouble>::value && CoDiIsDebugOutput()) {
		adjoint_petsc::ADVecDebugOutput(this->vector, message, CoDiGetUniqueID());
	}
#endif
}
/*}}}*/

template<typename doubletype>
void PetscVec<doubletype>::Assemble(void){/*{{{*/

	_assert_(this->vector);
	VecAssemblyBegin(this->vector);
	VecAssemblyEnd(this->vector);
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::SetValues(int ssize, int* list, doubletype* values, InsMode mode){/*{{{*/

	_assert_(this->vector);
	VecSetValues(this->vector,ssize,list,values,ISSMToPetscInsertMode(mode));

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::SetValue(int dof, doubletype value, InsMode mode){/*{{{*/

	SetValues(1, &dof, &value, mode);
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::GetValue(doubletype* pvalue,int dof){/*{{{*/

	_assert_(this->vector);
	VecGetValues(this->vector,1,&dof,pvalue);
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::GetSize(int* pM){/*{{{*/

	_assert_(this->vector);
	VecGetSize(this->vector,pM);
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::GetLocalSize(int* pm){/*{{{*/

	_assert_(this->vector);
	VecGetLocalSize(this->vector,pm);
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::GetLocalVector(doubletype** pvector,int** pindices){/*{{{*/

	_assert_(this->vector);

	/*First, check that vector size is not 0*/
	int vector_size;
	this->GetSize(&vector_size);
	if(vector_size==0){
		*pvector=NULL;
		*pindices=NULL;
		return;
	}

	/*Get Ownership range*/
	PetscInt lower_row,upper_row;
	VecGetOwnershipRange(this->vector,&lower_row,&upper_row);
	int range=upper_row-lower_row;

	/*return NULL if no range*/
	if(range==0){
		*pvector=NULL;
		*pindices=NULL;
		return;
	}

	/*Build indices*/
	int* indices=xNew<int>(range);
	for(int i=0;i<range;i++) indices[i]=lower_row+i;
	/*Get vector*/
	doubletype* values =xNew<doubletype>(range);
	VecGetValues(this->vector,range,indices,values);

	*pvector  = values;
	*pindices = indices;
} /*}}}*/
template<typename doubletype>
PetscVec<doubletype>* PetscVec<doubletype>::Duplicate(void){/*{{{*/

	_assert_(this->vector);

	/*Instantiate output Vector*/
	PetscVec* output=new PetscVec();

	/*Duplicate using "this" layout*/
	VecDuplicate(this->vector,&output->vector);

	/*Return new vector*/
	return output;
}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Set(doubletype value){/*{{{*/

	_assert_(this->vector);
	VecSet(this->vector,value);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::AXPY(PetscVec* X, doubletype a){/*{{{*/

	_assert_(this->vector);
	VecAXPY(this->vector,a,X->vector);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::AYPX(PetscVec* X, doubletype a){/*{{{*/

	_assert_(this->vector);
	VecAYPX(this->vector,a,X->vector);

}
/*}}}*/
template<typename doubletype>
doubletype* PetscVec<doubletype>::ToMPISerial(void){/*{{{*/

	doubletype* vec_serial=NULL;
	VecToMPISerial(&vec_serial, this->vector,IssmComm::GetComm(),true);
	return vec_serial;

}
/*}}}*/
template<typename doubletype>
doubletype* PetscVec<doubletype>::ToMPISerial0(void){/*{{{*/

	doubletype* vec_serial=NULL;
	VecToMPISerial(&vec_serial, this->vector,IssmComm::GetComm(),false);
	return vec_serial;

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Shift(doubletype shift){/*{{{*/

	if(this->vector) VecShift(this->vector,shift);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Copy(PetscVec* to){/*{{{*/

	if(this->vector) VecCopy(this->vector,to->vector);

}
/*}}}*/
template<typename doubletype>
doubletype PetscVec<doubletype>::Max(void){/*{{{*/

	_assert_(this->vector);

	doubletype max;
	VecMax(this->vector,NULL,&max);
	return max;

}
/*}}}*/
template<typename doubletype>
doubletype PetscVec<doubletype>::Norm(NormMode mode){/*{{{*/

	doubletype norm=0;
	_assert_(this->vector);
	VecNorm(this->vector,ISSMToPetscNormMode(mode),&norm);
	return norm;

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Scale(doubletype scale_factor){/*{{{*/

	_assert_(this->vector);
	VecScale(this->vector,scale_factor);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Pow(doubletype scale_factor){/*{{{*/

	_assert_(this->vector);
	VecPow(this->vector,scale_factor);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::Sum(doubletype* pvalue){/*{{{*/

	_assert_(this->vector);
	VecSum(this->vector,pvalue);

}
/*}}}*/
template<typename doubletype>
doubletype PetscVec<doubletype>::Dot(PetscVec* input){/*{{{*/

	doubletype dot;
	_assert_(this->vector);
	VecDot(this->vector,input->vector,&dot);
	return dot;

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::PointwiseDivide(PetscVec* x,PetscVec* y){/*{{{*/

	_assert_(this->vector);
	VecPointwiseDivide(this->vector,x->vector,y->vector);

}
/*}}}*/
template<typename doubletype>
void PetscVec<doubletype>::PointwiseMult(PetscVec* x,PetscVec* y){/*{{{*/

	_assert_(this->vector);
	VecPointwiseMult(this->vector,x->vector,y->vector);

}
/*}}}*/

// Explicit instantiations.
template class PetscVec<IssmDouble>;
#if _HAVE_CODIPACK_
template class PetscVec<IssmPDouble>;
#endif
