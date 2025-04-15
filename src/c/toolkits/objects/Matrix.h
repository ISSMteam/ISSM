/*!\file:  Matrix.h
 * \brief wrapper to matrix objects. The goal is to control which API (PETSc,Scalpack, Plapack?) 
 * implements our underlying matrix format.
 */ 

#ifndef _MATRIX_H_
#define _MATRIX_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif
#include <cstring>
#include "../../shared/Enum/Enum.h"
#include "../petsc/petscincludes.h"
#include "../issm/issmtoolkit.h"
/*}}}*/

enum matrixtype { PetscMatType, IssmMatType };

template <class doubletype> class Vector;

template <class doubletype> 
class Matrix{

	public:

		int       type;
		#ifdef _HAVE_PETSC_
		PetscMat<doubletype>  *pmatrix;
		#endif
		IssmMat<doubletype>   *imatrix;

		/*Matrix constructors, destructors*/
		Matrix(){/*{{{*/
			InitCheckAndSetType();
		}
		/*}}}*/
		Matrix(int M,int N){/*{{{*/

			InitCheckAndSetType();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix=new PetscMat<doubletype>(M,N);
				#endif
			}
			else{
				this->imatrix=new IssmMat<doubletype>(M,N);
			}

		}
		/*}}}*/
		Matrix(int m,int n,int M,int N,int* d_nnz,int* o_nnz){/*{{{*/

			InitCheckAndSetType();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix=new PetscMat<doubletype>(m,n,M,N,d_nnz,o_nnz);
				#endif
			}
			else{
				this->imatrix=new IssmMat<doubletype>(m,n,M,N,d_nnz,o_nnz);
			}

		}
		/*}}}*/
		Matrix(int M,int N,double sparsity){/*{{{*/

			InitCheckAndSetType();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix=new PetscMat<doubletype>(M,N,sparsity);
				#endif
			}
			else{
				this->imatrix=new IssmMat<doubletype>(M,N,sparsity);
			}
		}
		/*}}}*/
		Matrix(IssmPDouble* serial_mat,int M,int N,IssmPDouble sparsity){/*{{{*/

			InitCheckAndSetType();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix=new PetscMat<doubletype>(serial_mat,M,N,sparsity);
				#endif
			}
			else{
				this->imatrix=new IssmMat<doubletype>(serial_mat,M,N,sparsity);
			}

		}
		/*}}}*/
		Matrix(int M,int N,int connectivity,int numberofdofspernode){/*{{{*/

			InitCheckAndSetType();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix=new PetscMat<doubletype>(M,N,connectivity,numberofdofspernode);
				#endif
			}
			else{
				this->imatrix=new IssmMat<doubletype>(M,N,connectivity,numberofdofspernode);
			}

		}
		/*}}}*/
		~Matrix(){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				delete this->pmatrix;
				#endif
			}
			else delete this->imatrix;

		}
		/*}}}*/
		void InitCheckAndSetType(void){/*{{{*/

			#ifdef _HAVE_PETSC_
			pmatrix=NULL;
			#endif
			imatrix=NULL;

			/*retrieve toolkittype: */
			char* toolkittype=ToolkitOptions::GetToolkitType();

			/*set matrix type: */
			if (strcmp(toolkittype,"petsc")==0){
				#ifdef _HAVE_PETSC_
				type=PetscMatType; 
				#else
				_error_("cannot create petsc matrix without PETSC compiled!");
				#endif
			}
			else if(strcmp(toolkittype,"issm")==0){
				/*let this choice stand:*/
				type=IssmMatType;
			}
			else{
				_error_("unknow toolkit type ");
			}

			/*Free resources: */
			xDelete<char>(toolkittype);
		}
		/*}}}*/

		/*Matrix specific routines:*/
		void Echo(void){/*{{{*/
			_assert_(this);

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->Echo();
				#endif
			}
			else{
				this->imatrix->Echo();
			}

		}
		/*}}}*/
		void EchoDebug(std::string message){_assert_(this);/*{{{*/

			if(type==PetscMatType){
#ifdef _HAVE_PETSC_
				this->pmatrix->EchoDebug(message);
#endif
			}
			else this->imatrix->EchoDebug(message);
		}
		void AllocationInfo(void){/*{{{*/
			_assert_(this);
			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->AllocationInfo();
				#endif
			}
			else{
				//this->imatrix->AllocationInfo();
				_error_("not supported yet");
			}
		}/*}}}*/
		void Assemble(void){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->Assemble();
				#endif
			}
			else{
				this->imatrix->Assemble();
			}
		}
		/*}}}*/
		IssmDouble Norm(NormMode norm_type){/*{{{*/

			IssmDouble norm=0;

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				norm=this->pmatrix->Norm(norm_type);
				#endif
			}
			else{
				norm=this->imatrix->Norm(norm_type);
			}

			return norm;
		}
		/*}}}*/
		void GetSize(int* pM,int* pN){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->GetSize(pM,pN);
				#endif
			}
			else{
				this->imatrix->GetSize(pM,pN);
			}

		}
		/*}}}*/
		void GetLocalSize(int* pM,int* pN){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->GetLocalSize(pM,pN);
				#endif
			}
			else{
				this->imatrix->GetLocalSize(pM,pN);
			}

		}
		/*}}}*/
		void MatMult(Vector<doubletype>* X,Vector<doubletype>* AX){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->MatMult(X->pvector,AX->pvector);
				#endif
			}
			else{
				this->imatrix->MatMult(X->ivector,AX->ivector);
			}

		}
		/*}}}*/
		Matrix<doubletype>* Duplicate(void){/*{{{*/

			Matrix<doubletype>* output=new Matrix<doubletype>();

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				output->pmatrix=this->pmatrix->Duplicate();
				#endif
			}
			else{
				output->imatrix=this->imatrix->Duplicate();
			}

			return output;
		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/

			doubletype* output=NULL;

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				output=this->pmatrix->ToMPISerial0();
				#endif
			}
			else{
				output=this->imatrix->ToMPISerial0();
			}

			return output;
		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/

			doubletype* output=NULL;

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				output=this->pmatrix->ToMPISerial();
				#endif
			}
			else{
				_error_("not implemented yet!");
			}

			return output;
		}
		/*}}}*/
		void SetValues(int m,int* idxm,int n,int* idxn,IssmDouble* values,InsMode mode){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->SetValues(m,idxm,n,idxn,values,mode);
				#endif
			}
			else{
				this->imatrix->SetValues(m,idxm,n,idxn,values,mode);
			}
		}
		/*}}}*/
		void Convert(MatrixType newtype){/*{{{*/

			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->Convert(newtype);
				#endif
			}
			else{
				this->imatrix->Convert(newtype);
			}

		}
		/*}}}*/
		void SetZero(void) {/*{{{*/
			// sets all values to 0 but keeps the structure of a sparse matrix
			if(type==PetscMatType){
				#ifdef _HAVE_PETSC_
				this->pmatrix->SetZero();
				#endif
			}
			else{
				this->imatrix->SetZero();
			}
		}
		/*}}}*/
};

#endif //#ifndef _MATRIX_H_
