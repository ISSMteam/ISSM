/*!\file:  IssmMat.h
 * \brief Main Matrix class for the Issm toolkit. 
 */ 

#ifndef _ISSM_MAT_H_
#define _ISSM_MAT_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include "../ToolkitOptions.h"
#include "./IssmToolkitUtils.h"
#include "../mumps/mumpsincludes.h"
/*}}}*/

/*We need to template this class, in case we want to create Matrices that hold
  IssmDouble* matrix or IssmPDouble* matrix. 
  Such matrices are useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)
*/
template <class doubletype> class IssmVec;
template <class doubletype> class IssmDenseMat;
template <class doubletype> class IssmMpiDenseMat;
template <class doubletype> class IssmMpiSparseMat;
class Parameters;

template <class doubletype> 
class IssmMat{

	public:

		IssmAbsMat<doubletype>* matrix; /*abstract matrix, which implements object orientation*/

		/*IssmMat constructors, destructors*/
		IssmMat(){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>();
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>();
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>();
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}
		}
		/*}}}*/
		IssmMat(int M,int N){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>(M,N);
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>(M,N);
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>(M,N);
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}
		}
		/*}}}*/
		IssmMat(int M,int N, doubletype sparsity){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>(M,N,sparsity);
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>(M,N,sparsity);
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>(M,N,sparsity);
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}
		}
		/*}}}*/
		IssmMat(int m,int n,int M,int N,int* d_nnz,int* o_nnz){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>(m,n,M,N,d_nnz,o_nnz);
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>(m,n,M,N,d_nnz,o_nnz);
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>(m,n,M,N,d_nnz,o_nnz);
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}
		}
		/*}}}*/
		IssmMat(doubletype* serial_mat,int M,int N,doubletype sparsity){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>(serial_mat,M,N,sparsity);
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>(serial_mat,M,N,sparsity);
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>(serial_mat,M,N,sparsity);
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}

		}
		/*}}}*/
		IssmMat(int M,int N, int connectivity, int numberofdofspernode){ /*{{{*/

			switch(IssmMatTypeFromToolkitOptions()){

				case DenseEnum: 
					this->matrix=new IssmDenseMat<doubletype>(M,N,connectivity,numberofdofspernode);
					break;
				case MpiDenseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiDenseMat<doubletype>(M,N,connectivity,numberofdofspernode);
					#else
					_error_("MpiDense matrix requires compilation of MPI!");
					#endif
					break;
				case MpiSparseEnum:
					#ifdef _HAVE_MPI_
					this->matrix=new IssmMpiSparseMat<doubletype>(M,N,connectivity,numberofdofspernode);
					#else
					_error_("MpiSparse matrix requires compilation of MPI!");
					#endif
					break;
				default:
					_error_("matrix type not supported yet!");
			}
		}
		/*}}}*/
		~IssmMat(){ /*{{{*/
			delete matrix;
		} /*}}}*/

		/*Functionality: */
		void Echo(void){  /*{{{*/
			matrix->Echo();
		} /*}}}*/
		void EchoDebug(std::string message){  /*{{{*/
			matrix->EchoDebug(message);
		} /*}}}*/
		void Assemble(void){  /*{{{*/
			matrix->Assemble();
		} /*}}}*/
		doubletype Norm(NormMode mode){ /*{{{*/
			return matrix->Norm(mode);
		}
		/*}}}*/
		void GetSize(int* pM,int* pN){ /*{{{*/
			matrix->GetSize(pM,pN);
		} /*}}}*/
		void GetLocalSize(int* pM,int* pN){ /*{{{*/
			matrix->GetLocalSize(pM,pN);
		} /*}}}*/
		void MatMult(IssmVec<doubletype>* X,IssmVec<doubletype>* AX){ /*{{{*/
			matrix->MatMult(X->vector,AX->vector);
		} /*}}}*/
		IssmMat<doubletype>* Duplicate(void){ /*{{{*/

			IssmMat<doubletype>* issmmatrix=NULL;

			issmmatrix=new IssmMat<doubletype>();
			issmmatrix->matrix=this->matrix->Duplicate();

			return issmmatrix;
		} /*}}}*/
		doubletype* ToSerial(void){/*{{{*/
			return matrix->ToSerial();
		}/*}}}*/
		void SetValues(int m,int* idxm,int n,int* idxn,doubletype* values,InsMode mode){/*{{{*/
			matrix->SetValues(m,idxm,n,idxn,values,mode);
		}/*}}}*/
		void Convert(MatrixType type){/*{{{*/
			matrix->convert(type);
		}/*}}}*/
		void SetZero(void){/*{{{*/
			matrix->SetZero();
		}/*}}}*/
		#ifndef _HAVE_WRAPPERS_
		IssmVec<doubletype>* Solve(IssmVec<doubletype>* pf, Parameters* parameters){ /*{{{*/

			IssmVec<doubletype>* outvector=NULL;

			outvector=new IssmVec<doubletype>();

			outvector->vector=this->matrix->Solve(pf->vector,parameters);

			return outvector;

		}/*}}}*/
		#endif
};

#endif //#ifndef _ISSMMAT_H_
