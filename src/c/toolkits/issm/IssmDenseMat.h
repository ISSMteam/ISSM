/*!\file:  IssmDenseMat.h
 * \brief implementation of an ISSM matrix which run serially (1 cpu only), which is made of a fully dense 
 * matrix. Internally, this dense matrix is just a linear buffer of type doubletype. 
 * This object needs to answer the API defined by the virtual functions in IssmAbsMat, 
 * and the contructors required by IssmMat (see IssmMat.h)
 */ 

#ifndef _ISSM_DENSE_MAT_H_
#define _ISSM_DENSE_MAT_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "./IssmSeqVec.h"
#include "./IssmToolkitUtils.h"
#include "../../shared/shared.h"
#include "../gsl/gslincludes.h"

#include <math.h>

/*}}}*/

/*We need to template this class, in case we want to create Matrices that hold
  IssmDouble* matrix or IssmPDouble* matrix. 
  Such matrices would be useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)*/

template <class doubletype> class IssmAbsVec;
template <class doubletype> class IssmAbsMat;
template <class doubletype> class IssmSeqVec;

template <class doubletype> 
class IssmDenseMat: public IssmAbsMat<doubletype>{

	public:

		int M,N; 
		doubletype* matrix;  /*here, doubletype is either IssmDouble or IssmPDouble*/

		/*IssmDenseMat constructors, destructors*/
		/*IssmDenseMat(){{{*/
		IssmDenseMat(){

			this->M=0;
			this->N=0;
			this->matrix=NULL;
		}
		/*}}}*/
		/*IssmDenseMat(int M,int N){{{*/
		IssmDenseMat(int pM,int pN){

			this->M=pM;
			this->N=pN;
			this->matrix=NULL;
			if(M*N) this->matrix=xNewZeroInit<doubletype>(pM*pN);
		}
		/*}}}*/
		/*IssmDenseMat(int M,int N, doubletype sparsity){{{*/
		IssmDenseMat(int pM,int pN, doubletype sparsity){

			this->M=pM;
			this->N=pN;
			this->matrix=NULL;
			if(M*N) this->matrix=xNewZeroInit<doubletype>(pM*pN);
		}
		/*}}}*/
		/*IssmDenseMat(int m,int n,int M,int N,int* d_nnz,int* o_nnz){{{*/
		IssmDenseMat(int m,int n,int pM,int pN,int* d_nnz,int* o_nnz){

			this->M=pM;
			this->N=pN;
			this->matrix=NULL;
			if(pM*pN) this->matrix=xNewZeroInit<doubletype>(pM*pN);
		}
		/*}}}*/
		/*IssmDenseMat(doubletype* serial_mat,int M,int N,doubletype sparsity){{{*/
		IssmDenseMat(doubletype* serial_mat,int pM,int pN,doubletype sparsity){

			this->M=pM;
			this->N=pN;
			this->matrix=NULL;
			if(M*N){
				this->matrix=xNewZeroInit<doubletype>(pM*pN);
				xMemCpy<doubletype>(this->matrix,serial_mat,pM*pN);
			}

		}
		/*}}}*/
		/*IssmDenseMat(int M,int N, int connectivity, int numberofdofspernode){{{*/
		IssmDenseMat(int pM,int pN, int connectivity,int numberofdofspernode){

			this->M=pM;
			this->N=pN;
			this->matrix=NULL;
			if(M*N) this->matrix=xNewZeroInit<doubletype>(pM*pN);
		}
		/*}}}*/
		/*~IssmDenseMat(){{{*/
		~IssmDenseMat(){

			xDelete<doubletype>(this->matrix);
			M=0;
			N=0;
		}
		/*}}}*/

		/*IssmAbsMat virtual functions*/
		void Echo(void){/*{{{*/

			_printf_("IssmDenseMat size " << this->M << "-" << this->N << "\n");
			for(int i=0;i<M;i++){
				for(int j=0;j<N;j++){
					_printf_(this->matrix[N*i+j] << " ");
				}
				_printf_("\n");
			}
		}
		/*}}}*/
		void EchoDebug(std::string message) {
			std::cout << "Error: not implemented." << std::endl;
		}
		void Assemble(void){/*{{{*/

			/*do nothing*/

		}
		/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/

			doubletype norm;
			doubletype absolute;
			int i,j;

			switch(mode){
				case NORM_INF:
					norm=0.;
					for(i=0;i<this->M;i++){
						absolute=0;
						for(j=0;j<this->N;j++){
							absolute+=fabs(this->matrix[N*i+j]);
						}
						norm=max(norm,absolute);
					}
					return norm;
					break; 
				case NORM_FROB:
					norm=0.;
					for(i=0;i<this->M;i++){
						for(j=0;j<this->N;j++){
							norm+=this->matrix[N*i+j]*this->matrix[N*i+j];
						}
					}
					return sqrt(norm);
					break; 

				default:
					_error_("unknown norm !");
					break;
			}

			return 0.;
		}
		/*}}}*/
		void GetSize(int* pM,int* pN){/*{{{*/
			*pM=this->M;
			*pN=this->N;
		}
		/*}}}*/
		void GetLocalSize(int* pM,int* pN){/*{{{*/

			*pM=this->M;
			*pN=this->N;

		}
		/*}}}*/
		void MatMult(IssmAbsVec<doubletype>* Xin,IssmAbsVec<doubletype>* AXin){/*{{{*/

			/*We assume that the vectors coming in are of compatible type: */
			int i,j;
			int XM,AXM;
			doubletype dummy;
			IssmSeqVec<doubletype>* X=NULL;
			IssmSeqVec<doubletype>* AX=NULL;

			/*downcast X and AX: */
			X=(IssmSeqVec<doubletype>*)Xin;
			AX=(IssmSeqVec<doubletype>*)AXin;

			/*Some checks first: */
			X->GetSize(&XM);
			AX->GetSize(&AXM);

			if(M!=AXM)_error_("A and AX should have the same number of rows!");
			if(N!=XM)_error_("A and X should have the same number of columns!");

			for(i=0;i<M;i++){
				dummy=0;
				for(j=0;j<N;j++){
					dummy+= this->matrix[N*i+j]*X->vector[j];
				}
				AX->vector[i]=dummy;
			}

		}
		/*}}}*/
		IssmDenseMat<doubletype>* Duplicate(void){/*{{{*/

			doubletype dummy=0;

			return new IssmDenseMat<doubletype>(this->matrix,this->M,this->N,dummy);

		}
		/*}}}*/
		doubletype* ToSerial(void){/*{{{*/

			doubletype* buffer=NULL;

			if(this->M*this->N){
				buffer=xNew<doubletype>(this->M*this->N);
				xMemCpy<doubletype>(buffer,this->matrix,this->M*this->N);
			}
			return buffer;

		}
		/*}}}*/
		void SetValues(int m,int* idxm,int n,int* idxn,doubletype* values,InsMode mode){/*{{{*/

			int i,j;
			switch(mode){
				case ADD_VAL:
					for(i=0;i<m;i++){
						if(idxm[i]<0) continue;
						for(j=0;j<n;j++){
							if(idxn[j]<0) continue;
							this->matrix[N*idxm[i]+idxn[j]]+=values[n*i+j];
						}
					}
					break;
				case INS_VAL:
					for(i=0;i<m;i++){
						if(idxm[i]<0) continue;
						for(j=0;j<n;j++){
							if(idxn[j]<0) continue;
							this->matrix[N*idxm[i]+idxn[j]]=values[n*i+j];
						}
					}
					break;
				default:
					_error_("unknown insert mode!");
					break;
			}

		}
		/*}}}*/
		void Convert(MatrixType type){/*{{{*/

			/*do nothing*/

		}
		/*}}}*/		
		void SetZero(void){/*{{{*/
			for(int i=0;i<M;i++){
				for(int j=0;j<N;j++){
					this->matrix[N*i+j] = 0.;
				}
			}
		}/*}}}*/
		#ifndef _HAVE_WRAPPERS_
		IssmAbsVec<IssmDouble>* Solve(IssmAbsVec<IssmDouble>* pfin, Parameters* parameters){/*{{{*/

			/*First off, we assume that the type of IssmAbsVec is IssmSeqVec. So downcast: */
			IssmSeqVec<IssmDouble>* pf = NULL;
			IssmSeqVec<IssmDouble> *uf = NULL;
			IssmDouble* x=NULL;

			/*Assume we are getting an IssmMpiVec in input, downcast: */
			pf=(IssmSeqVec<IssmDouble>*)pfin;

			switch(IssmSolverTypeFromToolkitOptions()){
			#ifdef _HAVE_MUMPS_
				case MumpsEnum: {
					/*Assume we have a sequential vec, downcast*/
					uf=((IssmSeqVec<IssmDouble>*)pfin)->Duplicate();
					SeqDenseMumpsSolve(uf->vector,uf->M,uf->M, /*stiffness matrix:*/ this->matrix,this->M,this->N,this->M, /*right hand side load vector: */ pf->vector,pf->M,pf->M,parameters);
					return uf;
									 }
			#endif
			#ifdef _HAVE_GSL_
				case GslEnum: {
					DenseGslSolve(/*output*/ &x,/*stiffness matrix:*/ this->matrix,this->M,this->N, /*right hand side load vector: */ pf->vector,pf->M,parameters);
					uf=new IssmSeqVec<IssmDouble>(x,this->N); xDelete(x);
					return uf;
								  }
			#endif
				default: _error_("No solver available");
			}

			return NULL;

		}/*}}}*/
		#endif
};

#endif //#ifndef _ISSM_DENSE_MAT_H_
