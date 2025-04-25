/*!\file:  IssmMpiDenseMat.h
 * \brief implementation of parallel dense ISSM matrix. Internally, the parallel dense matrix is 
 * split in rows across each cpu. Each matrix (representing a subset of rows) on each cpu is fully 
 * dense, and is represented by a linear buffer of type doubletype. 
 * This object needs to answer the API defined by the virtual functions in IssmAbsMat, 
 * and the contructors required by IssmMat (see IssmMat.h)
 */ 

#ifndef _ISSM_MPI_DENSE_MAT_H_
#define _ISSM_MPI_DENSE_MAT_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/shared.h"
#include "../../datastructures/datastructures.h"
#include "../mumps/mumpsincludes.h"
#include "./Bucket.h"
#include "./IssmMpiVec.h"
#include <math.h>

/*}}}*/

/*We need to template this class, in case we want to create Matrices that hold
  IssmDouble* matrix or IssmPDouble* matrix. 
  Such matrices would be useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)*/
template <class doubletype> class IssmAbsMat;

template <class doubletype> 
class IssmMpiDenseMat:public IssmAbsMat<doubletype>{

	public:

		int M,N;  //global size
		int m;    //local number of rows
		doubletype* matrix;  /*here, doubletype is either IssmDouble or IssmPDouble*/
		DataSet*    buckets;  /*here, we store buckets of values that we will Assemble into a global matrix.*/

		/*IssmMpiDenseMat constructors, destructors*/
		IssmMpiDenseMat(){/*{{{*/
			this->M=0;
			this->N=0;
			this->m=0;
			this->matrix=NULL;
			this->buckets=new DataSet();
		}
		/*}}}*/
		IssmMpiDenseMat(int Min,int Nin){/*{{{*/
			this->Init(Min,Nin);
		}
		/*}}}*/
		IssmMpiDenseMat(int pM,int pN, doubletype sparsity){/*{{{*/
			/*no sparsity involved here, we are fully dense, so just use the previous constructor: */
			this->Init(pM,pN);
		}
		/*}}}*/
		IssmMpiDenseMat(int min,int nin,int Min,int Nin,int* d_nnz,int* o_nnz){/*{{{*/
			/*not needed, we are fully dense!: */

			this->buckets=new DataSet();

			this->M=Min;
			this->N=Nin;
			this->m=min;

			/*Initialize pointer: */
			this->matrix=NULL;

			/*Allocate: */
			if (m*N)this->matrix=xNewZeroInit<doubletype>(this->m*N);
		}
		/*}}}*/
		IssmMpiDenseMat(doubletype* serial_mat,int Min,int Nin,doubletype sparsity){/*{{{*/

			/*Here, we assume that the serial_mat is local to the local cpu, and that it has 
			 * the correct size (m rows by N colums), n determined by DetermineLocalSize: */
			this->buckets=new DataSet();
			this->M=Min;
			this->N=Nin;
			this->m=DetermineLocalSize(this->M,IssmComm::GetComm());

			this->matrix=NULL;
			if(m*N){
				this->matrix=xNewZeroInit<doubletype>(m*N);
				xMemCpy<doubletype>(this->matrix,serial_mat,m*N);
			}
		}
		/*}}}*/
		IssmMpiDenseMat(int pM,int pN, int connectivity,int numberofdofspernode){/*{{{*/
			/*not needed, we are fully dense!: */
			this->Init(pM,pN);
		}
		/*}}}*/
		~IssmMpiDenseMat(){/*{{{*/
			xDelete<doubletype>(this->matrix);
			M=0;
			N=0;
			m=0;
			delete this->buckets;
		}
		/*}}}*/
		void Init(int Min,int Nin){/*{{{*/

			this->buckets=new DataSet();

			this->M=Min;
			this->N=Nin;

			/*Figure out local number of rows: */
			this->m=DetermineLocalSize(this->M,IssmComm::GetComm());

			/*Initialize pointer: */
			this->matrix=NULL;

			/*Allocate: */
			if (m*N)this->matrix=xNewZeroInit<doubletype>(this->m*N);
		}
		/*}}}*/

		/*IssmMpiDenseMat specific routines */
		void Echo(void){/*{{{*/

			int my_rank;
			int i,j,k;

			/*Do a synchronized dump across all the rows: */
			my_rank=IssmComm::GetRank();
			for(i=0;i<IssmComm::GetSize();i++){
				if (my_rank==i){
					_printf_("cpu " << i << " #rows: " << this->m << "\n");
					for (j=0;j<this->m;j++){
						_printf_("row " << j << ":");
						for (k=0;k<this->N;k++){
							if(this->matrix[j*this->N+k]!=0)_printf_("(" << k << "," << this->matrix[j*this->N+k] << ") ");
						}
						_printf_("\n");
					}
				}
				ISSM_MPI_Barrier(IssmComm::GetComm());
			}

		}
		/*}}}*/
		void EchoDebug(std::string message) {
			_printf_("Error: not implemented.\n");
		}

		void Assemble(){/*{{{*/

			int           i,j;

			int         *RowRank            = NULL;
			int           num_procs;

			int        *row_indices_forcpu = NULL;
			int        *col_indices_forcpu = NULL;
			int        *modes_forcpu       = NULL;
			doubletype *values_forcpu      = NULL;
			int         *numvalues_forcpu   = NULL;
			DataSet     **bucketsforcpu       = NULL;

			int        **row_indices_fromcpu = NULL;
			int        **col_indices_fromcpu = NULL;
			int        **modes_fromcpu       = NULL;
			doubletype **values_fromcpu      = NULL;
			int         *numvalues_fromcpu   = NULL;

			int           lower_row;
			int           upper_row;
			int*          sendcnts            = NULL;
			int*          displs              = NULL;
			int           count               = 0;

			/*some communicator info: */
			num_procs=IssmComm::GetSize();
			ISSM_MPI_Comm comm=IssmComm::GetComm();

			/*First, make a vector of size M, which for each row between 0 and M-1, tells which cpu this row belongs to: */
			RowRank=DetermineRowRankFromLocalSize(M,m,comm);

			/*Now, sort out our dataset of buckets according to cpu ownership of rows: {{{*/
			bucketsforcpu=xNew<DataSet*>(num_procs);

			for(i=0;i<num_procs;i++){
				DataSet* bucketsofcpu_i=new DataSet();
				for (j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket=(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					bucket->SpawnBucketsPerCpu(bucketsofcpu_i,i,RowRank);
				}
				bucketsforcpu[i]=bucketsofcpu_i;
			}
			/*}}}*/

			/*Recap, each cpu has num_procs datasets of buckets. For a certain cpu j, for a given dataset i, the buckets this  {{{
			 * dataset owns correspond to rows that are owned by cpu i, not j!. Out of all the buckets we own, make row,col,value,insert_mode 
			 * vectors that will be shipped around the cluster: */
			this->BucketsBuildScatterBuffers(&numvalues_forcpu,&row_indices_forcpu,&col_indices_forcpu,&values_forcpu,&modes_forcpu,bucketsforcpu,num_procs);
			/*}}}*/

			/*Now, we need to allocate on each cpu arrays to receive data from all the other cpus. To know what we need to allocate, we need  {{{
			 *some scatter calls: */
			numvalues_fromcpu   = xNew<int>(num_procs);
			for(i=0;i<num_procs;i++){
				ISSM_MPI_Scatter(numvalues_forcpu,1,ISSM_MPI_INT,numvalues_fromcpu+i,1,ISSM_MPI_INT,i,comm);
			}

			row_indices_fromcpu=xNew<int*>(num_procs);
			col_indices_fromcpu=xNew<int*>(num_procs);
			values_fromcpu=xNew<doubletype*>(num_procs);
			modes_fromcpu=xNew<int*>(num_procs);
			for(i=0;i<num_procs;i++){
				int size=numvalues_fromcpu[i];
				if(size){
					row_indices_fromcpu[i]=xNew<int>(size);
					col_indices_fromcpu[i]=xNew<int>(size);
					values_fromcpu[i]=xNew<doubletype>(size);
					modes_fromcpu[i]=xNew<int>(size);
				}
				else{
					row_indices_fromcpu[i]=NULL;
					col_indices_fromcpu[i]=NULL;
					values_fromcpu[i]=NULL;
					modes_fromcpu[i]=NULL;
				}
			}
			/*}}}*/

			/*Scatter values around: {{{*/
			/*Now, to scatter values across the cluster, we need sendcnts and displs. Our sendbufs have been built by BucketsBuildScatterBuffers, with a stride given 
			 * by numvalues_forcpu. Get this ready to go before starting the scatter itslef. For reference, here is the ISSM_MPI_Scatterv prototype: 
			 * int ISSM_MPI_Scatterv( void *sendbuf, int *sendcnts, int *displs, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm) :*/
			sendcnts=xNew<int>(num_procs);
			displs=xNew<int>(num_procs);
			count=0;
			for(i=0;i<num_procs;i++){
				sendcnts[i]=numvalues_forcpu[i];
				displs[i]=count;
				count+=numvalues_forcpu[i];
			}

			for(i=0;i<num_procs;i++){
				ISSM_MPI_Scatterv( row_indices_forcpu, sendcnts, displs, ISSM_MPI_INT, row_indices_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
				ISSM_MPI_Scatterv( col_indices_forcpu, sendcnts, displs, ISSM_MPI_INT, col_indices_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
				ISSM_MPI_Scatterv( values_forcpu, sendcnts, displs, ISSM_MPI_DOUBLE, values_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_DOUBLE, i, comm);
				ISSM_MPI_Scatterv( modes_forcpu, sendcnts, displs, ISSM_MPI_INT, modes_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
			}
			/*}}}*/

			/*Plug values into global matrix: {{{*/
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,comm);
			for(i=0;i<num_procs;i++){
				int  numvalues=numvalues_fromcpu[i];
				int* rows=row_indices_fromcpu[i];
				int* cols=col_indices_fromcpu[i];
				doubletype* values=values_fromcpu[i];
				int* mods=modes_fromcpu[i];

				for(j=0;j<numvalues;j++){
					if(mods[j]==ADD_VAL) *(matrix+N*(rows[j]-lower_row)+cols[j])+=values[j];
					else *(matrix+N*(rows[j]-lower_row)+cols[j])=values[j];
				}
			}
			/*}}}*/

			/*Free resources:{{{*/
			xDelete<int>(RowRank);
			xDelete<int>(row_indices_forcpu);
			xDelete<int>(col_indices_forcpu);
			xDelete<int>(modes_forcpu);
			xDelete<doubletype>(values_forcpu);
			xDelete<int>(numvalues_forcpu);

			for(i=0;i<num_procs;i++){
				DataSet* buckets=bucketsforcpu[i];
				delete buckets;
			}
			xDelete<DataSet*>(bucketsforcpu);

			for(i=0;i<num_procs;i++){
				int* rows=row_indices_fromcpu[i];
				int* cols=col_indices_fromcpu[i];
				int* modes=modes_fromcpu[i];
				doubletype* values=values_fromcpu[i];

				xDelete<int>(rows);
				xDelete<int>(cols);
				xDelete<int>(modes);
				xDelete<doubletype>(values);
			}
			xDelete<int*>(row_indices_fromcpu);
			xDelete<int*>(col_indices_fromcpu);
			xDelete<int*>(modes_fromcpu);
			xDelete<doubletype*>(values_fromcpu);
			xDelete<int>(numvalues_fromcpu);

			xDelete<int>(sendcnts);
			xDelete<int>(displs);
			/*}}}*/

		}
		/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/

			doubletype norm,local_norm;
			doubletype absolute;
			int i,j;

			switch(mode){
				case NORM_INF:
					local_norm=0.;
					for(i=0;i<this->m;i++){
						absolute=0;
						for(j=0;j<this->N;j++){
							absolute+=fabs(this->matrix[N*i+j]);
						}
						local_norm=max(local_norm,absolute);
					}
					ISSM_MPI_Reduce(&local_norm, &norm, 1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, 0, IssmComm::GetComm());
					ISSM_MPI_Bcast(&norm,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
					return norm;
					break; 
				case NORM_FROB:
					local_norm=0.;
					for(i=0;i<this->m;i++){
						for(j=0;j<this->N;j++){
							local_norm+=this->matrix[N*i+j]*this->matrix[N*i+j];
						}
					}
					ISSM_MPI_Reduce(&local_norm, &norm, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
					ISSM_MPI_Bcast(&norm,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
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
			*pM=M;
			*pN=N;
		}
		/*}}}*/
		void GetLocalSize(int* pM,int* pN){/*{{{*/
			*pM=m;
			*pN=N;
		}
		/*}}}*/
		void MatMult(IssmAbsVec<doubletype>* Xin,IssmAbsVec<doubletype>* AXin){/*{{{*/

			int         i,j;
			doubletype *X_serial  = NULL;

			/*A check on the types: */
			if(IssmVecTypeFromToolkitOptions()!=MpiEnum)_error_("MatMult operation only possible with 'mpi' vectors");

			/*Now that we are sure, cast vectors: */
			IssmMpiVec<doubletype>* X=(IssmMpiVec<doubletype>*)Xin;
			IssmMpiVec<doubletype>* AX=(IssmMpiVec<doubletype>*)AXin;

			/*Serialize input Xin: */
			X_serial=X->ToMPISerial();

			/*Every cpu has a serial version of the input vector. Use it to do the Matrix-Vector multiply 
			 *locally and plug it into AXin: */
			for(i=0;i<this->m;i++){
				for(j=0;j<this->N;j++){
					AX->vector[i]+=this->matrix[i*N+j]*X_serial[j];
				}
			}

			/*Free resources: */
			xDelete<doubletype>(X_serial);
		}
		/*}}}*/
		IssmMpiDenseMat<doubletype>* Duplicate(void){/*{{{*/

			IssmMpiDenseMat<doubletype>* dup=new IssmMpiDenseMat<doubletype>(this->matrix,this->M,this->N,0);
			return dup;

		}
		/*}}}*/
		doubletype* ToSerial(void){/*{{{*/
			_error_("not supported yet!");
		}
		/*}}}*/
		void SetValues(int min,int* idxm,int nin,int* idxn,doubletype* values,InsMode mode){/*{{{*/

			/*we need to store all the values we collect here in order to Assemble later. 
			 * Indeed, the values we are collecting here most of the time will not belong 
			 * to us, but to another part of the matrix on another cpu: */
			_assert_(buckets);

			buckets->AddObject(new Bucket<doubletype>(min,idxm,nin,idxn,values,mode));

		}
		/*}}}*/
		void SetZero(void){/*{{{*/
			for(int i=0;i<this->m*this->N;i++) this->matrix[i] = 0.;
		}/*}}}*/
		void Convert(MatrixType type){/*{{{*/
			_error_("not supported yet!");
		}
		/*}}}*/		
		void BucketsBuildScatterBuffers(int** pnumvalues_forcpu,int** prow_indices_forcpu,int** pcol_indices_forcpu,doubletype** pvalues_forcpu,int** pmodes_forcpu,DataSet** bucketsforcpu,int num_procs){/*{{{*/

			/*intermediary: */
			int         i,j;
			int         count                   = 0;
			int         total_size              = 0;
			int        *temp_row_indices_forcpu = NULL;
			int        *temp_col_indices_forcpu = NULL;
			doubletype *temp_values_forcpu      = NULL;
			int        *temp_modes_forcpu       = NULL;

			/*output: */
			int        *numvalues_forcpu        = NULL;
			int        *row_indices_forcpu      = NULL;
			int        *col_indices_forcpu      = NULL;
			doubletype *values_forcpu           = NULL;
			int        *modes_forcpu            = NULL;

			/*figure out size of buffers per cpu: */

			numvalues_forcpu=xNew<int>(num_procs);
			for(i=0;i<num_procs;i++){
				DataSet    *buckets            = bucketsforcpu[i];

				count=0;
				for(j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket =(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					count+=bucket->MarshallSize();
				}

				numvalues_forcpu[i]=count;
			}

			/*now, figure out size of  total buffers (for all cpus!): */
			count=0;
			for(i=0;i<num_procs;i++){
				count+=numvalues_forcpu[i];
			}
			total_size=count;

			/*Allocate buffers: */
			row_indices_forcpu = xNew<int>(total_size);
			col_indices_forcpu = xNew<int>(total_size);
			values_forcpu = xNew<doubletype>(total_size);
			modes_forcpu = xNew<int>(total_size);

			/*we are going to march through the buffers, and marshall data onto them, so in order to not
			 *lose track of where these buffers are located in memory, we are going to work using copies 
			 of them: */
			temp_row_indices_forcpu=row_indices_forcpu;
			temp_col_indices_forcpu=col_indices_forcpu;
			temp_values_forcpu=values_forcpu;
			temp_modes_forcpu=modes_forcpu;

			/*Fill buffers: */
			for(i=0;i<num_procs;i++){
				DataSet    *buckets            = bucketsforcpu[i];
				for(j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket =(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					bucket->Marshall(&temp_row_indices_forcpu,&temp_col_indices_forcpu,&temp_values_forcpu,&temp_modes_forcpu); //pass in the address of the buffers, so as to have the Marshall routine increment them.
				}
			}

			/*sanity check: */
			if (temp_row_indices_forcpu!=row_indices_forcpu+total_size)_error_("problem with marshalling of buckets");
			if (temp_col_indices_forcpu!=col_indices_forcpu+total_size)_error_("problem with marshalling of buckets");
			if (temp_values_forcpu!=values_forcpu+total_size)_error_("problem with marshalling of buckets");
			if (temp_modes_forcpu!=modes_forcpu+total_size)_error_("problem with marshalling of buckets");

			/*output buffers: */
			*pnumvalues_forcpu   = numvalues_forcpu;
			*prow_indices_forcpu = row_indices_forcpu;
			*pcol_indices_forcpu = col_indices_forcpu;
			*pvalues_forcpu      = values_forcpu;
			*pmodes_forcpu       = modes_forcpu;
		}
		/*}}}*/		
		#ifndef _HAVE_WRAPPERS_
		/*Solve{{{*/
		IssmAbsVec<IssmDouble>* Solve(IssmAbsVec<IssmDouble>* pfin, Parameters* parameters){

			/*output: */
			IssmMpiVec<IssmDouble>* uf=NULL;
			IssmMpiVec<IssmDouble>* pf=NULL;

			/*Assume we are getting an IssmMpiVec in input, downcast: */
			pf=(IssmMpiVec<IssmDouble>*)pfin;

			switch(IssmSolverTypeFromToolkitOptions()){
				case MumpsEnum: {
					/*Initialize output: */
					uf=pf->Duplicate();
					#ifdef _HAVE_MUMPS_
					MpiDenseMumpsSolve(/*output*/ uf->vector,uf->M,uf->m, /*stiffness matrix:*/ this->matrix,this->M,this->N,this->m, /*right hand side load vector: */ pf->vector,pf->M,pf->m,parameters);
					#else
					_error_("IssmMpiDenseMat solver requires MUMPS solver");
					#endif
					return (IssmAbsVec<IssmDouble>*)uf;
									 }
				case GslEnum: {

					IssmDouble* x=NULL;
					#ifdef _HAVE_GSL_

					DenseGslSolve(/*output*/ &x,/*stiffness matrix:*/ this->matrix,this->M,this->N, /*right hand side load vector: */ pf->vector,pf->M,parameters);

					uf=new IssmMpiVec<IssmDouble>(x,this->N); xDelete(x);

					return (IssmAbsVec<IssmDouble>*)uf;
					#else
					_error_("GSL support not compiled in!");
					#endif
								  }
				default:
					_error_("solver type not supported yet!");
			}

		}/*}}}*/
		#endif
};

#endif //#ifndef _ISSM_MPI_DENSE_MAT_H_
