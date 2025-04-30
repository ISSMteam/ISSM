/*!\file:  IssmMpiSparseMat.h
 * \brief implementation of parallel sparse ISSM matrix. Internally, the parallel sparse matrix is 
 * split in rows across each cpu. Locally, on each cpu, the local matrix is represented by a vector of sparse rows.
 * This object needs to answer the API defined by the virtual functions in IssmAbsMat, 
 * and the contructors required by IssmMat (see IssmMat.h)
 */ 

#ifndef _ISSM_MPI_SPARSE_MAT_H_
#define _ISSM_MPI_SPARSE_MAT_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../datastructures/datastructures.h"
#include "../../shared/shared.h"
#include "../mumps/mumpsincludes.h"
#include "./Bucket.h"
#include "./IssmMpiVec.h"
#include "./SparseRow.h"
#include <math.h>

#ifdef _HAVE_CODIPACK_
#include "../codipack/CoDiPackDebug.h"
#endif
/*}}}*/

/*We need to template this class, in case we want to create Matrices that hold
  IssmDouble* matrix or IssmPDouble* matrix. 
  Such matrices would be useful for use without or with the matlab or python
  interface (which do not care for IssmDouble types, but only rely on
  IssmPDouble types)*/

template <class doubletype> class IssmAbsMat;

template <class doubletype> 
class IssmMpiSparseMat:public IssmAbsMat<doubletype>{

	public:

		int M,N;  //global size
		int m;    //local number of rows
		SparseRow<doubletype>** matrix;  /*here, doubletype is either IssmDouble or IssmPDouble*/
		DataSet*    buckets;  /*here, we store buckets of values that we will Assemble into a global matrix.*/
		/*IssmMpiSparseMat constructors, destructors*/
		IssmMpiSparseMat(){/*{{{*/
			this->M=0;
			this->N=0;
			this->m=0;
			this->matrix=NULL;
			this->buckets=new DataSet();
		}
		/*}}}*/
		IssmMpiSparseMat(int Min,int Nin){/*{{{*/
			this->Init(Min,Nin);
		}
		/*}}}*/
		IssmMpiSparseMat(int pM,int pN, doubletype sparsity){/*{{{*/
			/*no sparsity involved here, the sparsity pattern is resolve during the assemble phase: */
			this->Init(pM,pN);
		}
		/*}}}*/
		IssmMpiSparseMat(int min,int nin,int Min,int Nin,int* d_nnz,int* o_nnz){/*{{{*/

			int i;

			/*no sparsity involved here, the sparsity pattern is resolved at the assemble phase: */
			this->buckets=new DataSet();

			this->M=Min;
			this->N=Nin;
			this->m=min;

			/*Initialize pointer: */
			this->matrix=NULL;

			/*Allocate: */
			if (m*N){
				this->matrix=xNew<SparseRow<doubletype>*>(m);
				for(i=0;i<m;i++){
					this->matrix[i]=new SparseRow<doubletype>(N);
				}
			}
		}
		/*}}}*/
		IssmMpiSparseMat(int pM,int pN, int connectivity,int numberofdofspernode){/*{{{*/
			/*this is not needed, sparsity pattern is resolved at assemble phase: */
			this->Init(pM,pN);
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
			if (m*N){
				this->matrix=xNew<SparseRow<doubletype>*>(m);
				for(int i=0;i<m;i++){
					this->matrix[i]=new SparseRow<doubletype>(N);
				}
			}
		}
		/*}}}*/
		~IssmMpiSparseMat(){/*{{{*/
			if(m*N){
				for(int i=0;i<m;i++){
					delete this->matrix[i];
				}
				xDelete<SparseRow<doubletype>*>(this->matrix);
			}
			this->M=0;
			this->N=0;
			this->m=0;
			delete this->buckets;
		}
		/*}}}*/

		/*IssmMpiSparseMat specific routines */
		void Echo(void){/*{{{*/

			/*Do a synchronized dump across all the rows: */
			int my_rank=IssmComm::GetRank();
			for(int i=0;i<IssmComm::GetSize();i++){
				if(my_rank==i){
					_printf_("cpu " << i << " #rows: " << this->m << "\n");
					for(int j=0;j<this->m;j++){
						_printf_("row " << j << ":");
						this->matrix[j]->Echo();
						_printf_("\n");
					}
				}
				ISSM_MPI_Barrier(IssmComm::GetComm());
			}

		}/*}}}*/
		void EchoDebug(std::string message) {/*{{{*/
#if defined(_HAVE_CODIPACK_)
			if (std::is_same<doubletype, IssmDouble>::value) {
				void* h = MatDebugOutputStart(message, M, N);

				if(NULL != matrix) {
					for(int local_row = 0; local_row < m; local_row += 1) {
						SparseRow<doubletype>* cur_row = matrix[local_row];
						MatDebugOutputAddRow(h, local_row /*TODO: get global row*/, cur_row->ncols, cur_row->indices, cur_row->values);

					}
				}
				MatDebugOutputFinish(h);
			}
#endif
		}/*}}}*/
		void Assemble(){/*{{{*/

			int         *RowRank = NULL;
			int         *row_indices_forcpu  = NULL;
			int         *col_indices_forcpu  = NULL;
			int         *modes_forcpu        = NULL;
			doubletype  *values_forcpu       = NULL;
			int         *numvalues_forcpu    = NULL;
			DataSet    **bucketsforcpu       = NULL;
			int        **row_indices_fromcpu = NULL;
			int        **col_indices_fromcpu = NULL;
			int        **modes_fromcpu       = NULL;
			doubletype **values_fromcpu      = NULL;
			int         *numvalues_fromcpu   = NULL;
			int          lower_row;
			int          upper_row;
			int         *sendcnts = NULL;
			int         *displs   = NULL;
			int          this_row_numvalues;
			int         *this_row_cols       = NULL;
			int         *this_row_mods       = NULL;
			int         *numvalues_perrow    = NULL;

			doubletype **values_perrow       = NULL;
			int        **cols_perrow         = NULL;
			int        **mods_perrow         = NULL;
			int         *counters_perrow     = NULL;

			/*Early exit: */
			if(this->M*this->N==0) return;

			/*some communicator info: */
			int num_procs=IssmComm::GetSize();
			ISSM_MPI_Comm comm=IssmComm::GetComm();

			/*First, make a vector of size M, which for each row between 0 and M-1, tells which cpu this row belongs to: */
			RowRank=DetermineRowRankFromLocalSize(M,m,comm);

			/*Now, sort out our dataset of buckets according to cpu ownership of rows*/
			bucketsforcpu=xNew<DataSet*>(num_procs);
			for(int i=0;i<num_procs;i++){
				DataSet* bucketsofcpu_i=new DataSet();
				for(int j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket=(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					bucket->SpawnBucketsPerCpu(bucketsofcpu_i,i,RowRank);
				}
				bucketsforcpu[i]=bucketsofcpu_i;
			}

			/*Recap, each cpu has num_procs datasets of buckets. For a certain cpu j, for a given dataset i, the buckets this 
			 * dataset owns correspond to rows that are owned by cpu i, not j!. Out of all the buckets we own, make row,col,value,insert_mode 
			 * vectors that will be shipped around the cluster: */
			this->BucketsBuildScatterBuffers(&numvalues_forcpu,&row_indices_forcpu,&col_indices_forcpu,&values_forcpu,&modes_forcpu,bucketsforcpu,num_procs);

			/*Now, we need to allocate on each cpu arrays to receive data from all the other cpus. To know what we need to allocate, we need 
			 *some scatter calls: */
			numvalues_fromcpu   = xNew<int>(num_procs);
			for(int i=0;i<num_procs;i++) ISSM_MPI_Scatter(numvalues_forcpu,1,ISSM_MPI_INT,numvalues_fromcpu+i,1,ISSM_MPI_INT,i,comm);

			row_indices_fromcpu=xNew<int*>(num_procs);
			col_indices_fromcpu=xNew<int*>(num_procs);
			values_fromcpu=xNew<doubletype*>(num_procs);
			modes_fromcpu=xNew<int*>(num_procs);
			for(int i=0;i<num_procs;i++){
				int size=numvalues_fromcpu[i];
				if(size){
					row_indices_fromcpu[i]=xNew<int>(size);
					col_indices_fromcpu[i]=xNew<int>(size);
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
					values_fromcpu[i]=xNew<doubletype>(size,"t");
#else
					values_fromcpu[i]=xNew<doubletype>(size);
#endif
					modes_fromcpu[i]=xNew<int>(size);
				}
				else{
					row_indices_fromcpu[i] = NULL;
					col_indices_fromcpu[i] = NULL;
					values_fromcpu[i]      = NULL;
					modes_fromcpu[i]       = NULL;
				}
			}

			/*Scatter values around*/
			/*Now, to scatter values across the cluster, we need sendcnts and displs. Our sendbufs have been built by BucketsBuildScatterBuffers, with a stride given 
			 * by numvalues_forcpu. Get this ready to go before starting the scatter itslef. For reference, here is the ISSM_MPI_Scatterv prototype: 
			 * int ISSM_MPI_Scatterv( void *sendbuf, int *sendcnts, int *displs, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm) :*/
			sendcnts=xNew<int>(num_procs);
			displs=xNew<int>(num_procs);
			int count=0;
			for(int i=0;i<num_procs;i++){
				sendcnts[i]=numvalues_forcpu[i];
				displs[i]=count;
				count+=numvalues_forcpu[i];
			}

			for(int i=0;i<num_procs;i++){
				ISSM_MPI_Scatterv( row_indices_forcpu, sendcnts, displs, ISSM_MPI_INT, row_indices_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
				ISSM_MPI_Scatterv( col_indices_forcpu, sendcnts, displs, ISSM_MPI_INT, col_indices_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
				ISSM_MPI_Scatterv( values_forcpu, sendcnts, displs, ISSM_MPI_DOUBLE, values_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_DOUBLE, i, comm);
				ISSM_MPI_Scatterv( modes_forcpu, sendcnts, displs, ISSM_MPI_INT, modes_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
			}

			/*Plug values into global matrix. To do so, we are going to first figure out how many overall values each sparse row is going to get, then we fill up these values, and give it to each sparse row*/
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,comm);

			/*Figure out how many values each row is going to get: */
			numvalues_perrow=xNewZeroInit<int>(this->m);
			for(int i=0;i<num_procs;i++){ 
				int  numvalues=numvalues_fromcpu[i];
				int* rows=row_indices_fromcpu[i];
				for(int j=0;j<numvalues;j++)numvalues_perrow[rows[j]-lower_row]++;
			}

			/*Allocate all the values, cols and mods from each cpu: */
			values_perrow=xNew<doubletype*>(this->m);
			cols_perrow=xNew<int*>(this->m);
			mods_perrow=xNew<int*>(this->m);
			counters_perrow=xNewZeroInit<int>(this->m);

			for(int i=0;i<this->m;i++){
				values_perrow[i]=xNewZeroInit<doubletype>(numvalues_perrow[i]);
				cols_perrow[i]=xNewZeroInit<int>(numvalues_perrow[i]);
				mods_perrow[i]=xNewZeroInit<int>(numvalues_perrow[i]);
			}

			/*collect:*/
			for(int i=0;i<num_procs;i++){
				int  numvalues=numvalues_fromcpu[i];
				int* rows=row_indices_fromcpu[i];
				int* cols=col_indices_fromcpu[i];
				doubletype* values=values_fromcpu[i];
				int* mods=modes_fromcpu[i];
				for(int j=0;j<numvalues;j++){
					int row=rows[j]-lower_row;
					int counter=counters_perrow[row];
					values_perrow[row][counter]=values[j];
					cols_perrow[row][counter]=cols[j];
					mods_perrow[row][counter]=mods[j];
					counter=counters_perrow[row]++;
				}
			}

			/*Plug into matrix: */
			for(int i=0;i<this->m;i++) this->matrix[i]->SetValues(numvalues_perrow[i],cols_perrow[i],values_perrow[i],mods_perrow[i]);

			/*Free resources:*/
			xDelete<int>(numvalues_perrow);
			xDelete<int>(RowRank);
			xDelete<int>(row_indices_forcpu);
			xDelete<int>(col_indices_forcpu);
			xDelete<int>(modes_forcpu);
			xDelete<doubletype>(values_forcpu);
			xDelete<int>(numvalues_forcpu);
			for(int i=0;i<num_procs;i++){
				DataSet* buckets=bucketsforcpu[i];
				delete buckets;
			}
			xDelete<DataSet*>(bucketsforcpu);
			for(int i=0;i<num_procs;i++){
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
			for(int i=0;i<this->m;i++){
				doubletype* values=values_perrow[i]; xDelete<doubletype>(values);
				int* cols=cols_perrow[i]; xDelete<int>(cols);
				int* mods=mods_perrow[i]; xDelete<int>(mods);
			}
			xDelete<int>(counters_perrow);
			xDelete<doubletype*>(values_perrow);
			xDelete<int*>(cols_perrow);
			xDelete<int*>(mods_perrow);
		}/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/

			doubletype norm,local_norm;
			doubletype absolute;
			int i;

			switch(mode){
				case NORM_INF:
					local_norm=0;
					for(i=0;i<this->m;i++){
						local_norm=max(local_norm,this->matrix[i]->Norm(mode));
					}
					ISSM_MPI_Reduce(&local_norm, &norm, 1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, 0, IssmComm::GetComm());
					ISSM_MPI_Bcast(&norm,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
					return norm;
					break; 
				case NORM_FROB:
					local_norm=0;
					for(i=0;i<this->m;i++){
						local_norm+=this->matrix[i]->Norm(mode);
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

			/*A check on the types: */
			if(IssmVecTypeFromToolkitOptions()!=MpiEnum)_error_("MatMult operation only possible with 'mpi' vectors");

			/*Now that we are sure, cast vectors: */
			IssmMpiVec<doubletype>* X=(IssmMpiVec<doubletype>*)Xin;
			IssmMpiVec<doubletype>* AX=(IssmMpiVec<doubletype>*)AXin;

			/*Serialize input Xin: */
			doubletype* X_serial=X->ToMPISerial();

			/*Every cpu has a serial version of the input vector. Use it to do the Matrix-Vector multiply 
			 *locally and plug it into AXin: */
			for(int i=0;i<this->m;i++){
				AX->vector[i]=this->matrix[i]->Mult(X_serial);
			}

			/*Free resources: */
			xDelete<doubletype>(X_serial);
		}
		/*}}}*/
		IssmMpiSparseMat<doubletype>* Duplicate(void){/*{{{*/

			_error_("not supported yet!");

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

			/*Reset buckets*/
			delete this->buckets;
			this->buckets=new DataSet();

			/*reset matrix*/
			if(m*N){
				for(int i=0;i<m;i++) delete this->matrix[i];
				xDelete<SparseRow<doubletype>*>(this->matrix);

				this->matrix=xNew<SparseRow<doubletype>*>(m);
				for(int i=0;i<m;i++) this->matrix[i]=new SparseRow<doubletype>(N);
			}

			/*Reallocate matrix*/
		}/*}}}*/
		void Convert(MatrixType type){/*{{{*/
			_error_("not supported yet!");
		}
		/*}}}*/		
		void BucketsBuildScatterBuffers(int** pnumvalues_forcpu,int** prow_indices_forcpu,int** pcol_indices_forcpu,doubletype** pvalues_forcpu,int** pmodes_forcpu,DataSet** bucketsforcpu,int num_procs){/*{{{*/

			/*intermediary: */
			int         count;
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
			for(int i=0;i<num_procs;i++){
				DataSet    *buckets = bucketsforcpu[i];
				count=0;
				for(int j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket =(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					count+=bucket->MarshallSize();
				}
				numvalues_forcpu[i]=count;
			}

			/*now, figure out size of  total buffers (for all cpus!): */
			count=0;
			for(int i=0;i<num_procs;i++) count+=numvalues_forcpu[i];
			int total_size=count;

			/*Allocate buffers: */
			row_indices_forcpu = xNew<int>(total_size);
			col_indices_forcpu = xNew<int>(total_size);
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
			values_forcpu = xNew<doubletype>(total_size,"t");
#else
			values_forcpu = xNew<doubletype>(total_size);
#endif
			modes_forcpu = xNew<int>(total_size);

			/*we are going to march through the buffers, and marshall data onto them, so in order to not
			 *lose track of where these buffers are located in memory, we are going to work using copies 
			 of them: */
			temp_row_indices_forcpu = row_indices_forcpu;
			temp_col_indices_forcpu = col_indices_forcpu;
			temp_values_forcpu      = values_forcpu;
			temp_modes_forcpu       = modes_forcpu;

			/*Fill buffers: */
			for(int i=0;i<num_procs;i++){
				DataSet *buckets = bucketsforcpu[i];
				for(int j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket =(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					bucket->Marshall(&temp_row_indices_forcpu,&temp_col_indices_forcpu,&temp_values_forcpu,&temp_modes_forcpu); //pass in the address of the buffers, so as to have the Marshall routine increment them.
				}
			}

			/*sanity check: */
			if(temp_row_indices_forcpu!=row_indices_forcpu+total_size)_error_("problem with marshalling of buckets");
			if(temp_col_indices_forcpu!=col_indices_forcpu+total_size)_error_("problem with marshalling of buckets");
			if(temp_values_forcpu!=values_forcpu+total_size)_error_("problem with marshalling of buckets");
			if(temp_modes_forcpu!=modes_forcpu+total_size)_error_("problem with marshalling of buckets");

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
					MpiSparseMumpsSolve(/*output*/ uf->vector,uf->M,uf->m, /*stiffness matrix:*/ this->matrix,this->M,this->N,this->m, /*right hand side load vector: */ pf->vector,pf->M,pf->m,parameters);
					#else
					_error_("IssmMpiSparseMat solver requires MUMPS solver");
					#endif
					return (IssmAbsVec<IssmDouble>*)uf;
									 }
				case GslEnum: {

					IssmDouble* x=NULL;
					#ifdef _HAVE_GSL_

					_error_("not implemented yet!");
					//SparseGslSolve(/*output*/ &x,/*stiffness matrix:*/ this->matrix,this->M,this->N, /*right hand side load vector: */ pf->vector,pf->M,parameters);

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

#endif //#ifndef _ISSM_MPI_SPARSE_MAT_H_
