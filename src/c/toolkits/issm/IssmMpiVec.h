/*!\file:  IssmMpiVec.h
 * \brief implementation of parallel dense ISSM vector. Internally, the parallel dense vector is
 * split in rows across each cpu. Each vector (representing a subset of rows) on each cpu is fully
 * dense, and is represented by a linear buffer of type doubletype.
 * This object needs to answer the API defined by the virtual functions in IssmAbsVec,
 * and the contructors required by IssmVec (see IssmVec.h)
 */

#ifndef _ISSM_MPI_VEC_H_
#define _ISSM_MPI_VEC_H_

/*Headers:*/
/*{{{*/
#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../../shared/Exceptions/exceptions.h"
#include "../../shared/MemOps/MemOps.h"
#include "../../shared/io/io.h"
#include "../mpi/issmmpi.h"
#include <math.h>

#ifdef _HAVE_CODIPACK_
#include "../codipack/CoDiPackDebug.h"
#endif

/*}}}*/

/*We need to template this class, in case we want to create vectors that hold IssmDouble* vector or IssmPDouble* vector.
  Such vectors would be useful for use without or with the matlab or python interface (which do not care for IssmDouble types,
  but only rely on IssmPDouble types)*/
template <class doubletype> class IssmAbsVec;

template <class doubletype>
class IssmMpiVec:public IssmAbsVec<doubletype>{

	public:

		bool isassembled;
		int M; //global size
		int m; //local number of rows
		doubletype* vector;  /*here, doubletype is either IssmDouble or IssmPDouble*/
		DataSet*    buckets;  /*here, we store buckets of values that we will Assemble into a global vector.*/

		/*IssmMpiVec constructors, destructors*/
		IssmMpiVec(){/*{{{*/

			this->M=0;
			this->m=0;
			this->isassembled=false;
			this->vector=NULL;
			this->buckets=new DataSet();
		}
		/*}}}*/
		IssmMpiVec(int Min){/*{{{*/
			this->Init(Min,false);
		}
		/*}}}*/
		IssmMpiVec(int min, int Min){/*{{{*/
			this->Init(min,true);
		}
		/*}}}*/
		IssmMpiVec(int Min, bool fromlocalsize){/*{{{*/
			this->Init(Min,fromlocalsize);
		}
		/*}}}*/
		IssmMpiVec(doubletype* buffer,int Min){/*{{{*/

			this->Init(Min,false);

			if(this->M){
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
				this->vector=xNew<doubletype>(this->m,"t");
#else
				this->vector=xNew<doubletype>(this->m);
#endif

				xMemCpy<doubletype>(this->vector,buffer,this->m);
			}
		}
		/*}}}*/
		IssmMpiVec(doubletype* buffer,int Min,int min){/*{{{*/

			this->vector=NULL;
			this->buckets=new DataSet();
			this->M=Min;
			this->m=min;

			if(this->m){
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
				this->vector=xNew<doubletype>(this->m,"t");
#else
				this->vector=xNew<doubletype>(this->m);
#endif

				xMemCpy<doubletype>(this->vector,buffer,this->m);
			}
		}
		/*}}}*/
		void Init(int Min,bool fromlocalsize){/*{{{*/

			this->buckets=new DataSet();
			this->isassembled=false;

			if(fromlocalsize){
				this->m=Min;
				this->M=DetermineGlobalSize(this->m,IssmComm::GetComm());
			}
			else{
				this->M=Min;
				this->m=DetermineLocalSize(this->M,IssmComm::GetComm());
			}

			/*Initialize pointer: */
			this->vector=NULL;

			/*Allocate: */
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
			if (m)this->vector=xNewZeroInit<doubletype>(this->m,"t");
#else
			if (m)this->vector=xNewZeroInit<doubletype>(this->m);
#endif

		}
		/*}}}*/
		~IssmMpiVec(){/*{{{*/
			xDelete<doubletype>(this->vector);
			this->M=0;
			this->m=0;
			delete buckets;
		}
		/*}}}*/

		/*IssmMpiVec specific routines*/
		void Echo(void){/*{{{*/

			/*Do a synchronized dump across all the rows: */
			int my_rank=IssmComm::GetRank();
			for(int i=0;i<IssmComm::GetSize();i++){
				if(my_rank==i){
					if(i==0) _printf_("Vector of global size M="<<this->M<<"\n");
					_printf_("cpu " << i << " #rows: " << this->m << "\n");
					if(this->isassembled){
						for(int j=0;j<this->m;j++){
							_printf_("row " << j << ": "<<this->vector[j]);
							_printf_("\n");
						}
					}
					else{
						this->buckets->DeepEcho();
					}
				}
				ISSM_MPI_Barrier(IssmComm::GetComm());
			}
		}
		/*}}}*/

		void EchoDebug(std::string message){/*{{{*/
#if defined(_HAVE_CODIPACK_)
			if (std::is_same<doubletype, IssmDouble>::value) {
				VecDebugOutput(message, this->M, this->m, this->vector);
			}

#endif
		}/*}}}*/
		void Assemble(){/*{{{*/

			int           i,j;

			int         *RowRank            = NULL;
			int           num_procs;

			int        *row_indices_forcpu = NULL;
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
			this->BucketsBuildScatterBuffers(&numvalues_forcpu,&row_indices_forcpu,&values_forcpu,&modes_forcpu,bucketsforcpu,num_procs);
			/*}}}*/

			/*Now, we need to allocate on each cpu arrays to receive data from all the other cpus. To know what we need to allocate, we need  {{{
			 *some scatter calls: */
			numvalues_fromcpu   = xNew<int>(num_procs);
			for(i=0;i<num_procs;i++){
				ISSM_MPI_Scatter(numvalues_forcpu,1,ISSM_MPI_INT,numvalues_fromcpu+i,1,ISSM_MPI_INT,i,comm);
			}

			row_indices_fromcpu=xNew<int*>(num_procs);
			values_fromcpu=xNew<doubletype*>(num_procs);
			modes_fromcpu=xNew<int*>(num_procs);
			for(i=0;i<num_procs;i++){
				int size=numvalues_fromcpu[i];
				if(size){
					row_indices_fromcpu[i]=xNew<int>(size);
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
					row_indices_fromcpu[i]=NULL;
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
				sendcnts[i] = numvalues_forcpu[i];
				displs[i]   = count;
				count      += numvalues_forcpu[i];
			}

			for(i=0;i<num_procs;i++){
				ISSM_MPI_Scatterv( row_indices_forcpu, sendcnts, displs, ISSM_MPI_INT, row_indices_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
				ISSM_MPI_Scatterv( values_forcpu, sendcnts, displs, TypeToMPIType<doubletype>(), values_fromcpu[i], numvalues_fromcpu[i], TypeToMPIType<doubletype>(), i, comm);
				ISSM_MPI_Scatterv( modes_forcpu, sendcnts, displs, ISSM_MPI_INT, modes_fromcpu[i], numvalues_fromcpu[i], ISSM_MPI_INT, i, comm);
			}
			/*}}}*/

			/*Plug values into global vector: {{{*/
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,comm);
			for(i=0;i<num_procs;i++){
				int  numvalues=numvalues_fromcpu[i];
				int* rows=row_indices_fromcpu[i];
				doubletype* values=values_fromcpu[i];
				int* mods=modes_fromcpu[i];

				for(j=0;j<numvalues;j++){
					if(mods[j]==ADD_VAL) *(vector+(rows[j]-lower_row))+=values[j];
					else *(vector+(rows[j]-lower_row))=values[j];
				}
			}
			/*}}}*/
			this->isassembled=true;

			/*Free resources:{{{*/
			xDelete<int>(RowRank);
			xDelete<int>(row_indices_forcpu);
			xDelete<int>(modes_forcpu);
			xDelete<doubletype>(values_forcpu);
			xDelete<int>(numvalues_forcpu);

			for(i=0;i<num_procs;i++){
				DataSet* bucketsn=bucketsforcpu[i];
				delete bucketsn;
			}
			xDelete<DataSet*>(bucketsforcpu);

			for(i=0;i<num_procs;i++){
				int* rows=row_indices_fromcpu[i];
				int* modes=modes_fromcpu[i];
				doubletype* values=values_fromcpu[i];

				xDelete<int>(rows);
				xDelete<int>(modes);
				xDelete<doubletype>(values);
			}
			xDelete<int*>(row_indices_fromcpu);
			xDelete<int*>(modes_fromcpu);
			xDelete<doubletype*>(values_fromcpu);
			xDelete<int>(numvalues_fromcpu);

			xDelete<int>(sendcnts);
			xDelete<int>(displs);

			/*Get rid of all buckets, as we have already added them to the matrix!: */
			delete this->buckets;
			this->buckets=new DataSet();
			/*}}}*/

		}
		/*}}}*/
		void SetValues(int ssize, int* list, doubletype* values, InsMode mode){/*{{{*/

			/*we need to store all the values we collect here in order to Assemble later.
			 * Indeed, the values we are collecting here most of the time will not belong
			 * to us, but to another part of the vector on another cpu: */
			_assert_(buckets);

			buckets->AddObject(new Bucket<doubletype>(ssize, list, values, mode));

		}
		/*}}}*/
		void SetValue(int dof, doubletype value, InsMode mode){/*{{{*/

			/*we need to store the value we collect here in order to Assemble later.
			 * Indeed, the value we are collecting here most of the time will not belong
			 * to us, but to another part of the vector on another cpu: */
			_assert_(buckets);

			buckets->AddObject(new Bucket<doubletype>(1,&dof,&value, mode));
		}
		/*}}}*/
		void GetValue(doubletype* pvalue,int dof){/*{{{*/
			_error_("Get value on a MpiVec vector not implemented yet!");
		}
		/*}}}*/
		void GetSize(int* pM){/*{{{*/

			*pM=this->M;

		}
		/*}}}*/
		void GetLocalSize(int* pM){/*{{{*/

			*pM=this->m;

		}
		/*}}}*/
		void GetLocalVector(doubletype** pvector,int** pindices){/*{{{*/

			/*First, check that vector size is not 0*/
			int vector_size;
			this->GetSize(&vector_size);
			if(vector_size==0){
				*pvector=NULL;
				*pindices=NULL;
				return;
			}

			/*Get Ownership range*/
			int lower_row,upper_row;
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,m,IssmComm::GetComm());
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
			_assert_(range==this->m);
			doubletype* values =xNew<doubletype>(range);
			xMemCpy<doubletype>(values,this->vector,this->m);

			*pvector  = values;
			*pindices = indices;
		} /*}}}*/
		IssmMpiVec<doubletype>* Duplicate(void){/*{{{*/

			return new IssmMpiVec<doubletype>(this->vector,this->M,this->m);

		}
		/*}}}*/
		void Set(doubletype value){/*{{{*/

			for(int i=0;i<this->m;i++)this->vector[i]=value;

		}
		/*}}}*/
		void AXPY(IssmAbsVec<doubletype>* Xin, doubletype a){/*{{{*/

			int i;

			/*Assume X is of the correct type, and downcast: */
			IssmMpiVec* X=NULL;

			X=(IssmMpiVec<doubletype>*)Xin;

			/*y=a*x+y where this->vector is y*/
			for(i=0;i<this->m;i++)this->vector[i]=a*X->vector[i]+this->vector[i];

		}
		/*}}}*/
		void AYPX(IssmAbsVec<doubletype>* Xin, doubletype a){/*{{{*/
			int i;

			/*Assume X is of the correct type, and downcast: */
			IssmMpiVec* X=NULL;

			X=(IssmMpiVec<doubletype>*)Xin;

			/*y=x+a*y where this->vector is y*/
			for(i=0;i<this->m;i++)this->vector[i]=X->vector[i]+a*this->vector[i];

		}
		/*}}}*/
		doubletype* ToMPISerial(void){/*{{{*/

			/*communicator info: */
			ISSM_MPI_Comm comm;
			int num_procs;

			/*ISSM_MPI_Allgatherv info: */
			int  lower_row,upper_row;
			int* recvcounts=NULL;
			int* displs=NULL;

			/*output: */
			doubletype* buffer=NULL;

			/*initialize comm info: */
			comm=IssmComm::GetComm();
			num_procs=IssmComm::GetSize();

			/*Allocate: */
// AD performance is sensitive to calls to ensurecontiguous.
// Providing "t" will cause ensurecontiguous to be called.
#ifdef _HAVE_AD_
			buffer=xNew<doubletype>(M,"t");
#else
			buffer=xNew<doubletype>(M);
#endif

			recvcounts=xNew<int>(num_procs);
			displs=xNew<int>(num_procs);

			/*recvcounts:*/
			ISSM_MPI_Allgather(&this->m,1,ISSM_MPI_INT,recvcounts,1,ISSM_MPI_INT,comm);

			/*get lower_row: */
			GetOwnershipBoundariesFromRange(&lower_row,&upper_row,this->m,comm);

			/*displs: */
			ISSM_MPI_Allgather(&lower_row,1,ISSM_MPI_INT,displs,1,ISSM_MPI_INT,comm);

			/*All gather:*/
			ISSM_MPI_Allgatherv(this->vector, this->m, TypeToMPIType<doubletype>(), buffer, recvcounts, displs, TypeToMPIType<doubletype>(),comm);
			/*Free resources: */
			xDelete<int>(recvcounts);
			xDelete<int>(displs);

			/*return: */
			return buffer;

		}
		/*}}}*/
		doubletype* ToMPISerial0(void){/*{{{*/

			/*FIXME: Should not broadcast to every cpu*/
			return this->ToMPISerial();

		}
		/*}}}*/
		void Shift(doubletype shift){/*{{{*/
			for(int i=0;i<this->m;i++)this->vector[i]+=shift;
		}
		/*}}}*/
		void Copy(IssmAbsVec<doubletype>* toin){/*{{{*/

			int i;

			/*Assume toin is of the correct type, and downcast: */
			IssmMpiVec* to=NULL;

			to=(IssmMpiVec<doubletype>*)toin;

			to->M=this->M;
			for(i=0;i<this->m;i++)to->vector[i]=this->vector[i];

		}
		/*}}}*/
		doubletype Norm(NormMode mode){/*{{{*/

			doubletype local_norm;
			doubletype norm;
			int i;

			switch(mode){
				case NORM_INF:
					local_norm=0.; for(i=0;i<this->m;i++)local_norm=max(local_norm,fabs(this->vector[i]));
					//local_norm=0; for(i=0;i<this->m;i++)local_norm=max(local_norm,this->vector[i]);
					ISSM_MPI_Reduce(&local_norm, &norm, 1, ISSM_MPI_DOUBLE, ISSM_MPI_MAX, 0, IssmComm::GetComm());
					ISSM_MPI_Bcast(&norm,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
					return norm;
					break;
				case NORM_TWO:
					local_norm=0.;
					for(i=0;i<this->m;i++)local_norm+=this->vector[i]*this->vector[i];
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
		void Scale(doubletype scale_factor){/*{{{*/

			int i;
			for(i=0;i<this->M;i++)this->vector[i]=scale_factor*this->vector[i];

		}
		/*}}}*/
		doubletype Dot(IssmAbsVec<doubletype>* inputin){/*{{{*/

			int i;
			doubletype local_dot=0;
			doubletype dot=0;

			/*Assume inputin is of the correct type, and downcast: */
			IssmMpiVec* input=NULL;

			input=(IssmMpiVec<doubletype>*)inputin;

			for(i=0;i<this->m;i++)local_dot+=this->vector[i]*input->vector[i];

			/*ISSM_MPI_SUM all the dots across the cluster: */
			ISSM_MPI_Reduce(&local_dot, &dot, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
			ISSM_MPI_Bcast(&dot,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());

			return dot;
		}
		/*}}}*/
		void PointwiseDivide(IssmAbsVec<doubletype>* xin,IssmAbsVec<doubletype>* yin){/*{{{*/

			int i;

			/*Assume xin and yin are of the correct type, and downcast: */
			IssmMpiVec* x=NULL;
			IssmMpiVec* y=NULL;

			x=(IssmMpiVec<doubletype>*)xin;
			y=(IssmMpiVec<doubletype>*)yin;

			/*pointwise w=x/y where this->vector is w: */
			for(i=0;i<this->m;i++)this->vector[i]=x->vector[i]/y->vector[i];
		}
		/*}}}*/
		void PointwiseMult(IssmAbsVec<doubletype>* xin,IssmAbsVec<doubletype>* yin){/*{{{*/

			int i;

			/*Assume xin and yin are of the correct type, and downcast: */
			IssmMpiVec* x=NULL;
			IssmMpiVec* y=NULL;

			x=(IssmMpiVec<doubletype>*)xin;
			y=(IssmMpiVec<doubletype>*)yin;

			/*pointwise w=x*y where this->vector is w: */
			for(i=0;i<this->m;i++)this->vector[i]=x->vector[i]*y->vector[i];
		}
		/*}}}*/
		void Pow(doubletype scale_factor){/*{{{*/
			for(int i=0;i<this->M;i++)this->vector[i]=pow(this->vector[i],scale_factor);
		}
		/*}}}*/
		void Sum(doubletype* pvalue){/*{{{*/

			doubletype local_sum=0;
			doubletype sum;

			for(int i=0;i<this->m;i++) local_sum+=this->vector[i];
			ISSM_MPI_Reduce(&local_sum, &sum, 1, ISSM_MPI_DOUBLE, ISSM_MPI_SUM, 0, IssmComm::GetComm());
			ISSM_MPI_Bcast(&sum,1,ISSM_MPI_DOUBLE,0,IssmComm::GetComm());
			*pvalue=sum;

		}
		/*}}}*/
		void BucketsBuildScatterBuffers(int** pnumvalues_forcpu,int** prow_indices_forcpu,doubletype** pvalues_forcpu,int** pmodes_forcpu,DataSet** bucketsforcpu,int num_procs){/*{{{*/

			/*intermediary: */
			int         i,j;
			int         count                   = 0;
			int         total_size              = 0;
			int        *temp_row_indices_forcpu = NULL;
			doubletype *temp_values_forcpu      = NULL;
			int        *temp_modes_forcpu       = NULL;

			/*output: */
			int        *numvalues_forcpu        = NULL;
			int        *row_indices_forcpu      = NULL;
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
			temp_row_indices_forcpu=row_indices_forcpu;
			temp_values_forcpu=values_forcpu;
			temp_modes_forcpu=modes_forcpu;

			/*Fill buffers: */
			for(i=0;i<num_procs;i++){
				DataSet    *buckets            = bucketsforcpu[i];
				for(j=0;j<buckets->Size();j++){
					Bucket<doubletype>* bucket =(Bucket<doubletype>*)buckets->GetObjectByOffset(j);
					bucket->Marshall(&temp_row_indices_forcpu,&temp_values_forcpu,&temp_modes_forcpu); //pass in the address of the buffers, so as to have the Marshall routine increment them.
				}
			}

			/*sanity check: */
			if (temp_row_indices_forcpu!=row_indices_forcpu+total_size)_error_("problem with marshalling of buckets");
			if (temp_values_forcpu!=values_forcpu+total_size)_error_("problem with marshalling of buckets");
			if (temp_modes_forcpu!=modes_forcpu+total_size)_error_("problem with marshalling of buckets");

			/*output buffers: */
			*pnumvalues_forcpu   = numvalues_forcpu;
			*prow_indices_forcpu = row_indices_forcpu;
			*pvalues_forcpu      = values_forcpu;
			*pmodes_forcpu       = modes_forcpu;
		}
		/*}}}*/
};
#endif //#ifndef _ISSM_MPI_VEC_H_
