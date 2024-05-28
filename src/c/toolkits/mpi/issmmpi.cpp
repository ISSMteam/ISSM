/* \file issmmpi.cpp
 * \brief: implementation of all the mpi wrappers that ISSM requires. The goal is to control
 * which MPI layer we are using at compile time: the standard mpi, the autodiff mpi or no mpi at all.
 */

#ifdef HAVE_CONFIG_H
	#include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include <cassert>
#include <cstring> // for memcpy

#include "./issmmpi.h"

#if defined(_HAVE_MPI_) && defined(_HAVE_CODIPACK_)
MpiTypes* mpiTypes;
#endif

#include "../../shared/Numerics/types.h"

#ifndef _HAVE_MPI_
ISSM_MPI_Status ourIssmMPIStatusIgnore=0;
size_t sizeHelper(ISSM_MPI_Datatype type) { /*{{{*/

  switch(type) {
  case ISSM_MPI_CHAR:
    return sizeof(char);
    break;
  case ISSM_MPI_DOUBLE:
    return sizeof(double);
    break;
  case ISSM_MPI_INT:
    return sizeof(int);
    break;
  default:
    assert(0);
    break;
  }
  return 0;
}/*}}}*/
#endif

int ISSM_MPI_Allgather(void *sendbuf, int sendcount, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcount, ISSM_MPI_Datatype recvtype, ISSM_MPI_Comm comm) {  /*{{{*/
  int rc=0;
  assert(sendcount==recvcount || sendtype==recvtype); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Allgather(sendbuf,
		    sendcount,
		    sendtype,
		    recvbuf,
		    recvcount,
		    recvtype,
		    comm);
# else
  rc=MPI_Allgather(sendbuf,
		   sendcount,
		   sendtype,
		   recvbuf,
		   recvcount,
		   recvtype,
		   comm);
# endif
#else
# ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<sendcount;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,sendbuf,sizeHelper(sendtype)*sendcount);
#endif
  return rc;
} /*}}}*/
int ISSM_MPI_Allgatherv(void *sendbuf, int sendcount, ISSM_MPI_Datatype sendtype, void *recvbuf, int *recvcounts, int *displs, ISSM_MPI_Datatype recvtype, ISSM_MPI_Comm comm) {  /*{{{*/
  int rc=0;
  assert(sendtype==recvtype); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Allgatherv(sendbuf,
		     sendcount,
		     sendtype,
		     recvbuf,
		     recvcounts,
		     displs,
		     recvtype,
		     comm);
# else
  rc=MPI_Allgatherv(sendbuf,
		    sendcount,
		    sendtype,
		    recvbuf,
		    recvcounts,
		    displs,
		    recvtype,
		    comm);
# endif
#else
  assert(sendcount==recvcounts[0]); // we handle only identical representations
# ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)(recvbuf)+displs[0];
    for(int i=0;i<sendcount;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy((char*)recvbuf+(sizeHelper(recvtype)*displs[0]),sendbuf,sizeHelper(sendtype)*sendcount);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Allreduce(void *sendbuf, void *recvbuf, int count, ISSM_MPI_Datatype datatype, ISSM_MPI_Op op, ISSM_MPI_Comm comm){/*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Allreduce(sendbuf,
		    recvbuf,
		    count,
		    datatype,
		    op,
		    comm);
# else
  rc=MPI_Allreduce(sendbuf,
		   recvbuf,
		   count,
		   datatype,
		   op,
		   comm);
# endif
#else
#ifdef _HAVE_AD_
  if (datatype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<count;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,sendbuf,sizeHelper(datatype)*count);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Barrier(ISSM_MPI_Comm comm){  /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Barrier(comm);
# else
  rc=MPI_Barrier(comm);
# endif
#else
// do nothing
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Bcast(void *buffer, int count, ISSM_MPI_Datatype datatype, int root, ISSM_MPI_Comm comm){  /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Bcast(buffer,
		count,
		datatype,
		root,
		comm);
# else
  rc=MPI_Bcast(buffer,
	       count,
	       datatype,
	       root,
	       comm);
# endif
#else
// nothing to be done here
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Comm_free(ISSM_MPI_Comm *comm){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  assert(0); // to be implemented
# else
  rc=MPI_Comm_free(comm);
# endif
#else
// do nothing
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Comm_rank(ISSM_MPI_Comm comm, int *rank){  /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
  rc=MPI_Comm_rank(comm,
		   rank);
#else
  *rank=0;
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Comm_size( ISSM_MPI_Comm comm, int *size){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
  rc=MPI_Comm_size(comm,
		   size);
#else
  *size=1;
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Finalize(void){  /*{{{*/

	int rc=0;
	#ifdef _HAVE_MPI_
		#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
			#if defined(_HAVE_ADJOINTMPI_)
				rc=AMPI_Finalize();
			#elif defined(_HAVE_MEDIPACK_)
				/*Old implementation*/
				//TOOL::finalize();
				/*New implementation*/
				delete mpiTypes;
				rc=AMPI_Finalize();
			#else
				rc=AMPI_Finalize_NT();
			#endif
		#else
		  rc=MPI_Finalize();
		#endif
	#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Gather(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm){  /*{{{*/

  int rc=0;
  assert(sendtype==recvtype && sendcnt==recvcnt); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Gather(sendbuf,
		 sendcnt,
		 sendtype,
		 recvbuf,
		 recvcnt,
		 recvtype,
		 root,
		 comm);
# else
  rc=MPI_Gather(sendbuf,
		sendcnt,
		sendtype,
		recvbuf,
		recvcnt,
		recvtype,
		root,
		comm);
# endif
#else
# ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<sendcnt;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,sendbuf,sizeHelper(sendtype)*sendcnt);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Gatherv(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int *recvcnts, int *displs, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm){/*{{{*/

  int rc=0;
  assert(sendtype==recvtype); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Gatherv(sendbuf,
		  sendcnt,
		  sendtype,
		  recvbuf,
		  recvcnts,
		  displs,
		  recvtype,
		  root,
		  comm);
# else
  rc=MPI_Gatherv(sendbuf,
		 sendcnt,
		 sendtype,
		 recvbuf,
		 recvcnts,
		 displs,
		 recvtype,
		 root,
		 comm);
# endif
#else
  assert(sendcnt==recvcnts[0]); // we handle only identical representations
#ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)(recvbuf)+displs[0];
    for(int i=0;i<sendcnt;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy((char*)recvbuf+(sizeHelper(recvtype)*displs[0]),sendbuf,sizeHelper(sendtype)*sendcnt);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Init(int *argc, char ***argv){  /*{{{*/

	int rc=0;
	#ifdef _HAVE_MPI_
		#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
			#if defined(_HAVE_ADJOINTMPI_)
				rc=AMPI_Init(argc,argv);
			#elif defined(_HAVE_MEDIPACK_)
				rc=AMPI_Init(argc,argv);
				/*Old implementation of Medipack*/
				//TOOL::init();
				/*New*/
				//MpiTypes* mpiTypes;
				mpiTypes = new MpiTypes();
			#else
				rc=AMPI_Init_NT(argc,argv);
			#endif
		#else
			rc=MPI_Init(argc,argv);
		#endif
	#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Recv(void *buf, int count, ISSM_MPI_Datatype datatype, int source, int tag, ISSM_MPI_Comm comm, ISSM_MPI_Status *status){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Recv(buf,
	       count,
	       datatype,
	       source,
	       tag,
			 #if !defined(_HAVE_ADJOINTMPI_) && !defined(_HAVE_MEDIPACK_)
	       AMPI_FROM_SEND, // as long as there are no other variants
			 #endif
	       comm,
	       status);
# else
  rc=MPI_Recv(buf,
	      count,
	      datatype,
	      source,
	      tag,
	      comm,
	      status);
# endif
#else
// nothing to be done here
// as long as nobody tries to do anything with 'status'
// we won't do anything to it here either
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Reduce(void *sendbuf, void *recvbuf, int count, ISSM_MPI_Datatype datatype, ISSM_MPI_Op op, int root, ISSM_MPI_Comm comm){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Reduce(sendbuf,
		 recvbuf,
		 count,
		 datatype,
		 op,
		 root,
		 comm);
# else
  rc=MPI_Reduce(sendbuf,
		recvbuf,
		count,
		datatype,
		op,
		root,
		comm);
# endif
#else
# ifdef _HAVE_AD_
  if (datatype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<count;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,sendbuf,sizeHelper(datatype)*count);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Scatter(void *sendbuf, int sendcnt, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm){ /*{{{*/

  int rc=0;
  assert(sendtype==recvtype && sendcnt==recvcnt); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Scatter(sendbuf,
		  sendcnt,
		  sendtype,
		  recvbuf,
		  recvcnt,
		  recvtype,
		  root,
		  comm);
# else
  rc=MPI_Scatter(sendbuf,
		 sendcnt,
		 sendtype,
		 recvbuf,
		 recvcnt,
		 recvtype,
		 root,
		 comm);
# endif
#else
# ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)sendbuf;
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<recvcnt;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,sendbuf,sizeHelper(sendtype)*recvcnt);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Scatterv(void *sendbuf, int *sendcnts, int *displs, ISSM_MPI_Datatype sendtype, void *recvbuf, int recvcnt, ISSM_MPI_Datatype recvtype, int root, ISSM_MPI_Comm comm){ /*{{{*/

  int rc=0;
  assert(sendtype==recvtype); // we handle only identical representations
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Scatterv(sendbuf,
		   sendcnts,
		   displs,
		   sendtype,
		   recvbuf,
		   recvcnt,
		   recvtype,
		   root,
		   comm);
# else
  rc=MPI_Scatterv(sendbuf,
		  sendcnts,
		  displs,
		  sendtype,
		  recvbuf,
		  recvcnt,
		  recvtype,
		  root,
		  comm);
# endif
#else
  assert(sendcnts[0]==recvcnt); // we handle only identical representations
# ifdef _HAVE_AD_
  if (sendtype==ISSM_MPI_DOUBLE) {
    IssmDouble* activeSendBuf=(IssmDouble*)(sendbuf)+displs[0];
    IssmDouble* activeRecvBuf=(IssmDouble*)recvbuf;
    for(int i=0;i<recvcnt;++i) activeRecvBuf[i]=activeSendBuf[i];
  }
  else
# endif
    memcpy(recvbuf,(char*)sendbuf+(sizeHelper(sendtype)*displs[0]),sizeHelper(sendtype)*recvcnt);
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Send(void *buf, int count, ISSM_MPI_Datatype datatype, int dest, int tag, ISSM_MPI_Comm comm){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Send(buf,
	       count,
	       datatype,
	       dest,
	       tag,
			 #if !defined(_HAVE_ADJOINTMPI_) && !defined(_HAVE_MEDIPACK_)
	       AMPI_TO_RECV, // as long as there are no other variants
			 #endif
	       comm);
# else
  rc=MPI_Send(buf,
	      count,
	      datatype,
	      dest,
	      tag,
	      comm);
# endif
#else
// nothing to be done here
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Isend(void *buf, int count, ISSM_MPI_Datatype datatype, int dest, int tag, ISSM_MPI_Comm comm, ISSM_MPI_Request* req){ /*{{{*/

  int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
  rc=AMPI_Isend(buf,
	       count,
	       datatype,
	       dest,
	       tag,
			 #if !defined(_HAVE_ADJOINTMPI_) && !defined(_HAVE_MEDIPACK_)
	       AMPI_TO_RECV, // as long as there are no other variants
			 #endif
	       comm,
		   req);
# else
  rc=MPI_Isend(buf,
	      count,
	      datatype,
	      dest,
	      tag,
	      comm,
	      req);
# endif
#else
// nothing to be done here
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Wait(ISSM_MPI_Request *req, ISSM_MPI_Status *status){/*{{{*/
	int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
	rc=AMPI_Wait(req, status);
# else
	rc=MPI_Wait(req, status);
# endif
#else
// nothing to be done here
#endif
	return rc;
}/*}}}*/
double ISSM_MPI_Wtime(void){/*{{{*/

#ifdef _HAVE_MPI_
	return MPI_Wtime();
#else
	assert(0); // to be implemented
	return 0.0;
#endif
}/*}}}*/
void ISSM_MPI_ContiguousInAdolc(size_t aSize) { /*{{{*/

#ifdef _HAVE_ADOLC_
  ensureContiguousLocations(aSize);
#else
  fprintf(stderr, "*** Codipack ISSM_MPI_ContiguousInAdolc()\n");
#endif
}/*}}}*/
int ISSM_MPI_Comm_split(ISSM_MPI_Comm comm, int color, int key, ISSM_MPI_Comm *newcomm){ /*{{{*/

int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
rc=MPI_Comm_split(comm, color, key, newcomm);
#else
rc=MPI_Comm_split(comm, color, key, newcomm);
#endif
#else
// nothing to be done here
#endif
  return rc;
}/*}}}*/
int ISSM_MPI_Intercomm_create(ISSM_MPI_Comm comm,int local_leader,ISSM_MPI_Comm peer_comm, int remote_leader, int tag,ISSM_MPI_Comm *newintercomm){ /*{{{*/

	int rc=0;
#ifdef _HAVE_MPI_
#if defined(_HAVE_AMPI_) &&  !defined(_WRAPPERS_)
	rc=MPI_Intercomm_create(comm,local_leader,peer_comm,remote_leader,tag,newintercomm);
#else
	rc=MPI_Intercomm_create(comm,local_leader,peer_comm,remote_leader,tag,newintercomm);
#endif
#else
	// nothing to be done here
#endif
	return rc;
}/*}}}*/
