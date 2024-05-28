#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: CPL_ISSM
C     !INTERFACE:
      SUBROUTINE CPL_ISSM( myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE cpl_issm
C     | o Couple MITgcm ocean model with ISSM ice sheet model
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "DYNVARS.h"
#include "GRID.h"
#include "FFIELDS.h"
#include "SHELFICE_OPTIONS.h"
#include "SHELFICE.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
# include "EXF_FIELDS.h"
#endif

      LOGICAL  DIFFERENT_MULTIPLE
      EXTERNAL DIFFERENT_MULTIPLE

C     !LOCAL VARIABLES:
C     mytime - time counter for this thread (seconds)
C     myiter - iteration counter for this thread
C     mythid - thread number for this instance of the routine.
      _RL     mytime
      INTEGER myiter, mythid 
CEOP

#ifdef ALLOW_CPL_ISSM
#include "EESUPPORT.h"
      COMMON /CPL_MPI_ID/ mpiMyWid, toissmcomm
      INTEGER mpiMyWid, toissmcomm, mpiRC
      INTEGER mpistatus(MPI_STATUS_SIZE)
      INTEGER i, j, bi, bj, buffsize
      COMMON /CPL_ISSM_TIME/ CouplingTime
      _R8 CouplingTime, IceModelTime
      _R8 xfer_array(Nx,Ny)
      _R8 local(1:sNx,1:sNy,nSx,nSy)
      CHARACTER*(MAX_LEN_MBUF) suff

C Initialization steps I1, I2, and I3:
      IF( myTime .EQ. startTime ) THEN

C   I1. ISSM sends CouplingTime, the interval at which we couple
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Recv(CouplingTime,1,MPI_DOUBLE,0,10001000,
     &           toissmcomm,mpistatus,mpiRC)
            _END_MASTER( myThid )
         ENDIF
         _BEGIN_MASTER( myThid )
         CALL MPI_BCAST(CouplingTime,1,MPI_DOUBLE,0,
     &        MPI_COMM_MODEL,mpiRC)
         _END_MASTER( myThid )
#ifdef ALLOW_DEBUG
         IF ( debugLevel .GE. debLevD ) THEN
            print*, 'Ocean received CouplingTime: ', CouplingTime
         ENDIF
#endif

C   I2. MITgcm sends grid size (NX and NY)
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Send(Nx,1,MPI_INT,0,10001003,
     &           toissmcomm,mpistatus)
            call MPI_Send(Ny,1,MPI_INT,0,10001004,
     &           toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF

C   I3. MITgcm sends grid coordinates of center of cells
C       (longitude -180 <= XC < 180 and latitude YC)
C     Send longitude East of center of cell
         DO bj=1,nSy
            DO bi=1,nSx
               DO j=1,sNy
                  DO i=1,sNx
                     local(i,j,bi,bj) = xC(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL BAR2( myThid ) 
         CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            buffsize = Nx*Ny
            CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001005,toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid )
C     Send latitude North of center of cell
         DO bj=1,nSy
            DO bi=1,nSx
               DO j=1,sNy
                  DO i=1,sNx
                     local(i,j,bi,bj) = yC(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL BAR2( myThid ) 
         CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            buffsize = Nx*Ny
            CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001006,toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid )

      ENDIF
C End initialization steps I1, I2, and I3.

C Recurring steps C1 and C2:
      IF( MOD(myTime,CouplingTime) .LT. deltaT/2. ) THEN

C   C1. ISSM sends ice model time IceTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Recv(IceModelTime,1,MPI_DOUBLE,0,10001001,
     &           toissmcomm,mpistatus,mpiRC)
#ifdef ALLOW_DEBUG
            IF ( debugLevel .GE. debLevD ) THEN
               print*, 'Ocean received IceModelTime: ', IceModelTime
            ENDIF
#endif
            _END_MASTER( myThid )
         ENDIF

C   C2. MITgcm sends ocean model time OceanTimeTag
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            call MPI_Send(myTime,1,MPI_DOUBLE,0,10001002,
     &           toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF

      ENDIF
C End recurring steps C1 and C2.

C Recurring step C3 except during Initialization:
C  C3. MITgcm sends
C      (N-1)*CouplingTime <= OceanModelTime < N*CouplingTime
C      time-mean melt rate to ISSM
      IF( myTime .NE. startTime .AND.
     &     MOD(myTime,CouplingTime) .LT. deltaT/2. ) THEN
         DO bj=1,nSy
            DO bi=1,nSx
               DO j=1,sNy
                  DO i=1,sNx
                     local(i,j,bi,bj)=shelficeFreshWaterFlux(i,j,bi,bj)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         CALL BAR2( myThid ) 
         CALL GATHER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )
            buffsize = Nx*Ny
            CALL MPI_SEND(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001007,toissmcomm,mpistatus)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid )
#ifdef ALLOW_DEBUG
         IF ( debugLevel .GE. debLevD ) THEN
            print*,'Done Sending shelficeFreshWaterFlux array.'
         ENDIF
#endif
      ENDIF
C End recurring step C3.

C Recurring step C4 except during Termination:
C  C4. ISSM sends IceModelTime=(N-1)*CouplingTime base to MITgcm
      IF( myTime .NE. endtime .AND.
     &     MOD(myTime,CouplingTime) .LT. deltaT/2. ) THEN
         WRITE(suff,'(I10.10)') myIter
         CALL WRITE_FLD_XY_RS( 'R_shelfIce1_',suff,R_shelfIce,-1,myThid)
         IF( myProcId .EQ. 0 ) THEN
            _BEGIN_MASTER( myThid )         
            call MPI_Recv(xfer_array,buffsize,MPI_DOUBLE_PRECISION,
     &           0,10001008,toissmcomm,mpistatus,mpiRC)
            _END_MASTER( myThid )
         ENDIF
         CALL BAR2( myThid ) 
         CALL SCATTER_2D_R8( xfer_array, local, Nx, Ny,
     &        .FALSE., .FALSE., myThid )
         DO bj = myByLo(myThid), myByHi(myThid)
            DO bi = myBxLo(myThid), myBxHi(myThid)
               DO j=1,sNy
                  DO i=1,sNx
                     IF( local(i,j,bi,bj).LT.9998 ) THEN
                        R_shelfIce(i,j,bi,bj) = local(i,j,bi,bj)
                     ELSE
                        R_shelfIce(i,j,bi,bj) = 0. _d 0
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
C- fill in the overlap (+ BARRIER):
         _EXCH_XY_RS( R_shelfIce, myThid )
         CALL WRITE_FLD_XY_RS( 'R_shelfIce2_',suff,R_shelfIce,-1,myThid)
      ENDIF
C End recurring step C4.

#endif /* ALLOW_CPL_ISSM */

      RETURN
      END
