      subroutine qwise(t,ta,qjadon,xi0,xi1,xi2,xi3,xi4,slope,ycept,
     1decay,bhaq_w,bhaq_dwdt)
      implicit double precision (a-h,o-z)
      double precision decay(2)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c This subroutine retrieves the convolution for the J-th linear piece-wise
c q hat function (the load shape or Bessel function part having been removed)
c with the free-decay solution. (see notes of 3-27-97 "convo.ice" Mathematica
c session).  The convolution is returned as "bhaq".
c
c  THIS ROUTINE REPLACES pwise.f ONLY FOR t <  time(Ntime) *
c  (such that the load is still in place at time t).       *
c
c Note irate = 1 case has to be applied to the linear term only (freed.f applies
c this correction to exponential terms)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xg1 = xi1/(decay(1)*decay(1))
      xg2 = xi2/(decay(2)*decay(2))
      xg3 = xi3/(decay(1)*decay(1))
      xg4 = xi4/(decay(2)*decay(2))
      gb1 = decay(1)*ycept
      gb2 = decay(2)*ycept
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c define xit0 term:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit0_w = (xi0 + qjadon) * ( ( slope * t ) + ycept )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c define xit1 term:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit1_w = xg1 * (
     1              gb1 + slope * ( ( t * decay(1) ) - 1.0d0 )
     2          - ( gb1 + slope * ( ( ta * decay(1) ) - 1.0d0 ))
     3                                   * dexp( decay(1) * (ta - t) )
     4                       )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c define xit2 term:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit2_w = xg2 * (
     1              gb2 + slope * ( ( t * decay(2) ) - 1.0d0 )
     2          - ( gb2 + slope * ( ( ta * decay(2) ) - 1.0d0 ) )
     3                                   * dexp( decay(2) * (ta - t) )
     4                       )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c And the rate equivalents:
c (sign switch due to freed.f already
c having corrected in x1t, x2t pass).
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit0_dwdt = (xi0 + qjadon) * slope 
      xit1_dwdt =-xg3 * (
     1              slope  
     2     + ( gb1 + slope * ( ( ta * decay(1) ) - 1.0d0 ))
     3                                   * dexp( decay(1) * (ta - t) )
     4                       )
      xit2_dwdt =-xg4 * (
     1              slope 
     2     + ( gb2 + slope * ( ( ta * decay(2) ) - 1.0d0 ))
     3                                   * dexp( decay(2) * (ta - t) )
     4                       )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c add terms for the J-th (and final) interval contribution.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      bhaq_w = xit0_w + xit1_w + xit2_w
      bhaq_dwdt = xit0_dwdt + xit1_dwdt + xit2_dwdt
      return
      end
