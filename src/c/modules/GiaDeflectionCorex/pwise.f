      subroutine pwise(t,ta,tb,xi1,xi2,xi3,xi4,slope,ycept,decay,
     1bhaq_w,bhaq_dwdt)
      implicit double precision (a-h,o-z)
      double precision decay(2)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c This subroutine retrieves the convolution for the ith linear piece-wise
c q hat function (the load shape or Bessel function part having
c been removed) with the free-decay solution. (see notes of
c 12-31-96 "Convolution in time").  The convolution is returned as "bhaq".
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      tbt = tb - t
      tat = ta - t
      gat1 = tat * decay(1)
      gat2 = tat * decay(2)
      gbt1 = tbt * decay(1)
      gbt2 = tbt * decay(2)
      ea1 = dexp(gat1)
      ea2 = dexp(gat2)
      eb1 = dexp(gbt1)
      eb2 = dexp(gbt2)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c define xit1 term:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit1 =(ycept/decay(1)) * (eb1 - ea1) -
     1(slope/(decay(1)*decay(1))) *
     2                            ( (1.0d0 - tb*decay(1))*eb1 
     3                            - (1.0d0 - ta*decay(1))*ea1 )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c define xit2 term:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      xit2 =(ycept/decay(2)) * (eb2 - ea2) -
     1(slope/(decay(2)*decay(2))) *
     2                            ( (1.0d0 - tb*decay(2))*eb2 
     3                            - (1.0d0 - ta*decay(2))*ea2 )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c add terms for the i-th interval contribution. 
c ABOVE IS THE NON-DEGENERATE CASE
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      bhaq_w = (xi1 * xit1) + (xi2 * xit2)    
      bhaq_dwdt = (xi3 * xit1) + (xi4 * xit2)    
      return
      end
