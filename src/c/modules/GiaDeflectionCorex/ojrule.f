      subroutine ojrule(dk,bcin_w,bcin_dwdt)
      implicit double precision(a-h,o-z)
      parameter (nhank = 1024)
      double precision yvalue_w(nhank),yvalue_dwdt(nhank)
      double precision bcin_w(nhank),bcin_dwdt(nhank)
      double precision wok_w,wok_dwdt,rpos
      double precision swok_w,swok_dwdt
      double precision pset(7)
      double precision aswokm_w,aswokm_dwdt,distrad
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      common /blockrad/ distrad
      common /blockp/ pset
      common /blocks/ aswokm_w,aswokm_dwdt
      data zero /0.0d0/, one /1.0d0/, two /2.0d0/, three /3.0d0/,
     1rescal/ 1.0d0/
      data yearco /3.15576d7/
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      bath = dk / three
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c rpos should be normalized wrt lithosphere thickness 
c give r is normalized dist_rad :: r == dist_rad / h
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      r = distrad / (pset(1) / 1.0d3)
      rpos = r 
      ak = zero
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c form the yvalue's for the Simpson's rule formulas
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 425 ik = 1, nhank
      ak = ak + dk
      rak = ak * r
      rarg = dbesj0( rak )
      yvalue_w(ik) = bcin_w(ik) * rarg
      yvalue_dwdt(ik) = bcin_dwdt(ik) * rarg
  425 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c correct to end point val. in Simp. Rule
c      yvalue(nhank) = bcin(nhank) * rarg / two
c find the area under the curve using the Simpson's rule formulas
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      sumde_w = zero
      sumde_dwdt = zero
      do 300 int = 1, nhank
      intp1 = int + 1
      ide = 2 + ( (-1)**intp1 + 1 )
      fide = dfloat(ide)
      sumde_w = ( fide * yvalue_w(int) ) + sumde_w
      sumde_dwdt = ( fide * yvalue_dwdt(int) ) + sumde_dwdt
  300 continue
      wok_w = bath * sumde_w
      wok_dwdt = bath * sumde_dwdt
      
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      hscale = sngl(pset(1))
      hsckm = hscale / 1.0e3
      swok_w = hscale * sngl(wok_w)
      aswokm_w = swok_w
      swok_dwdt = (hscale * yearco * 1.0e3 * sngl(wok_dwdt))
     1                  * ( sngl(pset(4))/ sngl(pset(2)) )
      aswokm_dwdt = swok_dwdt
      return

      end
