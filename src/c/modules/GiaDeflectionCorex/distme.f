      subroutine distme(Ntime,Ntimp,Ntimm,time,bi,dmi,zhload)
      implicit double precision (a-h,o-y)
      integer Ntime,Ntimp,Ntimm
      parameter (Nafter=1)
      double precision pset(7)
      double precision time(Ntimp),dmi(Ntimm),bi(Ntimm),dumbt(Ntimp)
      double precision hload(Ntime),qpat(Ntime),qt(Ntime)
      double precision zhload(Ntime),rhoi,distrad
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      common /blockp/ pset
      common /blockrad/ distrad 
      common /blocko/ rhoi
      data g /9.832186d0/, yearco /3.15576d7/, eradm/6.371d6/
      data dpi /3.1415926535897932d0/, dzero/0.0d0/
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c The units of time(Ntimp) are ka and the height of the load in meters.
c The slope, then for example, is in units of meters per ka.
c Note that "dumbt( )" is designed to perserve the initial "time( )" variable.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 776 k = 1, Ntimp
      dumbt(k) = time(k)
  776 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 39 itime = 1, Ntime
      hload(itime) = dble( zhload(itime) )
   39 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c now set up a piece-wise history: bi() = y-intercept 
c                                 dmi() = slope 
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 70 i = 2, Ntime
      dmi(i-1) = ( hload(i) - hload(i-1) )/( dumbt(i)  - dumbt(i-1) )
      bi(i-1) = hload(i-1) - ( dmi(i)*dumbt(i-1) )  
   70 continue
c      write(6,*) zhload(1,1), zhload(1,2) 
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c With pset(6) in mks units, lets convert the piecewise linear formulas
c for the time-dependent ice load heights to dimensionless values w.r.t. time.
c (tfact is in seconds)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      tfact = pset(2)/pset(4)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c get all times as dimensionless 
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 20 jt = 1, Nafter
      time(Ntime + jt) = ( dumbt(Ntime + jt) * yearco * 1.0d3 ) / tfact
   20 continue
      do 75 ind = 1, Ntimm 
      dmi(ind) =  dmi(ind) / (( yearco * 1.0d3 ) / tfact )
   75 continue
      do 77 j = 1, Ntime 
      time(j) = ( dumbt(j) * yearco * 1.0d3 ) / tfact
   77 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c create an incremental load in Pa and non-dimensionalized:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 80 iq = 1, Ntime
      qpat(iq) = hload(iq)*rhoi*g
      qt(iq) = qpat(iq) / pset(4)
   80 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c As the final step in this routine, create a dimensionless stress load from
c qp.  Here we'll use bi( ) and dmi( ) vectors with dimensionless time.  Then
c qp (and it's piece-wise decomposition) is ready for the direct dimensionless
c integrals for the inverse Laplace transform and inverse Hankel transform
c without further mutiplicative factors.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 85 i = 2, Ntime
      dmi(i-1) = ( qt(i) - qt(i-1) )/( time(i)  - time(i-1) )
      bi(i-1) = qt(i-1) - ( dmi(i-1)*time(i-1) )  
   85 continue
  999 return
      end
