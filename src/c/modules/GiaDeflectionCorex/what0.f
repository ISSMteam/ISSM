      subroutine what0(iedge,Ntimp,Ntimm,time,bi,dmi)
      implicit double precision (a-h,o-z)
      integer Ntimp,Ntimm
      parameter (nhank = 1024)
      parameter (N = nhank/2)
      double precision dekay1(nhank),dekay2(nhank),amp0(nhank),
     1amp1(nhank),amp2(nhank),amp3(nhank),amp4(nhank),
     1zksam(nhank),zksamp(nhank)
      double precision decay(2),pset(7),amps(5),
     1decta(2),dyri1(nhank),dyri2(nhank),sna(nhank)
      double precision cinner_w(nhank),cinner_dwdt(nhank)
      double precision bcin_w(nhank),bcin_dwdt(nhank)
      double precision time(Ntimp),bi(Ntimm),dmi(Ntimm)
      integer maxk
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      common /blockp/ pset
      common /blockz/ zkp
      common /blockm/ dekay1,dekay2,amp0,amp1,amp2,amp3,amp4
      data yearco /3.15576d7/, pi /3.1415926535897932384d0/
      data g /9.832186d0/, four /4.d0/, two /2.0d0/,
     1 one /1.0d0/, zero/0.0d0/ , maxk/64/
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      twopi = two * pi
      r2 = pset(6)
      u2 = pset(4)
      r1 = pset(5)
      u1 = pset(3)
      h  = pset(1)
      urat = u1/u2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  alphap is dimensionless disk radius
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      alphap = pset(7)/pset(1)
      twoap = two * alphap
      rghm = ( r1 * g * h * alphap ) / (two * u2)
      taumx = pset(2)/pset(4)
      tmxyr = taumx / yearco
c
      dfac = dfloat(nhank)/dfloat(maxk)
      endk = dfloat(nhank)/dfac
      dk = endk/dfloat(nhank)
c
      ak = zero
      do 7000 ik = 1,nhank
      ak = ak + dk
      zkp = ak
      pikn = (6.371d6 * zkp) / h
      zkd = pikn / 6.371d6
c
      zkp2 = 2.0d0 * zkp
      zkp4 = 4.0d0 * zkp
      e1 = dexp(zkp)
      e2 = dexp(zkp2)
      e4 = dexp(zkp4)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      call freed(r2,u2,r1,u1,h,zkd,e1,e2,e4,b0,b1,a2,a1,a0,decay,amps)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      decta(1) = decay(1)/tmxyr
      decta(2) = decay(2)/tmxyr
      dyri1(ik) = decta(1)
      dyri2(ik) = decta(2)
      sna(ik) = pikn
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Form vectors for full construction in pwise.f and stot.f
c Note that freed will produce decay spectra defined as positive, 
c ie. negative decay must reinsert a minus sign.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      dekay1(ik) = decay(1)
      dekay2(ik) = decay(2)
      amp0(ik) = amps(3)  
      amp1(ik) = amps(1)  
      amp2(ik) = amps(2)  
      amp3(ik) = amps(4)  
      amp4(ik) = amps(5)  
      zksam(ik) = zkd
      zksamp(ik) = zkp
 7000 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c The following looped call sets up the free solution convolved with the
c load function q hat.  Note that the returned vector set "cinner" is the
c inner-most part of the arguement of the inverse Hankel trans. integral.
c It is time-dependent and the loop is for the k-dependancy. The time for
c calculation is given in the vector "time(Ntimp)" in the routine stot.f that is
c called below. Note that the sign on cinner(ik) below is for a load directed
c downward.   ** For iedge = 1 assume sq. edge load and for iedge = 2 assume an
c elliptical cross section.  Note loops 8500,8000 and 9500,9000 for the two
c cases, respectively.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   49 go to (8499,9499), iedge
 8499 do 8000 ik = 1, nhank
      xakap = zksamp(ik)*alphap
      diku = xakap * urat
      pref = diku / ( diku + rghm )
      qjadon = one / ( four * zksamp(ik) * urat )
      call stot(ik,qjadon,fltng_w,fltng_dwdt,Ntimp,Ntimm,time,bi,dmi)
      cinner_w(ik) = - fltng_w * pref * twoap
      cinner_dwdt(ik) = - fltng_dwdt * pref * twoap
      bcin_w(ik) = cinner_w(ik) * dbesj1(xakap)
      bcin_dwdt(ik) = cinner_dwdt(ik) * dbesj1(xakap)
 8000 continue
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c "ojrule.f" computes the inverse Hankel trasform with a simple
c Simpson's rule.  The routine "ojrule" is buliding a set of solutions stored
c in common "blocks" in r or "asrpos(nrv) ", and computed rate or displacement
c for each of N3G disks in "aswokm(nrv,N3G)" . 
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      call ojrule(dk,bcin_w,bcin_dwdt)
      go to 999
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 9499 do 9000 ik = 1, nhank
      xakap = zksamp(ik)*alphap
      oxakap = one/xakap
      diku = xakap * urat
      pref = diku / ( diku + rghm )
      qjadon = one / ( four * zksamp(ik) * urat )
      call stot(ik,qjadon,fltng_w,fltng_dwdt,Ntimp,Ntimm,time,bi,dmi)
      cinner_w(ik) = - fltng_w * pref * twoap
      cinner_dwdt(ik) = - fltng_dwdt * pref * twoap
      bcin_w(ik) = cinner_w(ik) * oxakap * ( dsin(xakap) * oxakap
     1 - dcos(xakap) )
      bcin_dwdt(ik) = cinner_dwdt(ik) * oxakap * ( dsin(xakap) * oxakap
     1 - dcos(xakap) )
 9000 continue
      call ojrule(dk,bcin_w,bcin_dwdt)
  999 return

      end
