      subroutine stot(ikval,qjadon,fltng_w,fltng_dwdt,Ntimp,Ntimm,
     1time,bi,dmi)
      implicit double precision (a-h,o-z)
      integer Ntimp,Ntimm
      parameter (Nafter = 1)
      parameter (nhank = 1024)
      double precision decay(2)
      double precision pset(7)
      double precision time(Ntimp),bi(Ntimm),dmi(Ntimm)
      double precision dekay1(nhank),dekay2(nhank),amp0(nhank),
     1amp1(nhank),amp2(nhank),amp3(nhank),amp4(nhank)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      common /blockm/ dekay1,dekay2,amp0,amp1,amp2,amp3,amp4
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  This subroutine returns the inverse Laplace transform to the
c  time-domain for the vertical displacement at time t for Hankel wavenumber
c  ikval.  (In general this routine needs to be called nhank times).
c  The main derivation uses the Faltung theorem of Laplace transforms.
c  (1-1-97)  NEW CASE OF 3-27-97 IS FOR t(Ntime + i) < t(Ntime) OR IN OTHER
c  WORDS, THE LOAD STILL IN PLACE AT t.  OPTION CALL to qwise.f
c  PERFORMS THIS. 
c 
c  A theory for the degenerate case was worked out but has been removed
c  as an option from this code.
c
c  Definition of tspan: nondimensional time span backwards form present
c                       when this routine is first called the dimensional
c                       equivalent might be say tspan = 12 ka, then 11 and
c                       then finally tspan = 0.
c  Additional note for r.s.l calculations: the routines qwise and pwise
c  are identical to the previous case for computations of present-day only
c  vertical deformation field.
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      decay(1) = dekay1(ikval)
      decay(2) = dekay2(ikval)
      xi0 = amp0(ikval)
      xi1 = amp1(ikval)
      xi2 = amp2(ikval)
      xi3 = amp3(ikval)
      xi4 = amp4(ikval)

c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c note that tspan must be updated in the calling routine "what0.f"
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      t = time(Ntimp)
      sumb_w = 0.0d0
      sumb_dwdt = 0.0d0
      ta = time(1)
      tb = time(2)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      do 97 i = 1,Ntimm
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c note that this "if" prevents adding load
c segments of "future" times when computing
c an r.s.l. history (10-06-98).
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      if(t.lt.ta) go to 97
      slope=dmi(i)
      ycept=bi(i)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      if( t . gt . ta . and . t . le . tb) go to 38
      call pwise(t,ta,tb,xi1,xi2,xi3,xi4,slope,ycept,decay,
     1bhaq_w,bhaq_dwdt)
      go to 39
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Note that qwise is employed only for the J-th Q hat term when t for
c evaluation still has to consider the load itself
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   38 call qwise(t,ta,qjadon,xi0,xi1,xi2,xi3,xi4,slope,ycept,decay,
     1bhaq_w,bhaq_dwdt)
   39 sumb_w = bhaq_w + sumb_w
      sumb_dwdt = bhaq_dwdt + sumb_dwdt
      ta = time(i + 1)
      tb = time(i + 2)
   97 continue
      fltng_w = sumb_w
      fltng_dwdt = sumb_dwdt
      return
      end
