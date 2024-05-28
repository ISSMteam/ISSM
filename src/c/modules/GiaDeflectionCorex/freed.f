      subroutine freed(r2,u2,r1,u1,h,zk,e1,e2,e4,b0,b1,a2,a1,a0,decay
     1,amps)        
      implicit double precision (a-h,o-z)
      double precision decay(2),amps(5)
      double precision ac0,ac1,ac2,ac3,ac4,ac5,ac6,ac7,ac8,ac9,ac10,
     1ac11
      common /blockz/ zkp
      data zero /0.0d0/, g /9.832186d0/
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Given the inputs to this subroutine(r2 through zk in the call
c statement above), the outputs are coefficients of "s" that
c are crucial to the Laplace transform inversion. From b0 and b1
c we can compute the decay poles (or eigenvalues).    
c 
c  This is NOT true in our case, though. - SA
c  NOTE IN THE CODE THAT A CALL TO THIS SUBROUTINE NEED NOT
c  BE MADE AT EACH TIME STEP --- BUT WILL HAVE TO BE CALLED
c  IN THE NUMERICAL INTEGRATION FOR COMPUTING THE INVERSE HANKEL
c  TRANSFORM  ****
c
c Each term should be returned as dimensionless 
c h => length   u2 => stress     taumx2 => time
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      zkp2 = zkp*zkp
      ur = u1/u2
      ghu2 = (g*h) / u2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac0 dimensional units are stress times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac0 = 4.0d0*ur*zkp2*( 1.0d0 + e4 +
     1    2.0d0*e2*(1.0d0 + 2.0d0*zkp2) )
c    DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac1 dimensional units are stress times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac1 = 2.0d0*r1*ghu2*zkp*(1.0d0 - e4 + 4.0d0*zkp*e2)
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac2 dimensional units are stress^2 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac2 = 8.0d0*ur*ur*(-1.0d0 + e1)*
     1                     (1.0d0 + e1)*(1.0d0 + e2)*zkp2
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac3  dimensional units are stress^2 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac3 =
     1   2.0d0*zkp*ghu2*ur*((r1 + r2)*(1.0d0 + e4) + 
     2       2.0d0*(r2 - r1)*e2*( 1.0d0 + 2.0d0*zkp2 ))
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac4 dimensional units are stress^2 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac4 = ghu2*ghu2*r1*(r2 - r1)*
     1      (1.0d0 - e4 + 4.0d0*zkp*e2)
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac5 dimensional units are stress^3 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac5 = 
     1    4.0d0*zkp2*ur*ur*ur*(1.0d0 - e2 - 2.0d0*e1*zkp)*
     2                (1.0d0 - e2 + 2.0d0*e1*zkp)
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac6 dimensional units are stress^3 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac6 =
     1   2.0d0*zkp*ur*ur*(1.0d0 - e4
     2               - 4.0d0*e2*zkp)*ghu2*r2
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac7 dimensional units are stress^3 times l^-2
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac7 =
     1 ur * ( ( (1.0d0 - e1)*(1.0d0 + e1) )**2)*r1*(r2 - r1)
     2 * ( ghu2*ghu2 )
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac8 dimensional units are stress^0 times l^-1
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac8 = -2.0d0*zkp*(1.0d0 + e2*(1.0d0 + 2.0d0*zkp*(1.0d0 + zkp))) 
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac9 dimensional units are stress^1 times l^-1
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac9 = 
     1  ( 4.0d0*zkp*u1 -
     2 g*h*(r2 - r1)*(1.0d0 + e2*(1.0d0 + 2.0d0*zkp*(1.0d0 + zkp)))
     3    ) / u2
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac10 dimensional units are stress^2 times l^-1
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac10 =
     1   -2.0d0*zkp*ur*ur*( 1.0d0 - e2
     2  - 2.0d0*zkp*e2*(1.0d0 + zkp) )
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ac11 dimensional units are stress^2 times l^-1
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      ac11 =
     1 ghu2*ur*(r2 - r1)*(1.0d0 - e2*(1.0d0 + 2.0d0*zkp))
c     DIMESIONLESS
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c Other functions may be found in file "apcw0.record"
c (Nov. 9 1996)
c The following is a Mathematica version of the isolation of the
c  coefficeints of the L transform variable s in the denominator.  
c  Here is where the set-up is performed to obtain the "free decay"
c  times (with the Hankel transform variable "zk" embedded.  Note that
c  some greater efficency could be achieved by further simplifying the
c  combinations of "acn" functions which are now a series of function
c  subroutines in the fortran code.  The corresponding Mathematica
c  session is "twolayer.Linversion" dated Nov. 23, 1996.
c
c In[59]:=
c Together[%]
c Out[59]=
c    ac2 + ac3 - ac4 + 2 ac5 - 2 ac6 + 2 ac7
c ---------------------------------------------
c ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c In[61]:=
c Simplify[Coefficient[els,s^2]]
c Out[61]=
c 1
c In[65]:=
c eslnos =
c ac5/(ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7) - 
c    ac6/(ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7) + 
c 
c   ac7/(ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7)
c    
c Out[65]=
c                      ac5
c --------------------------------------------- - 
c ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c 
c                        ac6
c  --------------------------------------------- + 
c   ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c 
c                        ac7
c  ---------------------------------------------
c  ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c In[66]:=
c Simplify[%]
c Out[66]=
c                ac5 - ac6 + ac7
c ---------------------------------------------
c ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c a common denominator factor is: bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      bc =  (   ac0 - ac1 +
     1                     ac2 + ac3
     2                                 - ac4 + ac5 -
     3           ac6 + ac7  )
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c b1: Denominator coefficent of s:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      b1 =
     1  (  ac2 + ac3
     2               - ac4 + ( 2.0d0 * ac5 )
     3                                       - ( 2.0d0 * ac6 )
     4                                       + ( 2.0d0 * ac7 )   ) / bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c b0: Denominator coefficent of s^0:
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      b0 =
     1 (  ac5 - ac6 +
     2                 ac7  ) / bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c the eigenvaules are just the solution of the quadratic in s:
c so return as "decay"
c  *** Note that the decay times are defined as positive ***
c      if a negative inverse decay time is returned there is an error!
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      sb1 = b1*b1
      fb0 = 4.0d0*b0
      diff =  sb1 - fb0          
      if(diff.le.zero) go to 25
      rs =  dsqrt( diff )          
      decay(1) = -( - b1 - rs ) / 2.0d0
      decay(2) = -( - b1 + rs ) / 2.0d0          
      go to 26
   25 idgen = 100
      go to 9990
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c and for the numerator part of the quadratic s dependence
c the Mathematica session is:
c
c Out[14]=
c                                                                      2
c ac10 + ac11 + (2 ac10 + 2 ac11 + ac9) s + (ac10 + ac11 + ac8 + ac9) s
c----------------------------------------------------------------------
c            ac0 - ac1 + ac2 + ac3 - ac4 + ac5 - ac6 + ac7
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
   26 a0 = ( ac10 + ac11 ) / bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      a1 = ( 2.0d0*( ac10 + ac11 )
     1                             + ac9 )  / bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      a2 =  (ac10 + ac11
     1                   + ac8 + ac9) / bc
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c and the following terms are the amplitudes of the inverse Laplace
c transform solution for the non-q part.  (See the boxed equation on
c page 4 of the Nov. 23 1996 notes.)
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      decdif = 1.0d0/(decay(2) - decay(1))
      amps(1) = -decdif*( decay(1) * ( a1 - a2*decay(1) ) - a0 )
      amps(2) =  decdif*( decay(2) * ( a1 - a2*decay(2) ) - a0 )
      amps(3) = a2
      amps(4) = - decay(1) * amps(1) 
      amps(5) = - decay(2) * amps(2)
      go to 999
 9990 write(6,998) idgen
  998 format(' idgen val ** fatal error ** degenerate e.v.'/1h ,1p,1i12) 
  999 return
      end
