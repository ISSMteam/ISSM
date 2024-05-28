#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #include "ddace_config.h"
#endif
/*

  These programs construct and manipulate orthogonal 
arrays.  They were prepared by

    Art Owen
    Department of Statistics
    Sequoia Hall
    Stanford CA 94305

  They may be freely used and shared.  This code comes
with no warranty of any kind.  Use it at your own
risk.

  I thank the Semiconductor Research Corporation and
the National Science Foundation for supporting this
work.

*/


/*     Construction of specific Galois fields.  See galois.c 
    for manipulation of generic Galois fields.

       To add a new Galois field GF(q)=GF(p^n) it is necessary 
    to find the polynomial representation of x^n, and use it to
    include another *xtn--- vector like those used below.

*/

#include <math.h>
#include <stdio.h>

#include "galois.h"

/* Fix for C99 */
extern int isprime(int p); /* primes.c */
extern void primepow(int q, int *p, int *n, int *isit); /* primes.c */
extern int GF_ready(struct GF *gf, int p, int n, int *xton); /* galois.c */

/*  Glossary:

       xtn2t2              Polynomial representation of x^2 in gf2t2
       xtn2t3              Polynomial representation of x^3 in gf2t3
       . . .               Polynomial representation of p^n in gfptn
       xtn13t2             Polynomial representation of x^2 in gf13t2

       GF_set_fields       Initialize polynomial representations
       GF_fields_are_set   Indicates initialization done
       GF_getfield         Construct and return GF(q) if possible
*/

#include "xtndeclare.h"
int *xtnpt1;

int GF_fields_are_set = 0;

GF_set_fields()
{
/* Brute force set up of defining vectors, from Carmichael */

/* Declare x-to-the-power-n vectors, for GFs p-to-the-n */

if(   GF_fields_are_set   )
  fprintf(stderr,"Warning: Fields being re-initialized.  Possible memory waste.\n");


#include "xtnset.h"

xtnpt1 = ivector(0,0);
xtnpt1[0] = 0; 
GF_fields_are_set = 1;
}


GF_getfield( q, gf )
int q;
struct GF *gf;
{
int *xtn;
int p,n,ispp;

if(  !GF_fields_are_set  )
  GF_set_fields();

if(  q<1  ){      /* Impossible argument */
  fprintf(stderr,"Field must have positive number of elements.\n");
  return 0; }
if(  q==1 ){      /* Pointless  argument */
  fprintf(stderr,"Field with 1 element was requested.\n");
  return 0; }

primepow( q,&p,&n,&ispp  );
if(  !ispp  ){
  fprintf(stderr,"q=%d is not a prime power.\n",q);
  return 0; }

xtn = NULL;

#include "xtndispatch.h"

if(  isprime(q)  )xtn = xtnpt1;  /* Could have tested p=q, or n=1 */

if(  xtn   ){
  if(  GF_ready( gf,p,n,xtn )  )
    return 1;
  else{
    fprintf(stderr,"Construction failed for GF(%d).\n",q);
    return 0;
  }
}
else {
  fprintf(stderr,"GF(%d) = GF(%d^%d) is not included in this program.\n",q,p,n);
  fprintf(stderr,"To add it, consider modifying gfields.c.\n");
  return 0;
}
}
