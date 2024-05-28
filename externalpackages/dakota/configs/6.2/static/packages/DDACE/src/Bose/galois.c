#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #include "ddace_config.h"
#endif
/**

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


/*     Manipulation of generic Galois fields.  See gfields.c
   for construction of specific Galois fields.  */

#include <math.h>
#include <stdio.h>

#include "galois.h"

/* Fix for C99 */
extern void free_ivector(int v, int nl, int nh); /* memory.c */
extern void free_imatrix(int m, int nrl, int nrh, int ncl, int nch); /* memory.c */

/*  Glossary:

       GF_poly_sum       Addition in polynomial representation
       GF_poly_prod      Multiplication in polynomial representation
       GF_poly2int       Convert polynomial to integer in 0..q-1
       GF_ready          Prepare (+,*,^-1) lookup tables
       GF_print          Print Galois field
       GF_free           Free storage

*/


/*---------------------------------------------------------------*/

GF_poly_sum( p,n,p1,p2,sum )
int  n,p,*p1,*p2,*sum;
{
int i;

for(  i=0; i<n; i++  )
  sum[i] = (p1[i]+p2[i]) % p;
}

/*---------------------------------------------------------------*/

GF_poly_prod( p,n,xton,p1,p2,prod )
/*
  Set prod = p1*p2 with coefficients modulo p, and x^n replaced
by polynomial xton.
*/
int  n,p,*xton,*p1,*p2,*prod;
{
int i,j, *longprod;

longprod = ivector(0,2*n-2);

for(  i=0; i<2*n-1; i++  )
  longprod[i] = 0;
for(  i=0; i<n; i++  )
for(  j=0; j<n; j++  )
  longprod[i+j] += p1[i]*p2[j];
for(  i=2*n-2; i>n-1; i--  )
for(  j=0; j<n; j++  )
  longprod[i-n+j] += xton[j]*longprod[i];
for(  i=0; i<n; i++  )
  prod[i] = longprod[i] % p;

free_ivector(longprod,0,2*n-2);
}

/*---------------------------------------------------------------*/

int GF_poly2int( p,n,poly )

int n,p,*poly;
{
int ans, i;

ans = 0;
for(  i=n-1; i>0; i--  )
  ans = (ans+poly[i])*p;
ans += poly[0];

return ans;
}

/*---------------------------------------------------------------*/

#define GFPUNT {fprintf(stderr,"Unable to allocate space for Galois field on %d elements.\n",q);return 0;}

GF_ready( gf,p,n,xton )
/* 
   Make ready the Galois Field
*/
struct GF *gf;
int  n,p,*xton;
{
int i,j,q,click,*poly;

poly = ivector(0,n-1);

gf->n = n;
gf->p = p;
q = 1;
for(  i=0; i<n; i++  )
  q *= p;
gf->q = q;
gf->xton = ivector(0,n-1);        if(  !gf->xton  )GFPUNT;
for(  i=0; i<n; i++  )
  gf->xton[i] = xton[i];
gf->plus = imatrix(0,q-1,0,q-1);  if(  !gf->plus  )GFPUNT;
gf->times= imatrix(0,q-1,0,q-1);  if(  !gf->times )GFPUNT;
gf->inv  = ivector(0,q-1);        if(  !gf->inv   )GFPUNT;
gf->neg  = ivector(0,q-1);        if(  !gf->neg   )GFPUNT;
gf->root = ivector(0,q-1);        if(  !gf->root  )GFPUNT;
gf->poly = imatrix(0,q-1,0,n-1);  if(  !gf->poly  )GFPUNT;

for(  i=0; i<n; i++  )
  gf->poly[0][i] = 0;

for( i=1; i<q; i++  ){
  for( click=0; gf->poly[i-1][click]==(p-1); click++  )
    gf->poly[i][click] = 0;
  gf->poly[i][click] = gf->poly[i-1][click]+1;
  for(  j=click+1; j<n; j++  )
    gf->poly[i][j] = gf->poly[i-1][j];
}

for(  i=0; i<q; i++  )
for(  j=0; j<q; j++  ){
  GF_poly_sum( p,n,gf->poly[i],gf->poly[j],poly );
  gf->plus[i][j] = GF_poly2int( p,n,poly );
  GF_poly_prod( p,n,xton,gf->poly[i],gf->poly[j],poly );
  gf->times[i][j] = GF_poly2int( p,n,poly );
}

for(  i=0; i<q; i++  ){
  gf->inv[i] = -1;
  for(  j=0; j<q;  j++  )
    if(  gf->times[i][j]==1  )
      gf->inv[i] = j;
  if(  i>0 && gf->inv[i] <= 0  ){
    fprintf(stderr,"There is something wrong with the Galois field\n");
    fprintf(stderr,"used for q=%d.  Element %d has no reciprocal.\n",q,i);
    return 0;
  }
}

for(  i=0; i<q; i++  ){
  gf->neg[i] = -1;
  for(  j=0; j<q;  j++  )
    if(  gf->plus[i][j]==0  )
      gf->neg[i] = j;
  if(  i>0 && gf->neg[i] <= 0  ){
    fprintf(stderr,"There is something wrong with the Galois field\n");
    fprintf(stderr,"used for q=%d.  Element %d has no negative.\n",q,i);
    return 0;
  }
}

for(  i=0; i<q; i++  ){
  gf->root[i] = -1;
  for(  j=0; j<q;  j++  )
    if(  gf->times[j][j]==i  )
      gf->root[i] = j;
}

/*printf("%d %d %d\n",q,p,n);*/

/* Added by J Cramp 07.02.04 to fix memory leak. */
free_ivector(poly,0,n-1);

return 1;
}

/*---------------------------------------------------------------*/

GF_print( gf )

struct GF *gf;
{
int i,j,n,p,q;

n=gf->n, p=gf->p, q=gf->q;

if( q>999 )fprintf(stderr,"Warning q=%d will overflow print field.\n",q);

printf("\nFor GF(%d) p=%d n=%d\n",q,p,n);
printf("x**n = (");
for( i=0; i<n-1; i++  )
  printf("%d,",gf->xton[i]);
printf("%d)\n",gf->xton[n-1]);
printf("\n\nGF(%d) Polynomial coefficients:\n",q);
for(  i=0; i<q; i++  ){
  printf("  %3d  ",i);
  for(  j=0; j<n; j++  )
    printf("%3d ",gf->poly[i][j]);
  printf("\n");
}
printf("\n\nGF(%d) Addition Table\n",q);
for(  i=0; i<q; i++  ){
  printf("  ");
  for(  j=0; j<q; j++  )
    printf(" %3d",gf->plus[i][j]);
  printf("\n");
}
printf("\n\nGF(%d) Multiplication table\n",q);
for(  i=0; i<q; i++  ){
  printf("  ");
  for(  j=0; j<q; j++  )
    printf(" %3d",gf->times[i][j]);
  printf("\n");
}
printf("\n\nGF(%d) Reciprocals\n",q);
for(  i=1; i<q; i++  )
  printf(" %3d %3d\n",i,gf->inv[i]);

printf("\n\nGF(%d) Negatives\n",q);
for(  i=0; i<q; i++  )
  printf(" %3d %3d\n",i,gf->neg[i]);

printf("\n\nGF(%d) Square roots\n",q);
for(  i=0; i<q; i++  )
  printf(" %3d %3d\n",i,gf->root[i]);
}

/*---------------------------------------------------------------*/

GF_free( gf )
/* 
   Deallocate the Galois Field
*/
struct GF *gf;
{
/*int q,p,n;
q = gf->q, p=gf->p, n=gf->n;*/

int q,n;
q = gf->q, n=gf->n; /* Shuts lint up */

free_imatrix(gf->poly,0,q-1,0,n-1);
free_ivector(gf->root,0,q-1);
free_ivector(gf->neg,0,q-1);
free_ivector(gf->inv,0,q-1);
free_imatrix(gf->times,0,q-1,0,q-1);
free_imatrix(gf->plus,0,q-1,0,q-1);
free_ivector(gf->xton,0,n-1);
}
