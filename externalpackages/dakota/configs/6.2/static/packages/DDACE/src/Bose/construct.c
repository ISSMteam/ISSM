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


/*  Constructions for designs using Galois fields */

#include <math.h>
#include <stdio.h>

#include "galois.h"

/* Fix for C99 */
extern int ipow(int a, int b); /* primes.c */
extern int isprime(int p); /* primes.c */
extern void free_ivector(int *v, int nl, int nh); /* memory.c */
extern void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch); /* memory.c */
extern int akodd(struct GF *gf, int *kay, int *b, int *c, int *k); /* akconst.c */
extern int akeven(struct GF *gf, int *kay, int *b, int *c, int *k); /* akconst.c */

/*  Glossary:

    bose:            OA( q^2, q+1, q, 2  )
                     R.C. Bose (1938) Sankhya Vol 3 pp 323-338
    bosecheck:       test input to bose

    bosebush:        OA( 2q^2, 2q+1, q, 2 ), only implemented for q=2^n
    bosebushcheck:   test input to bosebush

    polyeval:        evaluate a polynomial with coefficients, argument
                     and result in a Galois field

*/


bosecheck( q,ncol )
int q, ncol;
{
if(  ncol > q+1  ){
  fprintf(stderr,"Bose's design must have ncol <= q+1.\n");
  fprintf(stderr,"Had q=%d and ncol=%d.\n",q,ncol);
  return 0;
}
if(  ncol <= 0  ){
  fprintf(stderr,"Nonpositive number of columns requested for Bose's design\n");
  return 0;
}
return 1;
}


bose( gf, A, ncol )
struct GF *gf;
int    **A, ncol;
{
int i,j, icol,  q=gf->q, irow;

if(  !bosecheck(q,ncol)  )return 0;

irow = 0;
for(  i=0; i<q; i++  )
for(  j=0; j<q; j++  ){
  icol=0;
  A[ irow ][ icol++ ] = i;
  if(  ncol > 1 )A[ irow  ][ icol++ ] = j;
  for(  icol=2; icol<ncol; icol++ )
    A[ irow ][ icol ] = gf->plus[j][gf->times[i][icol-1]];
  irow++;
}


return 1;
}

itopoly( n,q,d,coef )
int  n,q,d,*coef;
{
int i;

for(  i=0; i<=d; i++  ){
  coef[i] = n % q;
  n = n/q;
}
}


polyeval( gf, d, poly, arg, value )
/*  find  value = poly(arg) where poly is a polynomial of degree d  
    and all the arithmetic takes place in the given Galois field.*/

struct GF *gf;
int    d, *poly, arg, *value;
{
int   i,ans;

ans = 0;
for(  i= d; i>=0; i--  )  /* Horner's rule */
  ans = gf->plus[  gf->times[ans][arg]  ][  poly[i] ];

*value = ans;
}

bushcheck(q,str,ncol)
int q,str,ncol;
{
if(  ncol > q+1  ){
  fprintf(stderr,"Bush designs require ncol <= q+1.\n");
  fprintf(stderr,"Cannot have q = %d and ncol = %d.\n",q,ncol);
  return 0;
}
if(  str > ncol  ){
  fprintf(stderr,"It doesn't make sense to have an array of strength\n");
  fprintf(stderr,"%d with only %d columns.\n",str,ncol);
  return 0;
}
if(  str >= q+1  ){
  fprintf(stderr,"Bush's (1952) theorem has a condition t<q where t\n");
  fprintf(stderr,"is the strength of the array and q is the number of symbols.\n");
  fprintf(stderr,"Here we have t = %d and q = %d.  The array may still\n",str,q);
  fprintf(stderr,"be useful, but a full factorial would have at least as\n");
  fprintf(stderr,"many columns.\n");
  return 1;
}

return 1;
}


bush( gf, A, str, ncol  )
struct GF *gf;
int       **A, str, ncol;
{
int   *coef;
int   q, i,j;

q = gf->q;
if(  !bushcheck(q,str,ncol)  )return 0;

coef = ivector( 0,str-1  );  
if( !coef  ){ /* Very unlikely */
  fprintf(stderr,"Could not allocate memory for Bush design.\n");
  return 0;
}

for(  i=0; i<ipow(q,str); i++  ){
  itopoly( i,q,str-1,coef );
  A[i][0] = coef[str-1];  
  for(  j=0; j<ncol-1; j++  )
    polyeval( gf, str-1, coef, j, &A[i][1+j] );
}
free_ivector( coef,0,str-1  );  
return 1;
}


addelkempcheck( q,p,ncol )
int  q,p,ncol;
{

if(  p==2 && q>4 ){
  fprintf(stderr,"This Addelman-Kempthorne OA(2q^2,ncol,q,2) is only\n");
  fprintf(stderr,"available for odd prime powers q and for even prime\n");
  fprintf(stderr,"powers q<=4.  q=%d is not available, but a Bose Bush\n",q);
  fprintf(stderr,"construction exists for that design.\n");
  return 0;
}

if(  ncol > 2*q+1  ){
  fprintf(stderr,"The Addelman-Kempthorne construction needs ncol <= 2q+1.\n");
  fprintf(stderr,"Can't have ncol = %d with q = %d,\n",ncol,q);
  return 0;
}

if(  ncol == 2*q+1  ){
  fprintf(stderr,"\nWarning: The Addelman-Kempthorne construction with ncol = 2q+1\n");
  fprintf(stderr,"has a defect.  While it is still an OA(2q^2,2q+1,q,2),\n");
  fprintf(stderr,"there exist some pairs of rows that agree in three columns.\n");
  fprintf(stderr,"The final column in the array is involved in all of these\n");
  fprintf(stderr,"triple coincidences.\n\n\n");
}
return 1;
}


addelkemp( gf, A, ncol )
/* Implement Addelman and Kempthorne's 1961 A.M.S. method with n=2 */
struct GF *gf;
int    ncol, **A;
{
int i,j,m,p,q;
int kay,*b,*c,*k;  /* A&K notation */
int row, col, square, ksquare, temp;

p=gf->p, q=gf->q;

if(  !addelkempcheck( q,p,ncol )  )return 0;

b = ivector( 0,q-1 );
c = ivector( 0,q-1 );
k = ivector( 0,q-1 );

for(  i=0; i<q; i++  ){           /* First q*q rows */
  square = gf->times[i][i];
  for(  j=0; j<q; j++  ){
    row = i*q+j;
    col = 0;
    if( col<ncol  )A[row][col++]=j;
    for(  m=1; m<q && col<ncol; m++  )
      A[row][col++] = gf->plus[i][gf->times[m][j]];
    for(  m=0; m<q && col<ncol; m++  ){
      temp = gf->plus[j][gf->times[m][i]];
      A[row][col++] = gf->plus[temp][square]; /* Rgt cols */
    }
    if( col<ncol  )A[row][col++]=i;
  }
}

if(  p !=2  )                    /* Constants kay,b,c,k for odd p */
  akodd(  gf,&kay,b,c,k );
else                             /* Constants kay,b,c,k for even p */
  akeven( gf,&kay,b,c,k );

for(  i=0; i<q; i++  ){           /* Second q*q rows */
  square = gf->times[i][i];
  ksquare = gf->times[kay][square];
  for(  j=0; j<q; j++  ){
    row = q*q+i*q+j;
    col = 0;
    if( col<ncol  )A[row][col++]=j;
    for(  m=1; m<q && col<ncol; m++,col++  )
      A[row][col] = gf->plus[A[row-q*q][col]][b[m]];
    if( col<ncol  )A[row][col++] = gf->plus[ksquare][j]; /* q+1 */
    for(  m=1; m<q && col<ncol; m++  ){
      temp = gf->times[i][k[m]];
      temp = gf->plus[ksquare][temp];
      temp = gf->plus[j][temp];
      A[row][col++] = gf->plus[temp][c[m]];
    }
    if( col<ncol  )A[row][col++]=i;
  }
}

/*for(  i=0; i<2*q*q; i++  )
for(  j=0; j<ncol; j++  )
  printf("%3d%s",A[i][j],j==(ncol-1)?"\n":" ");
*/

return 1;
}


bosebushcheck( q,p,ncol  )
int  q,p,ncol;
{

if(  p!=2  ){
  fprintf(stderr,"This version of Bose and Bush needs q=2^n for some n.\n");
  return 0;
}

if(  ncol > 2*q+1  ){
  fprintf(stderr,"The Bose-Bush construction needs ncol <= 2q+1.\n");
  fprintf(stderr,"Can't have ncol = %d with q = %d,\n",ncol,q);
  return 0;
}

if(  ncol == 2*q+1  ){
  fprintf(stderr,"\nWarning: The Bose-Bush construction with ncol = 2q+1\n");
  fprintf(stderr,"has a defect.  While it is still an OA(2q^2,2q+1,q,2),\n");
  fprintf(stderr,"there exist some pairs of rows that agree in three columns.\n\n\n");
}
return 1;
}

bosebush( gf, B, ncol )
/* Implement Bose and Bush's 1952 A.M.S. method with p=2, u=1 */
struct GF *gf;
int **B;
{
int p,q,s,irow;
int i,j,k,mul;
int **A;

p=gf->p,   /* GF(q) used to generate design with q/2 levels */
q=gf->q;
s=q/2;     /* number of levels in design */

if(  !bosebushcheck( s,p,ncol )  )
  return 0;

A = imatrix(0,s-1,0,q-1);
if(  !A  ){
  fprintf(stderr,"Unable to allocate scratch space for Bose-Bush array.\n");
  return 0;
}

irow = 0;
for(  i=0; i<q; i++  ){
  for(  j=0; j<q; j++  ){
    mul = gf->times[i][j];
    mul = mul % s;
    for( k=0; k<s; k++  )
/*      A[k][j] = gf->plus[mul][k];*/
      A[k][j] = gf->plus[mul][k];
  }
  for(  k=0; k<s; k++  ){
    for( j=0; j<ncol && j<2*s+1; j++ )
      B[irow][j] = A[k][j];
    if(  ncol==2*s+1  )
      B[irow][ncol-1] = i%s;
    irow++;
  }
}
free_imatrix(A,0,s-1,0,q-1);
return 1;
}
  
  
bosebushlcheck( s,p,lam,ncol  )
int  s,p,lam,ncol;
{

if(  !isprime(p)  ){
  fprintf(stderr,"Bose Bush routine given a nonprime.\n");
  return 0;
}

if(  ncol > lam*s+1  ){
  fprintf(stderr,"The Bose-Bush construction needs ncol <= lambda*q+1.\n");
  fprintf(stderr,"Can't have ncol = %d with lam = %d and q = %d,\n",ncol,lam,s);
  return 0;
}

if(  ncol == lam*s+1  ){
  fprintf(stderr,"\nWarning: The Bose-Bush construction with ncol = lambda*q+1\n");
  fprintf(stderr,"has a defect.  While it is still an OA(lambda*q^2,lambda*q+1,q,2),\n");
  fprintf(stderr,"it may have worse coincidence properties than\n");
  fprintf(stderr,"OA(lambda*q^2,lambda*q+1,q,2).\n");
}
return 1;
}

bosebushl( gf, lam, B, ncol )
/* Implement Bose and Bush's 1952 A.M.S. method with given lambda */
struct GF *gf;
int **B, lam;
{
int p,q,s,irow;
int i,j,k,mul;
int **A;

p=gf->p,   /* GF(q) used to generate design with q/lam levels */
q=gf->q;
s=q/lam;     /* number of levels in design */

if(  !bosebushlcheck( s,p,lam,ncol )  )
  return 0;

A = imatrix(0,s-1,0,q-1);
if(  !A  ){
  fprintf(stderr,"Unable to allocate scratch space for Bose-Bush array.\n");
  return 0;
}

irow = 0;
for(  i=0; i<q; i++  ){
  for(  j=0; j<q; j++  ){
    mul = gf->times[i][j];
    mul = mul % s;
    for( k=0; k<s; k++  )
/*      A[k][j] = gf->plus[mul][k];*/
      A[k][j] = gf->plus[mul][k];
  }
  for(  k=0; k<s; k++  ){
    for( j=0; j<ncol && j<lam*s+1; j++ )
      B[irow][j] = A[k][j];
    if(  ncol==lam*s+1  )
      B[irow][ncol-1] = i%s;
    irow++;
  }
}





free_imatrix(A,0,s-1,0,q-1);
return 1;
}
  
  
