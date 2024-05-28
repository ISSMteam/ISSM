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


#include <stdio.h>
/* for exit: */
#include <stdlib.h>

#include "oa.h"
  
int  **imatrix(), *ivector();

/* Fix for C99 */
extern int grow_imatrix_byrows(int ***imat, int oldrowsize, int newrowsize, int colsize); /* memory.c */
extern void free_ivector(int *v, int nl, int nh); /* memory.c */
extern int ipow(int a, int b); /* primes.c */

/**

   Glossary:

    OA_put            write OA to standard output
    OA_fput           write OA to stream
    OA_get            get OA from standard input
    OA_fget           get OA from stream

    OA_parsein        read arguments q,nrow,ncol to OA "filter programs"
    OA_strworkcheck   warn about large work loads in strength checking programs

*/


/*  OUTPUT    OUTPUT    OUTPUT    OUTPUT    OUTPUT    OUTPUT    OUTPUT  */

OA_fput( stream,A,n,k,q )
FILE *stream;
int **A, n, k, q;
{
int   i,j;
char* format;

format = "%d%s";
if(  q < 1000  )format = "%3d%s";
if(  q < 100   )format = "%2d%s";
if(  q < 10    )format = "%1d%s";
if(  q < 0     )format = "%d%s";

for(  i=0; i<n; i++  )
for(  j=0; j<k; j++  )
  fprintf(stream,format,A[i][j],(j==k-1)?"\n":" ");
}

OA_put( A,n,k,q )
int **A, n, k, q;
{
OA_fput( stdout,A,n,k,q );
}

/*  INPUT    INPUT    INPUT    INPUT    INPUT    INPUT    INPUT    INPUT  */

OA_fget( stream,A,n,k,q,eof_assert )
FILE *stream;
int **A, n, k, q, eof_assert;
{
int   i,j;
for( i=0; i<n; i++ )
for( j=0; j<k; j++ )
  if( fscanf(stream,"%d",&A[i][j]) == EOF  ){
    fprintf(stderr,"Unexpected end of input encountered.  Wanted to read\n");
    fprintf(stderr,"%d rows of %d cols.  Failed trying for row %d, col %d.\n",
      n,k,i,j);
    return 0;
  }else if(  A[i][j] >= q ){
    fprintf(stderr,"Invalid array element %d.  All elements should be\n",
      A[i][j]);
    fprintf(stderr,"strictly less than q = %d.\n",q);
    return 0;
  }else if(  A[i][j] <0   ){
    fprintf(stderr,"Invalid array element %d, should be >= 0.\n",A[i][j]);
    return 0;
  }

if(  eof_assert  &&  fscanf(stream,"%d",&eof_assert) != EOF  ){
    fprintf(stderr,"Input has more integers than expected.\n");
    fprintf(stderr,"Perhaps the number of rows and/or columns is incorrect.\n");
    return 0;
  }
return 1;
}

OA_get( A,n,k,q, eof_assert )
int **A, n, k, q, eof_assert;
{
return OA_fget( stdin,A,n,k,q,eof_assert );
}

/*  READ    READ    READ    READ    READ    READ    READ    READ    READ  */

/*  Read array and determine dimensions.   Inspired by approach
taken in Xgobi.  */

#define MAXK 5000
#define ROWINC 1000
int line0[ MAXK ];

OA_fread( stream,A,n,k,q )
FILE *stream;
int ***A, *n, *k, *q;
{
int   i, j;
char  c;

*k = 0;

while( (c = getc(stream)) != '\n'  ){
  if(  c=='\t' || c==' '  )
    ;
  else{
    ungetc(c,stream);
    if(  *k >= MAXK  ){
      fprintf(stderr,"Error: Input appears to have more than %d columns.\n",*k);
      fprintf(stderr,"If such large input is desired, increase MAXK in oa.c\n");
      fprintf(stderr,"and re-install the software.\n");
      return 0;
    }
    if(  fscanf(stream,"%d",&line0[(*k)++]) <= 0  ){
      fprintf(stderr,"Error: no newline character found.  Could be empty\n");
      fprintf(stderr,"input or matrix all on one line.\n");
      return 0;
    }
  }
}

*A = imatrix( 0,ROWINC-1, 0,*k-1 );
if(  !(*A)  ){
  fprintf(stderr,"Unable to allocate memory to read the array.\n");
  return 0;
}

*n = 0;
for(  j=0; j<*k; j++  )
  A[0][*n][j] = line0[j];

while( 1 ){
  (*n)++;
  if(  (*n % ROWINC)==0  )
    if(  !grow_imatrix_byrows( A, (*n), (*n)+ROWINC, *k )  ){
      fprintf(stderr,"Unable to allocate extra memory for row %d of the array.\n"
	      ,*n);
      return 0;
    }

  if(  fscanf(stream,"%d",&(A[0][*n][0])) == EOF  )
    break;

  for(  j=1; j<*k; j++  )
    if( fscanf(stream,"%d",&(A[0][*n][j])) == EOF  ){
      fprintf(stderr,"Unexpected end of input encountered.  Row %d only has\n",*n);
      fprintf(stderr,"%d elements of %d expected.\n",j,*k);
      return 0;
    }
}
*q = 0;                   /* Assume that q = max(A)+1  */

for(  i=0; i<*n; i++  )
for(  j=0; j<*k; j++  )
  if(  A[0][i][j] > *q  )
    *q = A[0][i][j];
*q = *q +1;

return 1;
}

OA_read( A,n,k,q )
int ***A, *n, *k, *q;
{
return OA_fread( stdin,A,n,k,q );
}

/*  PARSE    PARSE    PARSE    PARSE    PARSE    PARSE    PARSE    PARSE  */

OA_parsein( argc,argv, q,nrow,ncol, A )
int  argc;
char *argv[];
int *q, *nrow, *ncol, ***A;
{
if(  argc<=1  ){
  if(  !OA_read( A,nrow,ncol,q )  ){
    fprintf(stderr,"Fatal error while reading the array.\n");
    exit(1);
  }
}
else if( argc==2  ){
  sscanf(argv[1],"%d",q);
  scanf("%d %d",nrow,ncol);
}else if( argc==3  ){
  sscanf(argv[1],"%d",q);
  sscanf(argv[2],"%d",nrow);
  scanf("%d",ncol);
}else{
  sscanf(argv[1],"%d",q);
  sscanf(argv[2],"%d",nrow);
  sscanf(argv[3],"%d",ncol);
}

if(  *q<1  ){
  fprintf(stderr,"Array has only %d symbol(s).  At least one\n",*q);
  fprintf(stderr,"symbol is necessary in an orthogonal array.\n");
  exit(1);
}

if(  *ncol<1  ){
  fprintf(stderr,"Array has only %d column(s).  At least one\n",*ncol);
  fprintf(stderr,"column is necessary in an orthogonal array.\n");
  exit(1);
}

if(  *nrow<1  ){
  fprintf(stderr,"Array has only %d rows.  At least one\n",*nrow);
  fprintf(stderr,"row is necessary in an orthogonal array.\n");
  exit(1);
}


if(  argc >1  ){
  *A = imatrix( 0,*nrow-1,0,*ncol-1 );
  if(  !(*A)  ){
    fprintf(stderr,"The array is too large (%d by %d) to fit in memory.\n",nrow,ncol);
    exit(1);
  }
  
  if(  !OA_get( *A,*nrow,*ncol,*q, 1 )  ){/* Read 'em all */
    fprintf(stderr,"Read error getting the orthogonal array.\n");
    exit(1);
  }
}
}


/*  WORK    WORK    WORK    WORK    WORK    WORK    WORK    WORK    WORK  */

OA_strworkcheck( work,str )
double work;
int    str;
{
if(  work > BIGWORK  ){
  fprintf(stderr,"If the array has strength %d, %g comparisons will\n",
	  str, work);
  fprintf(stderr,"be required to prove it.  This might take a long time.\n");
  fprintf(stderr,"This warning is triggered when more than %d comparisons\n",
	  BIGWORK);
  fprintf(stderr,"are required.  To avoid this warning increase BIGWORK in\n");
  fprintf(stderr,"oa.h.  Intermediate results will be printed.\n\n");
  fflush(stderr);
}else if(  work > MEDWORK  ){
  fprintf(stderr,"Since more than %d comparisons may be required to\n",MEDWORK);
  fprintf(stderr,"to check whether the array has strength %d, intermediate\n",
	  str);
  fprintf(stderr,"results will be printed.  To avoid this warning increase\n");
  fprintf(stderr,"MEDWORK in oa.h\n\n");
  fflush(stderr);
}
}


/*  STRENGTH    STRENGTH    STRENGTH    STRENGTH    STRENGTH  */

/* Check strength 0 */
OA_str0( q,nrow,ncol,A,verbose   )
int      q,nrow,ncol,**A, verbose;
{
int  i, j1;

for(  j1=0; j1<ncol; j1++  )
for(  i=0; i<nrow; i++  )
  if(  A[i][j1] < 0  || A[i][j1] >= q  ){
    if(  verbose >= 2  ){
      printf("Array is not even of strength 0, that is there are elements\n");
      // BMA 20141204; added q-1 as argument to silence compiler warning
      //printf("other than integers 0 through %d inclusive in it.\n");
      printf("other than integers 0 through %d inclusive in it.\n", (q-1));
      printf("The first exception is A[%d][%d] = %d.\n",i,j1,A[i][j1]);
    }
    return 0;
  }
if(  verbose >=2  )
  printf("The array has strength (at least) 0.\n");
return 1;
}


/* Check strength 1 */
OA_str1( q,nrow,ncol,A,verbose   )
int      q,nrow,ncol,**A, verbose;
{
int     i, j1, q1;
int     lambda, count;
double  work;

if(  nrow%q  ){
  if(  verbose >= 2  ){
    printf("The array cannot have strength 1, because the number\n");
    printf("of rows %d is not a multiple of q = %d.\n",nrow,q);
  }
  return 0;
}

lambda = nrow/q;
work = nrow * ncol * q * 1.0;
OA_strworkcheck( work,1 );
for(  j1=0; j1<ncol; j1++  ){
  for(  q1=0; q1<q; q1++  ){
    count = 0;
    for( i=0; i<nrow; i++  )
      count += (A[i][j1]==q1);
    if(  count != lambda   ){
      if(  verbose >= 2  ){
	printf("Array is not of strength 1.  The first violation arises for\n");
	printf("the number of times A[,%d] = %d.\n",
	       j1, q1);
	printf("This happened in %d rows, it should have happened in %d rows.\n",
	       count, lambda);
      }
      return 0;
    }
  }
if(  work > MEDWORK && verbose > 0  )
  fprintf(stderr,"No violation of strength 1 involves column %d.\n",j1);
}
if(  verbose >=2  )
  printf("The array has strength (at least) 1.\n");
return 1;
}

/* Check strength 2  */
OA_str2( q,nrow,ncol,A,verbose   )
int      q,nrow,ncol,**A, verbose;
{
int  i, j1,j2, q1,q2;
int  lambda, count;
double  work;

if(  ncol<2  ){
  if(  verbose > 0 ){
    fprintf(stderr,"Array has only %d column(s).  At least two\n",ncol);
    fprintf(stderr,"columns are necessary for strength 2 to make sense.\n");
  }
  return 0;
}
if(  nrow % (q*q)  ){
  if(  verbose > 0 ){
    fprintf(stderr,"The array cannot have strength 2, because the number\n");
    fprintf(stderr,"of rows %d is not a multiple of q^2 = %d^2 = %d.\n",nrow,q,q*q);
  }
  return 0;
}

lambda = nrow/(q*q);
work = nrow*ncol*(ncol-1.0)*q*q/2.0;
OA_strworkcheck( work,2 );

for(  j1=0;    j1<ncol; j1++  ){
for(  j2=j1+1; j2<ncol; j2++  ){
  for(  q1=0; q1<q; q1++  )
  for(  q2=0; q2<q; q2++  ){
    count = 0;
    for( i=0; i<nrow; i++  )
      count += (A[i][j1]==q1)&&(A[i][j2]==q2);
    if(  count != lambda   ){
      if(  verbose >= 2 ){
	printf("Array is not of strength 2.  The first violation arises for\n");
	printf("the number of times (A[,%d],A[,%d]) = (%d,%d).\n",
	       j1,j2, q1,q2 );
	printf("This happened in %d rows, it should have happened in %d rows.\n",
	       count, lambda);
      }
      return 0;
    }
  }
}
if(  work > MEDWORK && verbose > 0 )
  fprintf(stderr,"No violation of strength 2 involves column %d.\n",j1);
}

if(  verbose >=2  )
  printf("The array has strength (at least) 2.\n");
return 1;
}



/* Check strength 3  */
OA_str3( q,nrow,ncol,A,verbose   )
int      q,nrow,ncol,**A, verbose;
{
int  i, j1,j2,j3, q1,q2,q3;
int  lambda, count;
double  work;

if(  ncol<3  ){
  if(  verbose > 0 ){
    fprintf(stderr,"Array has only %d column(s).  At least three\n",ncol);
    fprintf(stderr,"columns are necessary for strength 3 to make sense.\n");
  }
  return 0;
}
if(  nrow % (q*q*q)  ){
  if(  verbose > 0 ){
    fprintf(stderr,"The array cannot have strength 3, because the number\n");
    fprintf(stderr,"of rows %d is not a multiple of q^3 = %d^3 = %d.\n",nrow,q,q*q*q);
  }
  return 0;
}

lambda = nrow/(q*q*q);
work = nrow*ncol*(ncol-1.0)*(ncol-2.0)*q*q*q/6.0;
OA_strworkcheck( work,3 );

for(  j1=0;    j1<ncol; j1++  ){
for(  j2=j1+1; j2<ncol; j2++  )
for(  j3=j2+1; j3<ncol; j3++  ){
  for(  q1=0; q1<q; q1++  )
  for(  q2=0; q2<q; q2++  )
  for(  q3=0; q3<q; q3++  ){
    count = 0;
    for( i=0; i<nrow; i++  )
      count += (A[i][j1]==q1)&&(A[i][j2]==q2)&&(A[i][j3]==q3);
    if(  count != lambda   ){
      if(  verbose >= 2 ){
	printf("Array is not of strength 3.  The first violation arises for\n");
	printf("the number of times (A[,%d],A[,%d],A[,%d]) = (%d,%d,%d).\n",
	       j1,j2,j3,  q1,q2,q3 );
	printf("This happened in %d rows, it should have happened in %d rows.\n",
	       count, lambda);
      }
      return 0;
    }
  }
}
if(  work > MEDWORK && verbose > 0 )
  fprintf(stderr,"No violation of strength 3 involves column %d.\n",j1);
}
if(  verbose >=2  )
  printf("The array has strength (at least) 3.\n");
return 1;
}


/* Check strength 4  */
OA_str4( q,nrow,ncol,A,verbose   )
int      q,nrow,ncol,**A, verbose;
{
int  i, j1,j2,j3,j4, q1,q2,q3,q4;
int  lambda, count;
double  work;

if(  ncol<4  ){
  if(  verbose > 0 ){
    fprintf(stderr,"Array has only %d column(s).  At least four\n",ncol);
    fprintf(stderr,"columns are necessary for strength 4 to make sense.\n");
  }
  return 0;
}
if(  nrow % (q*q*q*q)  ){
  if(  verbose > 0 ){
    fprintf(stderr,"The array cannot have strength 4, because the number\n");
    fprintf(stderr,"of rows %d is not a multiple of q^4 = %d^4 = %d.\n",nrow,q,q*q*q*q);
  }
  return 0;
}

lambda = nrow/(q*q*q*q);
work = nrow*ncol*(ncol-1.0)*(ncol-2.0)*(ncol-3.0)*q*q*q*q/24.0;
OA_strworkcheck( work,4 );

for(  j1=0;    j1<ncol; j1++  ){
for(  j2=j1+1; j2<ncol; j2++  )
for(  j3=j2+1; j3<ncol; j3++  )
for(  j4=j3+1; j4<ncol; j4++  ){
  for(  q1=0; q1<q; q1++  )
  for(  q2=0; q2<q; q2++  )
  for(  q3=0; q3<q; q3++  )
  for(  q4=0; q4<q; q4++  ){
    count = 0;
    for( i=0; i<nrow; i++  )
      count += (A[i][j1]==q1)&&(A[i][j2]==q2)&&(A[i][j3]==q3)&&(A[i][j4]==q4);
    if(  count != lambda  ){
      if(  verbose >= 2  ){
	printf("Array is not of strength 4.  The first violation arises for\n");
	printf("the number of times (A[,%d],A[,%d],A[,%d],A[,%d]) = (%d,%d,%d,%d).\n",
	       j1,j2,j3,j4, q1,q2,q3,q4 );
	printf("This happened in %d rows, it should have happened in %d rows.\n",
	       count, lambda);
      }
      return 0;
    }
    }
  }
}
if(  work > MEDWORK && verbose > 0 )
  fprintf(stderr,"No violation of strength 4 involves column %d.\n",j1);

if(  verbose >=2  )
  printf("The array has strength (at least) 4.\n");
return 1;
}


/* Check strength t  */
OA_strt( q,nrow,ncol,A,t,verbose   )
int      q,nrow,ncol,**A,t,verbose;
{
int  row, i, ic, iq, *clist, *qlist, ctuples, qtuples;
int  lambda, count, match;
double  work;

if(  t<0  ){
  if(  verbose > 0 ){
    fprintf(stderr,"Don't know how to verify strength %d.  It doesn't\n",t);
    fprintf(stderr,"make sense.\n");
  }
  return 0;
}
if(  ncol<t  ){
  if(  verbose > 0 ){
    fprintf(stderr,"Array has only %d column(s).  At least %d\n",ncol,t);
    fprintf(stderr,"columns are necessary for strength %d to make sense.\n",t);
  }
  return 0;
}
if(  t==0  )
  return OA_str0( q,nrow,ncol,A,verbose );
if(  nrow % ipow(q,t)  ){
  if(  verbose > 0 ){
    fprintf(stderr,"The array cannot have strength %d, because the number\n",t);
    fprintf(stderr,"of rows %d is not a multiple of q^%d = %d^%d = %d.\n",
	    nrow,t,q,t,ipow(q,t));
  }
  return 0;
}

lambda = nrow/ipow(q,t);
work   = nrow*ipow(q,t);
ctuples = 1;

clist = ivector( 0,t-1 );
qlist = ivector( 0,t-1 );

for(  i=0; i<t; i++  ){
  work *= (ncol-i)/(i+1.0);
  ctuples *= ncol-i;
  qlist[i] = 0;
  clist[i] = i;
}
for(  i=0; i<t; i++  )
  ctuples /= (i+1);
qtuples = ipow(q,t);

OA_strworkcheck( work,t );

for(  ic=0; ic<ctuples; ic++  ){   /* Loop over ordered tuples of columns */
/*  for( i=0; i<t; i++  )
    printf("%s %d%s",(i==0)?"Col":"",clist[i],(i==t-1)?"\n":" ");*/

  for(  iq=0; iq<qtuples;   iq++  ){ /* Loop over unordered tuples of symbols */
/*    for( i=0; i<t; i++  )
      printf("  %s %d%s",(i==0)?"Sym":"",qlist[i],(i==t-1)?"\n":" ");*/
    count = 0;
    for( row=0; row<nrow; row++  ){
      match = 1;
      for(  i=0; i<t && match; i++  )
	match *= A[row][clist[i]] == qlist[i];
      count += match;
    }
    if(  count != lambda  ){
      if(  verbose >= 2  ){
	printf("Array is not of strength %d.  The first violation arises for\n",t);
	printf("the number of times (");
	for(  i=0; i<t; i++  )
	  printf("A[,%d]%s",clist[i],(i==t-1)?")":",");
	printf(" = (");
	for(  i=0; i<t; i++  )
	  printf("%d%s",qlist[i],(i==t-1)?").\n":",");
	printf("This happened in %d rows, it should have happened in %d rows.\n",
	       count, lambda);
      }
      return 0;
    }
    for( i=t-1; i>=0; i--  ){
      qlist[i] = (qlist[i]+1) % q;
      if(  qlist[i]  )break;
    }
  }

  for( i= t-1; i>=0; i--  ){
    clist[i] = (clist[i]+1) % (ncol+i-t+1);
    if(  clist[i]  )break;
  }

  if(  work > MEDWORK && verbose > 0 && (t==1||t>1 && clist[1]==0)  )
    fprintf(stderr,"No violation of strength %d involves column %d.\n",
	    t,(clist[0]+ncol-1)%ncol);

  for( i=1; i< t; i++  )
    if(  clist[i] <= clist[i-1]  )
      clist[i] = clist[i-1]+1;
}

if(  verbose >=2  )
  printf("The array has strength (at least) %d.\n",t);

/* Added by J Cramp on 07.02.04 to
 * fix memory leak problems. */
free_ivector(clist,0,t-1);
free_ivector(qlist,0,t-1);

return 1;
}

void OA_strength( q,nrow,ncol,A,str,verbose )
int  q,nrow,ncol,**A,*str, verbose;
/*
     Calculate and return the strength of the array A.

verbose = 0   =>   No printed output
verbose = 1   =>   Only stderr output
verbose = 2   =>   Output to both stdout and stderr

*/
{
*str = -1;

if(  OA_str0( q,nrow,ncol,A,verbose)   )
  *str = 0;
else
  return;
if(  OA_str1( q,nrow,ncol,A,verbose)   )
  *str = 1;
else
  return;
while( OA_strt( q,nrow,ncol,A,*str+1,verbose )  )
  (*str)++;
return;
}
