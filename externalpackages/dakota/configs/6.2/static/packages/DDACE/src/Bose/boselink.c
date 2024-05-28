#if defined(HAVE_CONFIG_H) && !defined(DISABLE_DAKOTA_CONFIG_H)
  #include "ddace_config.h"
#endif
#include <math.h>
#include <stdio.h>
#include "galois.h"

/* Fix for C99*/
extern int bose(struct GF *gf, int **A, int ncol); /* construct.c */
extern int bush(struct GF *gf, int **A, int str, int ncol); /* construct.c */
extern void GF_free(struct GF *gf); /* galois.c */
extern int GF_getfield(int q, struct GF *gf); /* gfields.c */

int bose_link( int n, int ninputs, int str, int ***AA )
{
   int    i, nsym, nsamples, **A; 
   double dtemp, dpower;
   struct GF gf;

   if ( n <= 0 ) return -1;

   /* only one method for strength of 3 */

   if ( str == 3 )
   {
      dtemp = (double) n;
      dtemp = pow(dtemp, 0.333333334);
      nsym  = (int) dtemp;
      if( ninputs > (nsym+1) ){
         fprintf(stderr,"Only q+1 = %d cols given in Bush design.\n",nsym+1);
         fprintf(stderr,"Columns requested was %d.\n",ninputs);
	 return -1;
      }
      if(  !GF_getfield(nsym, &gf)  ) {
         fprintf(stderr,"Could not construct the Galois field needed\n");
         fprintf(stderr,"for the strength 3 Bush design.\n");
	 GF_free( &gf );
	 return -1;
      }
      nsamples = nsym * nsym * nsym;
      A = imatrix( 0, nsamples, 0, ninputs-1  );
      if(  !A  )
      {
         fprintf(stderr,"Could not allocate array for Bush design.\n");
	 GF_free( &gf );
	 return -1;
      }
      if(  ! bush( &gf, A, 3, ninputs )  )
      {
         fprintf(stderr,"Unable to construct the strength 3 ");
         fprintf(stderr,"Bush design nsym=%d, ninputs=%d.\n", nsym,ninputs);
	 GF_free( &gf );
	 return -1;
      }
      else
      {
         (*AA) = A;
         return nsamples;
      }
   }

   /* only one method for strength greater than 3 */

   if ( str > 3 )
   {
      dtemp = (double) n;
      dpower = 1.0 / (double) str + 0.00000001;
      dtemp = pow(dtemp, dpower);
      nsym  = (int) dtemp;
      if( ninputs > (nsym+1) ){
         fprintf(stderr,"Only q+1 = %d cols given in Bush design.\n",nsym+1);
         fprintf(stderr,"Columns requested was %d.\n",ninputs);
	 return -1;
      }
      nsamples = nsym;
      for ( i = 1; i < str; i++ ) nsamples = nsamples * nsym;
      if(  !GF_getfield(nsym, &gf)  ){
         fprintf(stderr,"Could not construct the Galois field needed\n");
         fprintf(stderr,"for the strength %d Bush design\n",str);
         fprintf(stderr,"on %d levels.\n", nsym);
	 GF_free( &gf );
	 return -1;
      }
      A = imatrix( 0, nsamples-1, 0, ninputs-1  );
      if(  !A  )
      {
         fprintf(stderr,"Could not allocate array for Bush design.\n");
	 GF_free( &gf );
	 return -1;
      } 
      if(  ! bush( &gf, A, str, ninputs )  )
      {
         fprintf(stderr,"Unable to construct the strength %d \n", str);
         fprintf(stderr,"Bush design nsym=%d, ninputs=%d.\n", nsym,ninputs);
	 GF_free( &gf );
	 return -1;
      }
      else
      {
         (*AA) = A;
	 GF_free( &gf );
         return nsamples;
      }
   }

   /* for the cases when strength = 2 */

   if ( str == 2 )
   {
      dtemp = (double) n;
      dtemp = pow(dtemp, 0.500001);
      nsym  = (int) dtemp;

      if ( ninputs > ( nsym + 1 ) ) 
      {
         fprintf(stderr,"Number of samples too small to construct OA.\n");
         fprintf(stderr,"Need at least %d. \n", (ninputs-1)*(ninputs-1));
         return -1;
      } 
      if(  !GF_getfield( nsym, &gf)  )
      {
         fprintf(stderr,"Could not construct Galois field needed\n");
         fprintf(stderr,"for Bose design.\n");
	 GF_free( &gf );
	 return -1;
      }
      nsamples = nsym * nsym;
      A = imatrix( 0, nsamples-1, 0, ninputs-1  );
      if(  !A  )
      {
         fprintf(stderr,"Could not allocate array for Bose design.\n");
	 GF_free( &gf );
	 return -1;
      }
      if(  ! bose( &gf, A, ninputs )  )
      {
         fprintf(stderr,"Unable to construct Bose design q=%d,",nsym);
         fprintf(stderr," ninputs=%d.\n", ninputs);
	 GF_free( &gf );
	 return -1;
      }
      else
      {
         (*AA) = A;
	 GF_free( &gf );
         return nsamples;
      }
   }
}

