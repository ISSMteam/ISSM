#ifndef _INTERPOLATION_H_
#define _INTERPOLATION_H_

#ifdef HAVE_CONFIG_H
   #include <config.h>
#else
#error "Cannot compile with HAVE_CONFIG_H symbol! run configure first!"
#endif

#include "../Exceptions/exceptions.h"
/*findindices {{{*/
template <typename doubleT>
bool findindices(int* pn,int* pm,doubleT* x,int x_rows, doubleT* y,int y_rows, doubleT xgrid, doubleT ygrid){
   /*Find indices m and n into y_grid and x_grid, for which  y_grid(m)<=y<=y_grid(m+1) and x_grid(n)<=x<=x_grid(n+1)
    *    Q12             Q22
    * y2 x---------+-----x
    *    |         |     |
    *    |         |P    |
    *    |---------+-----|
    *    |         |     |
    *    |         |     |
    * y1 x---------+-----x Q21
    *    x1                 x2
    */

   bool foundx=false,foundy=false;
   int m=-1,n=-1;
   int i;

   for (i=0;i<x_rows-1;i++){
      if ((x[i]<=xgrid) && (xgrid<x[i+1])){
         n=i;
         foundx=true;
         break;
      }
   }
   if(xgrid==x[x_rows-1]){
      n=x_rows-2;
      foundx=true;
   }

   for (i=0;i<y_rows-1;i++){
      if ((y[i]<=ygrid) && (ygrid<y[i+1])){
         m=i;
         foundy=true;
         break;
      }
   }
   if(ygrid==y[y_rows-1]){
      m=y_rows-2;
      foundy=true;
   }

   /*Assign output pointers:*/
   *pm=m; *pn=n;
   return (foundx && foundy);
}/*}}}*/
/*bilinearinterp{{{*/
template <typename doubleT> doubleT bilinearinterp(doubleT x1,doubleT x2,doubleT y1,doubleT y2,doubleT Q11,doubleT Q12,doubleT Q21,doubleT Q22,doubleT x,doubleT y){
   /*Bilinear  interpolation: (http://en.wikipedia.org/wiki/Bilinear_interpolation) */

   /*    Q12    R2        Q22
    * y2 x------x---------x
    *    |      |         |
    *    |      |         |
    *    |      +P        |
    *    |      |         |
    *    |Q11   R1        Q21
    * y1 x------x---------x
    *    x1               x2
    *
    */

   /*Checks*/
   _assert_(x2>x1 && y2>y1);
   _assert_(x<=x2 && x>=x1 && y<=y2 && y>=y1);

   return
     +Q11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1))
     +Q21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1))
     +Q12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1))
     +Q22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1));
}
/*}}}*/
#endif
