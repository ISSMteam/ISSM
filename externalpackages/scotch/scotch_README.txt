8/06/09:

following INSTALL.txt:
- GNU Make 3.81, lex, yacc all present


[jschierm@astrid src]$ diff Makefile_old.inc Makefile.inc

- started with Makefile.inc.i686_sun_solaris5
- changed CCS from cc to gcc
- removed SCOTCH_PTHREAD, since MPICH2 1.0.2p1 is prior to 1.0.7
- added -std=c99 to CFLAGS for "restrict" attribute (based on various internet pages)
- added -DCOMMON_TIMING_OLD to CFLAGS for undeclared CLOCK_REALTIME (based on common2.c)
- added MEX, CCM, and MFLAGS for Matlab mex modules (as well as Matlab libraries and -largeArrayDims)
- removed -DCOMMON_PTHREAD from CFLAGS and MFLAGS to eliminate fatal exception at pthread_exit in gmap_mex.c
- removed -DCOMMON_FILE_COMPRESS_GZ from CFLAGS and MFLAGS since gzip not used
- added -Wl,-rpath-link,/usr/local/pkgs/matlab-7.6/bin/glnxa64 to LDFLAGS to eliminate missing libhdf5.so.0 in Matlab directories (based on Cielo make)

3a4
> MEX   = .mexa64
12,13c13,16
< CFLAGS        = -m64 -O3 -std=c99 -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DINTSIZE64 -Dintptr_t="long int" -DCOMMON_TIMING_OLD
< LDFLAGS       = -lz -lm -lrt
---
> CCM   = mex
> CFLAGS        = -m64 -O3 -std=c99 -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DCOMMON_TIMING_OLD -DMATLAB -fPIC -I/usr/local/pkgs/matlab-7.6/extern/include
> LDFLAGS       = -lz -lm -lrt -Wl,-rpath-link,/usr/local/pkgs/matlab-7.6/bin/glnxa64 -L/usr/local/pkgs/matlab-7.6/bin/glnxa64 -lmex -lmat
> MFLAGS        = -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -Du_int32_t=uint32_t -Du_int64_t=uint64_t -DCOMMON_TIMING_OLD -DMATLAB -I/usr/local/pkgs/matlab-7.6/extern/include -largeArrayDims


[jschierm@astrid libscotch]$ diff common_old.c common.c

- redirect usagePrint to within Matlab
102a103
> #ifndef MATLAB
105a107,111
> #else /* MATLAB */
>   mexPrintf ("Usage is:\n");
>   for (cptr = data; *cptr != NULL; cptr ++)
>     mexPrintf ("  %s\n", *cptr);
> #endif /* MATLAB */


[jschierm@astrid libscotch]$ diff common_old.h common.h

- fix undeclared type intptr_t
68a69
> #include            <stdint.h>                    /* added for intptr_t */

- for running as a matlab mex function
- add macro for fprintf to capture missing output
- fix exits to exit within Matlab
90a92,105
> #ifdef MATLAB
>     #include "mat.h"
>     #include "mex.h"
>     #include "matrix.h"
> 
>     #define printf mexPrintf
>     #define fprintf(file,...) (file == stdout || file == stderr ? mexPrintf(__VA_ARGS__) : fprintf(file,__VA_ARGS__))
>     #define malloc mxMalloc
>     #define calloc mxCalloc
>     #define realloc mxRealloc
>     #define free mxFree
>     #define exit(status) mexErrMsgTxt("exit=" #status)
> #endif
> 


[jschierm@astrid libscotch]$ diff dummysizes_old.c dummysizes.c

- dummysizes must run by itself during compilation and exit cleanly
269a270
> #ifndef MATLAB
270a272
> #endif /* MATLAB */


[jschierm@astrid libscotch]$ diff library_error_exit_old.c library_error_exit.c

- redirect errorPrint and errorPrintW to within Matlab
116a117
> #ifndef MATLAB
133a135,152
> #else /* MATLAB */
>   mexPrintf  ("%s", _SCOTCHerrorProgName);
> #ifdef SCOTCH_PTSCOTCH
>   if ((MPI_Initialized (&proclocnum) == MPI_SUCCESS) &&
>       (proclocnum != 0)                              &&
>       (MPI_Comm_rank (MPI_COMM_WORLD, &proclocnum) == MPI_SUCCESS))
>     mexPrintf ("(%d): ", proclocnum);
> #else /* SCOTCH_PTSCOTCH */
>   if (_SCOTCHerrorProgName[0] != '\0')
>     mexPrintf  (": ");
> #endif /* SCOTCH_PTSCOTCH */
>   mexPrintf  ("ERROR: ");
>   va_start (errlist, errstr);
>   mexPrintf (errstr, errlist);             /* Print arguments */
>   va_end   (errlist);
>   mexPrintf  ("\n");
> #endif /* MATLAB */
> 
157a177
> #ifndef MATLAB
173a194,211
> 
> #else /* MATLAB */
>   mexPrintf  ("%s", _SCOTCHerrorProgName);
> #ifdef SCOTCH_PTSCOTCH
>   if ((MPI_Initialized (&proclocnum) == MPI_SUCCESS) &&
>       (proclocnum != 0)                              &&
>       (MPI_Comm_rank (MPI_COMM_WORLD, &proclocnum) == MPI_SUCCESS))
>     mexPrintf ("(%d): ", proclocnum);
> #else /* SCOTCH_PTSCOTCH */
>   if (_SCOTCHerrorProgName[0] != '\0')
>     mexPrintf  (": ");
> #endif /* SCOTCH_PTSCOTCH */
>   mexPrintf  ("WARNING: ");
>   va_start (errlist, errstr);
>   mexPrintf (errstr, errlist);             /* Print arguments */
>   va_end   (errlist);
>   mexPrintf  ("\n");
> #endif /* MATLAB */


[jschierm@astrid scotch]$ diff Makefile_old Makefile
- add MEX rule
- add gmap_mex into scotch, clean, and install rules
- add gmap_mex object
51a52,54
> %$(MEX)       :       %.c
>                               $(CCM) $(MFLAGS) -I$(includedir) -I../libscotch -DSCOTCH_VERSION=\"$(VERSION)\" $(<) -o $(@) -L$(libdir) -l$(SCOTCHLIB) -l$(SCOTCHLIB)errexit $(LDFLAGS)
> 
70a74
>                                       gmap_mex$(MEX) \
98c102
<                                       -$(CP) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) gmap$(EXE) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) gord$(EXE) gotst$(EXE) gout$(EXE) *gtst$(EXE) gscat$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE) $(bindir)
---
>                                       -$(CP) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) gmap$(EXE) gmap_mex$(MEX) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) gord$(EXE) gotst$(EXE) gout$(EXE) *gtst$(EXE) gscat$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE) $(bindir)
108c112
<                                       -$(RM) *~ *$(OBJ) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) *gmap$(EXE) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) *gord$(EXE) gotst$(EXE) gout$(EXE) *gpart$(EXE) *gscat$(EXE) *gtst$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE)
---
>                                       -$(RM) *~ *$(OBJ) acpl$(EXE) amk_ccc$(EXE) amk_fft2$(EXE) amk_grf$(EXE) amk_hy$(EXE) amk_m2$(EXE) amk_p2$(EXE) atst$(EXE) gbase$(EXE) gcv$(EXE) *gmap$(EXE) gmap_mex$(MEX) gmk_hy$(EXE) gmk_m2$(EXE) gmk_m3$(EXE) gmk_msh$(EXE) gmk_ub2$(EXE) gmtst$(EXE) *gord$(EXE) gotst$(EXE) gout$(EXE) *gpart$(EXE) *gscat$(EXE) *gtst$(EXE) mcv$(EXE) mmk_m2$(EXE) mmk_m3$(EXE) mord$(EXE) mtst$(EXE)
239a244,251
> gmap_mex$(MEX)                :       gmap_mex.c \
>                                       ../libscotch/module.h \
>                                       ../libscotch/common.h \
>                                       $(includedir)/scotch.h \
>                                       $(libdir)/libscotch$(LIB) \
>                                       $(libdir)/libscotcherrexit$(LIB) \
>                                       gmap.h
> 


[jschierm@astrid scotch]$ diff gmap.c gmap_mex.c
- convert gmap to gmap_mex mex function with variable argument list
117,120c117,120
< int
< main (
< int                         argc,
< char *                      argv[])
---
> void mexFunction( int nlhs,
>                   mxArray *plhs[],
>                   int nrhs,
>                   const mxArray *prhs[] )
121a122,123
>   int                         argcm;
>   char                        argvm[21][257];
130a133,152
> /*  check static variables from previous runs  */
> 
>   if (C_paraNum > 0 || C_fileNum > 0)
>     mexErrMsgTxt("gmap_mex still in memory -- clear gmap_mex and try again.\n");
> 
> /*  load matlab argument list  */
> 
>   argcm=nrhs+1;
>   mexPrintf("argcm=%d\n",argcm);
>   strcpy(argvm[0],"gmap");
>   for (i=0; i<nrhs; i++)
>     if (!mxIsChar(prhs[i])) {
>       mexPrintf("%s -- prhs[%d] must be character.\n","gmap",i);
>       mexErrMsgTxt(" ");
>     }
>     else
>       mxGetString(prhs[i],argvm[i+1],256);
>   for (i=0; i<nrhs+1; i++)
>     mexPrintf("argvm[%d]=%s\n",i,argvm[i]);
> 
132,133c154,155
<   i = strlen (argv[0]);
<   if ((i >= 5) && (strncmp (argv[0] + i - 5, "gpart", 5) == 0)) {
---
>   i = strlen (argvm[0]);
>   if ((i >= 5) && (strncmp (argvm[0] + i - 5, "gpart", 5) == 0)) {
144c166
<   if ((argc >= 2) && (argv[1][0] == '?')) {       /* If need for help */
---
>   if ((argcm >= 2) && (argvm[1][0] == '?')) {       /* If need for help */
151a174,175
>   printf("point 0: C_FILENBR=%d, C_fileNbr=%d, C_paraNbr=%d\n",
>          C_FILENBR,C_fileNbr,C_paraNbr);
154,155c178,183
<   for (i = 1; i < argc; i ++) {                   /* Loop for all option codes                        */
<     if ((argv[i][0] != '-') || (argv[i][1] == '\0') || (argv[i][1] == '.')) { /* If found a file name */
---
>   for (i = 1; i < argcm; i ++) {                   /* Loop for all option codes                        */
>   printf("point 1: i=%d; C_fileNbr=%d, C_fileNum=%d, C_paraNbr=%d, C_paraNum=%d\n",
>          i,C_fileNbr,C_fileNum,C_paraNbr,C_paraNum);
>     if ((argvm[i][0] != '-') || (argvm[i][1] == '\0') || (argvm[i][1] == '.')) { /* If found a file name */
>   printf("point 2: i=%d; C_fileNbr=%d, C_fileNum=%d, C_paraNbr=%d, C_paraNum=%d\n",
>          i,C_fileNbr,C_fileNum,C_paraNbr,C_paraNum);
157,158c185,186
<         if ((C_partNbr = atoi (argv[i])) < 1)     /* Get the number of parts */
<           errorPrint ("main: invalid number of parts (\"%s\")", argv[i]);
---
>         if ((C_partNbr = atoi (argvm[i])) < 1)     /* Get the number of parts                          */
>           errorPrint ("main: invalid number of parts (\"%s\")", argvm[i]);
161a190,191
>   printf("point 3: i=%d; C_fileNbr=%d, C_fileNum=%d, C_paraNbr=%d, C_paraNum=%d\n",
>          i,C_fileNbr,C_fileNum,C_paraNbr,C_paraNum);
163c193
<         C_fileTab[C_fileNum ++].name = argv[i];
---
>         C_fileTab[C_fileNum ++].name = argvm[i];
165a196,197
>   printf("point 4: i=%d; C_fileNbr=%d, C_fileNum=%d, C_paraNbr=%d, C_paraNum=%d\n",
>          i,C_fileNbr,C_fileNum,C_paraNbr,C_paraNum);
168c200
<       switch (argv[i][1]) {
---
>       switch (argvm[i][1]) {
177c209
<           SCOTCH_stratGraphMap (&stradat, &argv[i][2]);
---
>           SCOTCH_stratGraphMap (&stradat, &argvm[i][2]);
181,182c213,214
<           for (j = 2; argv[i][j] != '\0'; j ++) {
<             switch (argv[i][j]) {
---
>           for (j = 2; argvm[i][j] != '\0'; j ++) {
>             switch (argvm[i][j]) {
192c224
<                 errorPrint ("main: invalid source graph option (\"%c\")", argv[i][j]);
---
>                 errorPrint ("main: invalid source graph option (\"%c\")", argvm[i][j]);
202,203c234,235
<           for (j = 2; argv[i][j] != '\0'; j ++) {
<             switch (argv[i][j]) {
---
>           for (j = 2; argvm[i][j] != '\0'; j ++) {
>             switch (argvm[i][j]) {
217c249
<                 errorPrint ("main: unprocessed parameter \"%c\" in \"%s\"", argv[i][j], argv[i]);
---
>                 errorPrint ("main: unprocessed parameter \"%c\" in \"%s\"", argvm[i][j], argvm[i]);
222c254
<           errorPrint ("main: unprocessed option (\"%s\")", argv[i]);
---
>           errorPrint ("main: unprocessed option (\"%s\")", argvm[i]);
225a258
>   printf("point 5\n");
230a264
>   printf("point 6\n");
232a267
>   printf("point 7\n");
235a271
>   printf("point 8\n");
239a276
>   printf("point 9\n");
244a282
>   printf("point 10\n");
251a290
>   printf("point 11\n");
256a296
>   printf("point 12\n");
264a305
>   printf("point 13\n");
271a313
>   printf("point 14\n");
273a316
>   printf("point 15\n");
276a320
>   printf("point 16\n");
277a322
>   printf("point 16a\n");
278a324
>   printf("point 16b\n");
279a326
>   printf("point 16c\n");
280a328
>   printf("point 17\n");
284a333
>   printf("point 18\n");

