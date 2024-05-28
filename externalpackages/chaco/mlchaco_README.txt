mlchaco (from meshpartdist.tar, 2/08/02)
http://www.cerfacs.fr/algor/Softs/MESHPART/


[jschierm@astrid chaco]$ diff mlchaco_old.c mlchaco.c
34a35  [CLK_TCK undefined, so use CLOCKS_PER_SEC]
> #include <time.h>      /*  CLOCKS_PER_SEC  */
56,61c57,62  [update for current Matlab]
< void mexFunction(    
<     int         nlhs,           /* number of expected outputs */
<     Matrix      *plhs[],        /* matrix pointer array returning outputs */
<     int         nrhs,           /* number of inputs */
<     Matrix      *prhs[]         /* matrix pointer array for inputs */
<     )
---
> void mexFunction(
>     int           nlhs,           /* number of expected outputs */
>     mxArray       *plhs[],        /* array of pointers to output arguments
>     */
>     int           nrhs,           /* number of inputs */
>     const mxArray *prhs[]         /* array of pointers to input arguments */
> )
87a89  [update for current Matlab]
>     mwIndex *mwstart,*mwadjacency;
104c106,110  [update for current Matlab]
<     start = mxGetJc(A_in);
---
> /*    start = mxGetJc(A_in);*/
>     mwstart = mxGetJc(A_in);
>     start = mxMalloc((mxGetN(A_in)+1)*sizeof(int));
>     for (i=0; i<(mxGetN(A_in)+1); i++)
>         start[i]= (int)mwstart[i];
106c112,116  [update for current Matlab]
<     adjacency = mxGetIr(A_in);
---
> /*    adjacency = mxGetIr(A_in);*/
>     mwadjacency = mxGetIr(A_in);
>     adjacency = mxMalloc(mxGetNzmax(A_in)*sizeof(int));
>     for (i=0; i<mxGetNzmax(A_in); i++)
>         adjacency[i]= (int)mwadjacency[i];
131,132c141,144  [provide default filenames, since no way to input]
<     outassignname = NULL;
<     outfilename = NULL;
---
> /*    outassignname = NULL;
>     outfilename = NULL;*/
>     outassignname = "chaco_assign.txt";
>     outfilename = "chaco_out.txt";
169c181  [update for current Matlab]
<       plhs [1] = mxCreateFull (1, 1, REAL) ;
---
>       plhs [1] = mxCreateDoubleMatrix (1, 1, mxREAL) ;
173c185  [CLK_TCK undefined, so use CLOCKS_PER_SEC]
<       ((double) CLK_TCK) ;
---
>       ((double) CLOCKS_PER_SEC) ;
179c191  [update for current Matlab]
<         map_out = mxCreateFull(1,nvtxs,REAL);
---
>         map_out = mxCreateDoubleMatrix(1,nvtxs,mxREAL);
185a198,199  [update for current Matlab]
>     if (start != NULL) mxFree((char *) start);
>     if (adjacency != NULL) mxFree((char *) adjacency);


[jschierm@astrid chaco]$ diff Makefile_old Makefile
27c27,28  [current Matlab location]
< MATLAB =        /usr/local/libexec/matlab
---
> #MATLAB =        /usr/local/libexec/matlab
> MATLAB =        /usr/local/pkgs/matlab-7.6
30c31,32  [current Chaco 2.2 location]
< CHACO =         ../../Chaco-2.0/code
---
> #CHACO =         ../../Chaco-2.0/code
> CHACO =         ../../Chaco-2.2/code
34,35c36,37  [use gcc instead of cc]
< #CC =            gcc 
< CC =            cc 
---
> CC =            gcc 
> #CC =            cc 
37c39,40  [add MATLAB flag to compile for mex-function]
< CFLAGS =        -Xa -G -xO4 -xcg92
---
> #CFLAGS =        -Xa -G -xO4 -xcg92
> CFLAGS =        -Xa -G -xO4 -xcg92 -DMATLAB
39,40c42,45  [current function locations]
< AR =             /usr/ccs/bin/ar rcv   # for solaris 2
< RANLIB =         /usr/ccs/bin/ranlib   # for solaris 2
---
> #AR =             /usr/ccs/bin/ar rcv   # for solaris 2
> AR =             /usr/bin/ar rcv   # for solaris 2
> #RANLIB =         /usr/ccs/bin/ranlib   # for solaris 2
> RANLIB =         /usr/bin/ranlib   # for solaris 2
45,47d49  [use Chaco versions with MATLAB switch rather than local versions]
<               ${MLCHACO}/chaco_check_graph.c \
<               ${MLCHACO}/check_input.c \
<                 ${MLCHACO}/smalloc.c \
48a51,53
> #             ${MLCHACO}/chaco_check_graph.c \
> #             ${MLCHACO}/check_input.c \
> #                ${MLCHACO}/smalloc.c \
55a61  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/input/check_input.c \
61a68  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/graph/check_graph.c \
91a99  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/klvspiff/flow.c \
97a106  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/klvspiff/flatten.c \
100c109,110  [update CHLIST to match Chaco 2.2]
<               ${CHACO}/coarsen/makecgraph.c \
---
>               ${CHACO}/coarsen/makefgraph.c \
>               ${CHACO}/coarsen/makeccoords.c \
102d111  [update CHLIST to match Chaco 2.2]
<               ${CHACO}/coarsen/countcedges.c \
108a118  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/coarsen/maxmatch5.c \
234a245  [update CHLIST to match Chaco 2.2]
>               ${CHACO}/util/smalloc.c \
245a257,258  [use local versions rather than Chaco versions (added for clarity)]
> #             ${CHACO}/main/user_params.c \
> #             ${CHACO}/util/bail.c \
254a268,271  [update for current mex]
> #mlchaco:     ${MLFILES.c} chaco.a Makefile
> #             mex -V4 -output mlchaco ${MLFILES.c} chaco.a -I${CHACO}/main
> #             mv mlchaco.mex* ${DEST_DIR}
> 
256c273  [update for current mex]
<               mex -V4 -output mlchaco ${MLFILES.c} chaco.a -I${CHACO}/main
---
>               mex -output mlchaco -largeArrayDims -DMATLAB ${MLFILES.c}
>               chaco.a -I${CHACO}/main

