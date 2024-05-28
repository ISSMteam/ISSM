Chaco-2.2.tar (5/03/00? -- 2.0 was 2/95)
http://www.sandia.gov/~bahendr/chaco.html


Note that meshpart pulls the object files from here, compiled with
the MATLAB flag for use in the mlchaco mex function, and assembles
them into a library in its own directory.  This means that any
objects compiled for the mlchaco mex function and left here should
not be used for the stand-alone chaco executable, and vice versa.


[jschierm@astrid main]$ diff defs_old.h defs.h
11a12,27  [for running as a matlab mex function]
> 
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
> 
>     #define check_graph chaco_check_graph
> #endif
> 


added #include "defs.h" to the following:
assign/assign_out.c
eigen/get_extval.c
klvspiff/matching.c
misc/timing.c
util/bail.c
util/checkpnt.c
util/doubleout.c
util/smalloc.c
util/strout.c


[jschierm@astrid code]$ diff Makefile_old Makefile
4,5c4,9  [fPIC required, CFLAGS and OFLAGS copied from Cielo gccopts_v75.sh for glnxa64]
< IFLAG =               -Imain
< CFLAGS =      -O2
---
> IFLAG =               -Imain -I/usr/local/pkgs/matlab-7.6/extern/include
> #CFLAGS =     -O2
> #OFLAGS =     -O2
> #CFLAGS =     -fPIC -fno-omit-frame-pointer -D_GNU_SOURCE -pthread -fexceptions
> CFLAGS =      -fPIC -fno-omit-frame-pointer -pthread -fexceptions -DMATLAB
> #CFLAGS =     -fPIC -fno-omit-frame-pointer -pthread -fexceptions

