/*This is needed, because the chaco library defines a "C" function seconds that conflicts with the Metis version.: */

#if defined(_INTEL_WIN_) || defined(_MSYS2_)
#include   <time.h>
#else
#include   <sys/time.h>
#include   <sys/resource.h>
#endif

double chaco_seconds(void){

	double    curtime;

#ifdef RUSAGE_SELF

/* This timer is faster and more robust (if it exists). */
    struct rusage rusage;
    /*int getrusage(); commenting this out. not sure why it's there anymore
	 *as it clobbers the prototype int getrusag(int target,rusage* results) which 
	 *is defined in the <sys/time.h> and <sys/resource.h> header files. Leaving it 
	 *for reference in case we have a problem here in the future*/

    getrusage(RUSAGE_SELF, &rusage);
    curtime = ((rusage.ru_utime.tv_sec + rusage.ru_stime.tv_sec) +
	    1.0e-6 * (rusage.ru_utime.tv_usec + rusage.ru_stime.tv_usec));

#else

/* ANSI timer, but lower resolution & wraps around after ~36 minutes. */

    curtime = clock()/((double) CLOCKS_PER_SEC);

#endif

    return (curtime);
}
