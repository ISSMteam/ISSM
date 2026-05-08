/*!\file: m1qn3.h
 * \brief: C++ implementation of the m1qn3 L-BFGS optimizer.
 *
 * Translated from the Fortran m1qn3 (Version 3.3, October 2009) by
 * Jean Charles Gilbert and Claude Lemarechal, INRIA.
 *
 * This translation covers only the subset used by ISSM:
 *   - Direct communication mode (reverse = 0)
 *   - DIS (Diagonal Initial Scaling) mode (imode[0] = 0)
 *   - Cold start (imode[1] = 0)
 *   - In-memory (y,s) pair storage
 *   - "dfn" norm (Euclidean scalar product)
 */

#ifndef _M1QN3_CPP_H_
#define _M1QN3_CPP_H_

/*! Callback type for the cost function / simulator.
 *  indic = 4 on entry: compute f and g.
 *  indic = 0 on return: user requests stop.
 *  indic < 0 on return: computation failed at this point.
 */
typedef void (*M1qn3SimulFunc)(long* indic, long* n, double* x, double* pf,
                               double* g, long izs[], float rzs[], void* dzs);

/*! m1qn3_cpp - L-BFGS minimizer (C++ implementation).
 *
 * Parameters (matching the Fortran m1qn3_ interface, minus unused args):
 *   simul   [in]     cost-function callback (same signature as Fortran simul)
 *   n       [in]     problem dimension
 *   x       [in/out] initial guess on entry, solution on exit
 *   f       [in/out] f(x) on entry (already evaluated), final f on exit
 *   g       [in/out] grad f on entry (already evaluated), final grad on exit
 *   dxmin   [in]     minimum step-size for line search
 *   df1     [in]     expected decrease in f during the first iteration
 *   epsg    [in/out] relative gradient tolerance on entry; achieved value on exit
 *   impres  [in]     verbosity level (0 = silent)
 *   io      [in]     (unused, kept for interface compatibility)
 *   omode   [out]    exit code (1=converged, 4=max iter, 5=max sim, ...)
 *   niter   [in/out] max iterations on entry; iterations performed on exit
 *   nsim    [in/out] max simulations on entry; simulations performed on exit
 *   iz      [work]   integer working array of size 5
 *   dz      [work]   double working array of size ndz
 *   ndz     [in]     size of dz; must be >= 4*n + m*(2*n+1) for some m >= 1
 *   izs     [in]     passed through to simul (not used internally)
 *   rzs     [in]     passed through to simul (not used internally)
 *   dzs     [in]     passed through to simul (user data pointer)
 */
void m1qn3_cpp(M1qn3SimulFunc simul,
               long* n, double* x, double* f, double* g,
               double* dxmin, double* df1, double* epsg,
               long* impres, long* io,
               long* omode, long* niter, long* nsim,
               long* iz, double* dz, long* ndz,
               long izs[], float rzs[], void* dzs,
               bool return_best = true);

#endif /* _M1QN3_CPP_H_ */
