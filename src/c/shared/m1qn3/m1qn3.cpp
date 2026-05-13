/*!\file: m1qn3.cpp
 * \brief: C++ implementation of the m1qn3 L-BFGS optimizer.
 *
 * Translated from Fortran m1qn3 (Version 3.3, October 2009) by
 * Jean Charles Gilbert and Claude Lemarechal, INRIA.
 * Original copyright 2008-2009 INRIA, GPL v3+.
 *
 * This translation covers only the subset used by ISSM:
 *   - Direct communication mode (reverse = 0)
 *   - DIS (Diagonal Initial Scaling), cold start
 *   - In-memory (y,s) pair storage (inmemo = true)
 *   - "dfn" norm with Euclidean scalar product
 *
 * Naming follows the Fortran original closely for auditability.
 * Fortran 1-based indices -> C++ 0-based throughout.
 * 2-D Fortran arrays ybar(n,m) / sbar(n,m) are stored as flat row-major
 * arrays ybar[j*n + i]  (j = column/update index, i = component).
 */

#include <cmath>
#include <cstring>
#include <cfloat>

#include "../io/Print/Print.h"
#include "./m1qn3.h"

/* =========================================================================
 * Internal helpers
 * ========================================================================= */

/* Euclidean dot product (replaces Fortran euclid / prosca). */
static double dot(long n, const double* x, const double* y){
   double ps = 0.0;
   for(long i = 0; i < n; ++i) ps += x[i] * y[i];
   return ps;
}

/* ctonb / ctcab: identity transformation (Euclidean metric). */
static void ctonb(long n, const double* u, double* v){
   for(long i = 0; i < n; ++i) v[i] = u[i];
}
static void ctcab(long n, const double* u, double* v){
   for(long i = 0; i < n; ++i) v[i] = u[i];
}

/* =========================================================================
 * ecube  (Fortran: ecube)
 *
 * Cubic interpolation: given f and f' at t and ta, compute a new t
 * safeguarded to [tlower, tupper].
 * ========================================================================= */
static void ecube(double& t, double f, double fp, double ta, double fa, double fpa, double tlower, double tupper){

   double z1     = fp + fpa - 3.0 * (fa - f) / (ta - t);
   double b      = z1 + fp;
   double discri;

   if (fabs(z1) <= 1.0) {
      discri = z1 * z1 - fp * fpa;
   }
   else {
      discri = fp / z1;
      discri = discri * fpa;
      discri = z1 - discri;

      if (z1 >= 0.0 && discri >= 0.0) {
         discri = sqrt(z1) * sqrt(discri);
      }
      else if (z1 <= 0.0 && discri <= 0.0) {
         discri = sqrt(-z1) * sqrt(-discri);
      }
      else {
         discri = -1.0;
      }
   }

   if (discri < 0.0) {
      t = (fp < 0.0) ? tupper : tlower;
   }
   else {
      discri = sqrt(discri);
      if (t - ta < 0.0) discri = -discri;

      double sign = (t - ta) / fabs(t - ta);
      if (b * sign > 0.0) {
         t = t + fp * (ta - t) / (b + discri);
      }
      else {
         double den  = z1 + b + fpa;
         double anum = b - discri;
         if (fabs((t - ta) * anum) < (tupper - tlower) * fabs(den)) {
            t = t + anum * (ta - t) / den;
         }
         else {
            t = tupper;
         }
      }
   }

   t = fmax(t, tlower);
   t = fmin(t, tupper);
}

/* =========================================================================
 * mupdts  (Fortran: mupdts)
 *
 * Compute m (number of L-BFGS updates) for DIS mode with in-memory storage.
 * DIS requires ndz >= 4*n + m*(2*n+1), so m = (ndz - 4*n) / (2*n+1).
 * ========================================================================= */
static long mupdts(long n, long ndz){
   return (ndz - 4 * n) / (2 * n + 1);
}

/* =========================================================================
 * dd  (Fortran: dd)
 *
 * L-BFGS two-loop recursion: compute H * depl in-place.
 * DIS mode (sscale = false).
 *
 * ybar, sbar: flat arrays of shape [m][n] (update j at offset j*n).
 * alpha:      scratch array of length m.
 * jmin, jmax: 0-based circular buffer pointers (inclusive range in the ring).
 * precos:     unused here (DIS mode uses diag scaling).
 * ========================================================================= */
static void dd(long n, long m, double* depl, double* aux, long jmin, long jmax, double /*precos*/, double* diag, double* alpha, double* ybar, double* sbar, long izs[], float rzs[], void* dzs){
   /* jfin = jmax, but if the ring has wrapped, jfin = jmax + m */
   long jfin = jmax;
   if (jfin < jmin) jfin = jmax + m;

   /* --- descent phase --- */
   for (long j = jfin; j >= jmin; --j) {
      long jp = j % m;   /* 0-based circular index */
      double ps = dot(n, depl, &sbar[jp * n]);
      alpha[jp] = ps;
      for (long i = 0; i < n; ++i)
         depl[i] -= ps * ybar[jp * n + i];
   }

   /* --- DIS preconditioning: depl = D * ctonb(depl), then ctcab --- */
   ctonb(n, depl, aux);
   for (long i = 0; i < n; ++i) aux[i] *= diag[i];
   ctcab(n, aux, depl);

   /* --- ascent phase --- */
   for (long j = jmin; j <= jfin; ++j) {
      long jp = j % m;
      double ps = dot(n, depl, &ybar[jp * n]);
      double r  = alpha[jp] - ps;
      for (long i = 0; i < n; ++i)
         depl[i] += r * sbar[jp * n + i];
   }
}

/* =========================================================================
 * mlis3  (Fortran: mlis3)
 *
 * Wolfe-condition line search with cubic interpolation.
 * Direct communication mode only.
 *
 * On exit: logic =
 *   0  serious step (both Wolfe conditions satisfied)
 *   1  blocked on tmax
 *   4  napmax exceeded
 *   5  user requested stop (indic = 0)
 *   6  stopped on dxmin / inconsistent
 *  <0  simulator failure
 * ========================================================================= */
static void mlis3(long n, M1qn3SimulFunc simul, double* x, double& f, double& fpn, double& t, double tmin, double tmax, double* d, double* g, double amd, double amf, long impres, int& logic, long& nap, long napmax, double* xn, long izs[], float rzs[], void* dzs){
   /* --- sanity check --- */
   if (!(n > 0 && fpn < 0.0 && t > 0.0 && tmax > 0.0
         && amf > 0.0 && amd > amf && amd < 1.0)) {
      logic = 6;
      return;
   }

   double tesf = amf * fpn;
   double tesd = amd * fpn;

   const double barmin  = 0.01;
   const double barmul  = 5.0;
   const double barmax  = 0.3;
   double       barr    = barmin;

   double td = 0.0, tg = 0.0;
   double fn = f, fg = fn, fpg = fpn;
   double ta = 0.0, fa = fn, fpa = fpn;
   double d2 = dot(n, d, d);

   /* eliminate ridiculously small initial t */
   if (t < tmin) {
      t = tmin;
      if (t > tmax) tmin = tmax;
   }
   bool t_increased = false;
   while (fn + t * fpn >= fn + 0.9 * t * fpn) {
      t_increased = true;
      t *= 2.0;
   }

   int  indica = 1;
   logic = 0;
   if (t > tmax) { t = tmax; logic = 1; }

   /* trial point */
   for (long i = 0; i < n; ++i) { xn[i] = x[i]; x[i] = xn[i] + t * d[i]; }

   /* ---- main loop ---- */
   bool done        = false;
   bool napmax_exit = false;   /* true when we exit via napmax (no simul called) */
   while (!done) {
      nap++;
      if (nap > napmax) {
         /* napmax exceeded: restore fn/xn to left bracket, but leave x at the
          * current trial position (matching Fortran: x is not restored here,
          * only xn/fn are updated; the final x = xn assignment below is
          * skipped via napmax_exit so that x stays at the last trial point
          * where g was evaluated). */
         logic = 4;
         fn = fg;
         for (long i = 0; i < n; ++i) xn[i] += tg * d[i];
         napmax_exit = true;
         break;
      }

      /* call simulator: ask for f and g */
      bool need_advance = false;
      {
         long indic = 4;
         simul(&indic, &n, x, &f, g, izs, rzs, dzs);

         if (indic == 0) {
            /* user stop */
            logic = 5;
            fn = f;
            for (long i = 0; i < n; ++i) xn[i] = x[i];
            done = true;
         }
         else if (indic < 0) {
            /* computation failed: shrink interval and re-enter advance.
             * Fortran goes to label 905 (skipping fa=f;fpa=fp update). */
            td    = t;
            logic = 0;
            t     = tg + 0.1 * (td - tg);
            /* indica = indic (goes to 905, NOT 900 — do not update fa/fpa) */
            indica = (int)indic;
            need_advance = true;
         }
      }

      if (!done && !need_advance) {
         double fp  = dot(n, d, g);
         double ffn = f - fn;

         if (ffn > t * tesf) {
            /* first Wolfe condition violated → bracket from above, interpolate.
             * Fortran: td=t; go to 500; then fa=f; fpa=fp at label 900.
             * ecube must be called with the OLD fa/fpa (previous anchor),
             * so do NOT overwrite fa/fpa before calling ecube. */
            td     = t;
            logic  = 0;

            /* interpolation (label 500 in Fortran) */
            if (indica <= 0) {
               ta = t;
               t  = 0.9 * tg + 0.1 * td;
            }
            else {
               double test   = barr * (td - tg);
               double gauche = tg + test;
               double droite = td - test;
               double taa    = t;
               ecube(t, f, fp, ta, fa, fpa, gauche, droite);
               ta = taa;
               if (t > gauche && t < droite)
                  barr = fmax(barmin, barr / barmul);
               else
                  barr = fmin(barmul * barr, barmax);
            }
            /* label 900: update anchor AFTER ecube */
            fa = f; fpa = fp;
            need_advance = true;
         }
         else if (fp > tesd) {
            /* first Wolfe ok, curvature condition satisfied: serious step → accept.
             * Fortran: fp > tesd → logic=0; go to 320 (accept). */
            logic = 0;
            fn = f;
            for (long i = 0; i < n; ++i) xn[i] = x[i];
            done = true;
         }
         else {
            /* first Wolfe ok, curvature condition NOT satisfied.
             * Fortran label 350: tg=t; fg=f; fpg=fp; then extrapolate (td==0)
             * or interpolate (td!=0) — regardless of whether logic==0 or logic==1.
             * The old "else if (logic==0) accept" branch was WRONG: it caused the
             * line search to return a step that doesn't satisfy the strong Wolfe
             * curvature condition, corrupting (y,s) pairs and producing <y,s><=0. */
            tg = t; fg = f; fpg = fp;

            if (td != 0.0) {
               /* interpolation (label 500 in Fortran) */
               if (indica <= 0) {
                  ta = t;
                  t  = 0.9 * tg + 0.1 * td;
               }
               else {
                  double test   = barr * (td - tg);
                  double gauche = tg + test;
                  double droite = td - test;
                  double taa    = t;
                  ecube(t, f, fp, ta, fa, fpa, gauche, droite);
                  ta = taa;
                  if (t > gauche && t < droite)
                     barr = fmax(barmin, barr / barmul);
                  else
                     barr = fmin(barmul * barr, barmax);
               }
               /* label 900: update anchor AFTER ecube */
               fa = f; fpa = fp;
               need_advance = true;
            }
            else {
               /* extrapolation */
               double taa    = t;
               double gauche = (1.0 + barmin) * t;
               double droite = 10.0 * t;
               ecube(t, f, fp, ta, fa, fpa, gauche, droite);
               ta = taa;
               /* label 900: update anchor AFTER ecube */
               fa = f; fpa = fp;
               if (t >= tmax) {
                  logic = 1;
                  t = tmax;
               }
               need_advance = true;
            }
         }
      }

      if (!done && need_advance) {
         indica = 1; /* treated as "computed" for next cubic */

         /* check stopping: td - tg < tmin? */
         if (td != 0.0) {
            bool stop_dxmin = false;
            if (td - tg < tmin) {
               stop_dxmin = true;
            }
            else {
               /* machine-precision guard */
               bool no_change = true;
               for (long i = 0; i < n; ++i) {
                  double z = xn[i] + t * d[i];
                  if (z != xn[i] && z != x[i]) { no_change = false; break; }
               }
               if (no_change) stop_dxmin = true;
            }

            if (stop_dxmin) {
               logic = 6;
               if (tg != 0.0) {
                  fn = fg;
                  for (long i = 0; i < n; ++i) xn[i] += tg * d[i];
               }
               done = true;
            }
         }

         if (!done) {
            /* update trial point and loop */
            for (long i = 0; i < n; ++i) x[i] = xn[i] + t * d[i];
         }
      }
   }

   /* restore x = best point found.
    * For normal exits (accepted step, dxmin, user-stop): x <- xn (best safe point).
    * For napmax exit: leave x as the last trial point where simul was called,
    * matching Fortran m1qn3 behaviour (x is not explicitly restored there). */
   f = fn;
   if (!napmax_exit) {
      for (long i = 0; i < n; ++i) x[i] = xn[i];
   }
}

/* =========================================================================
 * m1qn3a  (Fortran: m1qn3a)
 *
 * Core L-BFGS minimizer.  DIS mode, cold start, direct comms, in-memory.
 * ========================================================================= */
static void m1qn3a(M1qn3SimulFunc simul,
                   long n, double* x, double& f, double* g,
                   double dxmin, double df1, double& epsg,
                   long impres,
                   long& omode, long& niter, long& nsim,
                   long m,
                   long& jmin, long& jmax,
                   double* d,    /* length n */
                   double* gg,   /* length n */
                   double* diag, /* length n */
                   double* aux,  /* length n */
                   double* alpha,/* length m */
                   double* ybar, /* length m*n */
                   double* sbar, /* length m*n */
                   long izs[], float rzs[], void* dzs,
                   bool return_best = true)
{
   const double rm1  = 0.0001;   /* Armijo constant */
   const double rm2  = 0.99;     /* curvature constant */
   const double rmin = 1.0e-20;

   /* --- initialisation --- */
   long itmax = niter;
   niter = 0;
   long isim  = 1;      /* already called once before entering here */
   double eps1 = 1.0;

   /* initial gradient norm (dfn norm = Euclidean) */
   double gnorm  = sqrt(dot(n, g, g));
   double gnorms = gnorm;  /* dfn-norm = Euclidean norm */

   if (impres >= 1) {
      _printf0_("     f             = " << f << "\n");
      _printf0_("     dfn-norm of g = " << gnorms << "\n");
   }
   if (gnorms < rmin) {
      omode = 2;
      if (impres >= 1) _printf0_("   >>> m1qn3a: initial gradient is too small\n");
      return;
   }

   /* best-point tracking: keep the (f,x,g) triplet with the lowest f seen */
   double  f_best = f;
   double* x_best = new double[n];
   double* g_best = new double[n];
   for(long i = 0; i < n; ++i){ x_best[i] = x[i]; g_best[i] = g[i]; }

   /* cold start: jmin=0, jmax=-1 (empty ring, 0-based) */
   jmin = 0;
   jmax = -1;
   long jcour = jmax;   /* current slot (unused for cold start) */

   /* Fletcher scaling of first descent direction, initialise diag = 1 */
   double precos = 2.0 * df1 / (gnorm * gnorm);
   for (long i = 0; i < n; ++i) {
      d[i]    = -g[i] * precos;
      diag[i] = 1.0;
   }

   /* check descent */
   double tmax = 1.0e20;
   double hp0  = dot(n, d, g);
   if (hp0 >= 0.0) {
      omode = 7;
      if (impres >= 1)
         _printf0_("   >>> m1qn3 (iter 0): search direction is not descent: (g,d)=" << hp0 << "\n");
      goto m1qn3a_exit;
   }

   /* ---- main iteration loop ---- */
   for (;;) {
      niter++;

      if (impres >= 4)
         _printf0_("   m1qn3: iter " << niter << ", simul " << isim
                   << ", f=" << f << ", h'(0)=" << hp0 << "\n");

      /* save current gradient */
      for (long i = 0; i < n; ++i) gg[i] = g[i];
      double ff = f;

      /* compute tmin */
      double tmin_ls = 0.0;
      for (long i = 0; i < n; ++i)
         tmin_ls = fmax(tmin_ls, fabs(d[i]));
      tmin_ls = dxmin / tmin_ls;

      double t   = 1.0;
      double d1  = hp0;
      int    moderl;

      /* --- line search --- */
      mlis3(n, simul, x, f, d1, t, tmin_ls, tmax, d, g,
            rm2, rm1, impres, moderl, isim, nsim, aux,
            izs, rzs, dzs);

      /* handle line-search exit codes */
      if (moderl != 0) {
         if (moderl < 0) {
            omode = moderl; goto m1qn3a_exit;
         }
         else if (moderl == 4) {
            omode = 5; goto m1qn3a_exit;
         }
         else if (moderl == 5) {
            omode = 0; goto m1qn3a_exit;
         }
         else if (moderl == 6) {
            omode = 6; goto m1qn3a_exit;
         }
         /* moderl == 1 (blocked on tmax): skip L-BFGS update, fall through */
      }

      /* --- L-BFGS matrix update (skipped when blocked on tmax) --- */
      if (moderl != 1 && m > 0) {
         /* advance ring pointer */
         jmax = jmax + 1;
         if (jmax >= m) jmax -= m;
         /* advance jmin only when ring is full (cold start: Fortran niter > m) */
         if (niter > m) {
            jmin++;
            if (jmin >= m) jmin -= m;
         }
         jcour = jmax;

         /* y = g_new - g_old, s = t * d */
         for (long i = 0; i < n; ++i) {
            sbar[jcour * n + i] = t * d[i];
            ybar[jcour * n + i] = g[i] - gg[i];
         }

         /* ys = (y, s) — must be positive */
         double ys = dot(n, &ybar[jcour * n], &sbar[jcour * n]);
         if (ys <= 0.0) {
            omode = 7;
            if (impres >= 1)
               _printf0_("   >>> m1qn3 (iter " << niter
                         << "): (y,s) = " << ys << " is not positive\n");
            goto m1qn3a_exit;
         }

         /* normalise ybar and sbar by 1/sqrt(ys) */
         double inv_sqrt_ys = 1.0 / sqrt(ys);
         for (long i = 0; i < n; ++i) {
            sbar[jcour * n + i] *= inv_sqrt_ys;
            ybar[jcour * n + i] *= inv_sqrt_ys;
         }

         /* DIS: update diagonal preconditioner.
          * Scale to Rayleigh ellipsoid, then rank-1 update.
          * aux = ctonb(ybar_j), then scale by 1/(ybar_j^T D ybar_j)
          * then update: D <- D - (D*sbar_j)(D*sbar_j)^T/(sbar_j^T D^{-1} sbar_j)
          *                      + ybar_j ybar_j^T  (in orthonormal basis)
          * Following the Fortran exactly:
          */
         {
            /* ctonb(ybar_j) -> aux */
            ctonb(n, &ybar[jcour * n], aux);
            double ps = 0.0;
            for (long i = 0; i < n; ++i) ps += diag[i] * aux[i] * aux[i];
            double d1_scale = 1.0 / ps;
            for (long i = 0; i < n; ++i) diag[i] *= d1_scale;

            /* ctonb(sbar_j) -> gg (used as temp) */
            ctonb(n, &sbar[jcour * n], gg);
            ps = 0.0;
            for (long i = 0; i < n; ++i) ps += gg[i] * gg[i] / diag[i];
            double den = ps;
            for (long i = 0; i < n; ++i) {
               diag[i] = 1.0 / (1.0 / diag[i] + aux[i] * aux[i]
                                 - (gg[i] / diag[i]) * (gg[i] / diag[i]) / den);
               if (diag[i] <= 0.0) diag[i] = rmin;
            }
            /* restore gg to current gradient (was overwritten by ctonb) */
            for (long i = 0; i < n; ++i) gg[i] = g[i];
         }
      }

      /* --- update best point seen so far (after line search + L-BFGS update) --- */
      if (return_best && f < f_best) {
         f_best = f;
         for (long i = 0; i < n; ++i) { x_best[i] = x[i]; g_best[i] = g[i]; }
      }

      /* --- stopping tests --- */
      gnorm = sqrt(dot(n, g, g));
      eps1  = gnorm / gnorms;

      if (impres >= 3)
         _printf0_("   m1qn3: iter " << niter << ", simul " << isim
                   << ", step=" << t << ", f=" << f
                   << ", |g|=" << gnorm << ", |g|/|g0|=" << eps1 << "\n");

      if (eps1 < epsg) {
         omode = 1;
         epsg  = eps1;
         nsim  = isim;
         goto m1qn3a_exit;
      }
      if (niter == itmax) {
         omode = 4;
         if (impres >= 1)
            _printf0_("   >>> m1qn3 (iter " << niter << "): max iterations reached\n");
         epsg = eps1; nsim = isim;
         goto m1qn3a_exit;
      }
      if (isim > nsim) {
         omode = 5;
         if (impres >= 1)
            _printf0_("   >>> m1qn3 (iter " << niter << "): max simulations reached\n");
         epsg = eps1; nsim = isim;
         goto m1qn3a_exit;
      }

      /* --- compute new descent direction d = -H * g --- */
      for (long i = 0; i < n; ++i) d[i] = -g[i];
      if (m > 0 && jmax >= 0) {
         dd(n, m, d, aux, jmin, jmax, precos, diag, alpha, ybar, sbar,
            izs, rzs, dzs);
      }
      else {
         /* m=0 fallback: steepest descent scaled by preco */
         double ps   = dot(n, g, g);
         double prec = 2.0 * (ff - f) / ps;
         for (long i = 0; i < n; ++i) d[i] = -g[i] * prec;
      }

      /* check descent */
      hp0 = dot(n, d, g);
      if (hp0 >= 0.0) {
         omode = 7;
         if (impres >= 1)
            _printf0_("   >>> m1qn3 (iter " << niter
                      << "): d is not a descent direction: (g,d)=" << hp0 << "\n");
         epsg = eps1; nsim = isim;
         goto m1qn3a_exit;
      }
   } /* end main loop */

m1qn3a_exit:
   /* Restore the best (f, x, g) triplet seen during the run */
   if (return_best && f_best < f) {
      f = f_best;
      for (long i = 0; i < n; ++i){
			x[i] = x_best[i];
			g[i] = g_best[i];
		}
   }
   delete [] x_best;
   delete [] g_best;
}

/* =========================================================================
 * m1qn3_cpp  (public entry point, Fortran: m1qn3)
 * ========================================================================= */
void m1qn3_cpp(M1qn3SimulFunc simul,
               long* n_ptr, double* x, double* f_ptr, double* g,
               double* dxmin, double* df1, double* epsg,
               long* impres, long* /*io*/,
               long* omode, long* niter, long* nsim,
               long* iz, double* dz, long* ndz,
               long izs[], float rzs[], void* dzs,
               bool return_best)
{
   long   n      = *n_ptr;
   double f      = *f_ptr;
   long   impres_ = *impres;

   /* --- input checks --- */
   if (n <= 0)          { *omode = 2; return; }
   if (*niter <= 0)     { *omode = 2; return; }
   if (*nsim  <= 0)     { *omode = 2; return; }
   if (*dxmin <= 0.0)   { *omode = 2; return; }
   if (*epsg  <= 0.0)   { *omode = 2; return; }
   if (*epsg  >= 1.0)   { *omode = 1; *niter = 0; *nsim = 0; return; }

   /* --- compute m (number of L-BFGS updates) --- */
   long m = mupdts(n, *ndz);
   if (m < 1) {
      *omode = 2;
      if (impres_ >= 1)
         _printf0_("   >>> m1qn3_cpp: not enough memory (ndz too small)\n");
      return;
   }
   long ndz_needed = 4 * n + m * (2 * n + 1);
   if (*ndz < ndz_needed) { *omode = 2; return; }

   if (impres_ >= 1) {
      _printf0_("   m1qn3_cpp (C++ implementation): entry point\n");
      _printf0_("     n=" << n << ", dxmin=" << *dxmin << ", df1=" << *df1
                << ", epsg=" << *epsg << "\n");
      _printf0_("     max iter=" << *niter << ", max sim=" << *nsim << "\n");
      _printf0_("     allocated ndz=" << *ndz << ", used=" << ndz_needed
                << ", updates m=" << m << "\n");
   }

   /* --- split working array dz (DIS mode layout):
    *
    *   dz[0 .. n-1]            : diag  (diagonal preconditioner)
    *   dz[n .. n+n*m-1]        : ybar  (m blocks of n: ybar[j*n+i])
    *   dz[n+n*m .. n+2*n*m-1]  : sbar
    *   dz[n+2*n*m .. 2*n+2*n*m-1] : d  (descent direction)
    *   dz[2*n+2*n*m .. 3*n+2*n*m-1] : gg (old gradient)
    *   dz[3*n+2*n*m .. 4*n+2*n*m-1] : aux
    *   dz[4*n+2*n*m .. 4*n+2*n*m+m-1]: alpha
    *
    * Total = 4*n + m*(2*n+1) = ndz_needed.
    --- */
   double* diag  = dz;
   double* ybar  = diag + n;
   double* sbar  = ybar + n * m;
   double* d_vec = sbar + n * m;
   double* gg    = d_vec + n;
   double* aux   = gg + n;
   double* alpha = aux + n;

   /* store problem size in iz (for potential warm-restart info) */
   iz[0] = n;
   iz[1] = 0;   /* DIS */
   iz[2] = m;

   long jmin = 0, jmax = -1;

   /* Call the core optimizer */
   m1qn3a(simul, n, x, f, g,
          *dxmin, *df1, *epsg,
          impres_,
          *omode, *niter, *nsim,
          m, jmin, jmax,
          d_vec, gg, diag, aux, alpha, ybar, sbar,
          izs, rzs, dzs,
          return_best);

   /* store final pointers */
   iz[3] = jmin;
   iz[4] = jmax;

   *f_ptr = f;

   if(impres_ >= 1){
      _printf0_("   m1qn3_cpp: exit code " << *omode << ", iter=" << *niter << ", sim=" << *nsim << ", |g|/|g0|=" << *epsg << "\n");
   }
}
