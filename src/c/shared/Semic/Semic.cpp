/*!\file semic.cpp
 * \brief C++ port of the SEMIC surface-physics module.
 *
 * Original Fortran source:
 *   externalpackages/semic/src/surface_physics.f90   (Mario Krapp)
 *   src/c/modules/SurfaceMassBalancex/run_semic.f90
 *   src/c/modules/SurfaceMassBalancex/run_semic_transient.f90
 *
 * Every subroutine, function, and literal constant is translated verbatim.
 * Names follow the Fortran originals as closely as C++ allows.
 */

#include <config.h>
#include "./Semic.h"
#include "../io/Print/Print.h"
#include "../Exceptions/exceptions.h"
#include "../Numerics/recast.h"

#include <algorithm>   /* std::min, std::max                    */
#include <cmath>       /* exp, log, sqrt, tanh, acos, fabs ...  */
#include <limits>      /* machine epsilon*/
#include <cstring>     /* memset                                 */
#include <stdexcept>
#include <string>
#include <vector>

using namespace semic_const;

/* ======================================================================== */
/*  Helpers – scalar physics (identical to Fortran ELEMENTAL routines)       */
/* ======================================================================== */

/* Saturation water-vapour pressure over water (Magnus formula) */
static inline IssmDouble ew_sat(IssmDouble t) {
	return 611.2 * std::exp(17.62 * (t - t0) / (243.12 + t - t0));
}

/* Saturation water-vapour pressure over ice */
static inline IssmDouble ei_sat(IssmDouble t) {
	return 611.2 * std::exp(22.46 * (t - t0) / (272.62 + t - t0));
}

/* ---------------------------------------------------------------------- */
IssmDouble semic_sensible_heat_flux(IssmDouble ts, IssmDouble ta, IssmDouble wind,
                                IssmDouble rhoatm, IssmDouble csh, IssmDouble cap_air)
{
	return csh * cap_air * rhoatm * wind * (ts - ta);
}

/* ---------------------------------------------------------------------- */
void semic_latent_heat_flux(IssmDouble ts, IssmDouble wind, IssmDouble shum, IssmDouble sp,
                            IssmDouble rhoatm, int mask, IssmDouble clh,
                            IssmDouble& lhf, IssmDouble& subl, IssmDouble& evap)
{
	subl = 0.0;
	evap = 0.0;
	lhf  = 0.0;
	if (ts < t0) {
		IssmDouble esat_sur = ei_sat(ts);
		IssmDouble shum_sat = esat_sur * eps / (esat_sur * (eps - 1.0) + sp);
		subl = clh * wind * rhoatm * (shum_sat - shum);
		lhf  = subl * cls;
	} else {
		IssmDouble esat_sur = ew_sat(ts);
		IssmDouble shum_sat = esat_sur * eps / (esat_sur * (eps - 1.0) + sp);
		evap = clh * wind * rhoatm * (shum_sat - shum);
		lhf  = evap * clv;
	}
}

/* ---------------------------------------------------------------------- */
IssmDouble semic_longwave_upward(IssmDouble ts) {
	return sigm * ts * ts * ts * ts;
}

/* ---------------------------------------------------------------------- */
void semic_diurnal_cycle(IssmDouble amp, IssmDouble tmean, IssmDouble& above, IssmDouble& below)
{
	IssmDouble tmp1 = 0.0, tmp2 = 0.0;
	if (std::fabs(tmean / amp) < 1.0) {
		tmp1 = std::acos(tmean / amp);
		tmp2 = std::sqrt(1.0 - tmean * tmean / (amp * amp));
	}
	if (tmean + amp < 0.0) {
		below = tmean;
		above = 0.0;
	} else {
		above = tmean;
		below = 0.0;
		if (std::fabs(tmean) < amp) {
			above = (-tmean * tmp1 + amp * tmp2 + pi * tmean) / (pi - tmp1);
			below = (tmean * tmp1 - amp * tmp2) / tmp1;
		}
	}
}

/* ---------------------------------------------------------------------- */
/*  Albedo schemes                                                          */
/* ---------------------------------------------------------------------- */

IssmDouble semic_albedo_slater(IssmDouble alb_snow, IssmDouble tsurf, IssmDouble tmin,
                           IssmDouble tmax, IssmDouble alb_smax, IssmDouble alb_smin)
{
	(void)alb_snow; /* the Fortran subroutine ignores the input alb_snow   */
	(void)tmax;     /* tmax is not actually used inside the formulation    */
	IssmDouble tm = 0.0;
	IssmDouble f  = 1.0 / (t0 - tmin);
	if (tsurf >= tmin && tsurf <= t0)
		tm = f * (tsurf - tmin);
	if (tsurf > t0)
		tm = 1.0;
	return alb_smax - (alb_smax - alb_smin) * tm * tm * tm;
}

IssmDouble semic_albedo_denby(IssmDouble melt, IssmDouble alb_smax,
                          IssmDouble alb_smin, IssmDouble mcrit)
{
	return alb_smin + (alb_smax - alb_smin) * std::exp(-melt / mcrit);
}

IssmDouble semic_albedo_isba(IssmDouble alb, IssmDouble sf, IssmDouble melt, IssmDouble tstic,
                         IssmDouble tau_a, IssmDouble tau_f, IssmDouble w_crit,
                         IssmDouble mcrit, IssmDouble alb_smin, IssmDouble alb_smax)
{
	/* dry case: linear decline                                            */
	IssmDouble alb_dry = alb - tau_a * tstic / tstic; /* tau [1/day], step [s] */
	/* wet case: exponential decline                                       */
	IssmDouble alb_wet = (alb - alb_smin) * std::exp(-tau_f * tstic / tstic) + alb_smin;
	IssmDouble alb_new = sf * tstic / (w_crit / rhow) * (alb_smax - alb_smin);

	IssmDouble w_alb = 0.0;
	if (melt > 0.0) w_alb = 1.0 - melt / mcrit;
	w_alb = std::min(1.0, std::max(w_alb, 0.0));

	IssmDouble alb_out = (1.0 - w_alb) * alb_dry + w_alb * alb_wet + alb_new;
	return std::min(alb_smax, std::max(alb_out, alb_smin));
}

/* ======================================================================== */
/*  SemicParam constructor – default values from run_semic_transient.f90    */
/* ======================================================================== */
SemicParam::SemicParam()
	: nx(1), n_ksub(3),
	  ceff(2.0e6),
	  albi(0.41), albl(0.07),
	  alb_smax(0.79), alb_smin(0.6),
	  hcrit(0.028), rcrit(0.85),
	  amp(3.0),
	  csh(2.0e-3), clh(5.0e-4),
	  tmin(-999.0), tmax(273.15),
	  tstic(86400.0), tsticsub(86400.0 / 3.0),
	  tau_a(0.008), tau_f(0.24),
	  w_crit(15.0), mcrit(6.0e-8),
	  afac(0.0), tmid(273.15),
	  alb_scheme(SEMIC_ALB_NONE)
{}

/* ======================================================================== */
/*  SemicBnd constructor                                                     */
/* ======================================================================== */
SemicBnd::SemicBnd()
	: t2m(false), tsurf(false), hsnow(false), alb(false),
	  melt(false), refr(false), smb(false), acc(false),
	  lhf(false), shf(false), subl(false), amp(false)
{}

/* ======================================================================== */
/*  Memory management                                                        */
/* ======================================================================== */
void semic_alloc(SemicState& s, int npts) {
	s.t2m        .assign(npts, 0.0);
	s.tsurf      .assign(npts, 0.0);
	s.hsnow      .assign(npts, 0.0);
	s.hice       .assign(npts, 0.0);
	s.alb        .assign(npts, 0.0);
	s.alb_snow   .assign(npts, 0.0);
	s.melt       .assign(npts, 0.0);
	s.melted_snow.assign(npts, 0.0);
	s.melted_ice .assign(npts, 0.0);
	s.refr       .assign(npts, 0.0);
	s.smb        .assign(npts, 0.0);
	s.acc        .assign(npts, 0.0);
	s.lhf        .assign(npts, 0.0);
	s.shf        .assign(npts, 0.0);
	s.lwu        .assign(npts, 0.0);
	s.subl       .assign(npts, 0.0);
	s.evap       .assign(npts, 0.0);
	s.smb_snow   .assign(npts, 0.0);
	s.smb_ice    .assign(npts, 0.0);
	s.runoff     .assign(npts, 0.0);
	s.qmr        .assign(npts, 0.0);
	s.qmr_res    .assign(npts, 0.0);
	s.amp        .assign(npts, 0.0);
	/* forcing */
	s.sf  .assign(npts, 0.0);
	s.rf  .assign(npts, 0.0);
	s.sp  .assign(npts, 0.0);
	s.lwd .assign(npts, 0.0);
	s.swd .assign(npts, 0.0);
	s.wind.assign(npts, 0.0);
	s.rhoa.assign(npts, 0.0);
	s.qq  .assign(npts, 0.0);
	s.mask.assign(npts, 0);
}

void semic_dealloc(SemicState& s) {
	/* std::vector destructor handles memory; just clear to release */
	s.t2m        .clear();  s.tsurf      .clear();
	s.hsnow      .clear();  s.hice       .clear();
	s.alb        .clear();  s.alb_snow   .clear();
	s.melt       .clear();  s.melted_snow.clear();
	s.melted_ice .clear();  s.refr       .clear();
	s.smb        .clear();  s.acc        .clear();
	s.lhf        .clear();  s.shf        .clear();
	s.lwu        .clear();  s.subl       .clear();
	s.evap       .clear();  s.smb_snow   .clear();
	s.smb_ice    .clear();  s.runoff     .clear();
	s.qmr        .clear();  s.qmr_res    .clear();
	s.amp        .clear();
	s.sf  .clear();  s.rf  .clear();  s.sp  .clear();
	s.lwd .clear();  s.swd .clear();  s.wind.clear();
	s.rhoa.clear();  s.qq  .clear();  s.mask.clear();
}

/* ======================================================================== */
/*  Boundary helper                                                          */
/* ======================================================================== */
void semic_boundary_define(SemicBnd& bnd, const std::vector<std::string>& names)
{
	bnd = SemicBnd();   /* reset all to false */
	for (const auto& nm : names) {
		if      (nm == "t2m"  ) bnd.t2m   = true;
		else if (nm == "tsurf") bnd.tsurf = true;
		else if (nm == "hsnow") bnd.hsnow = true;
		else if (nm == "alb"  ) bnd.alb   = true;
		else if (nm == "melt" ) bnd.melt  = true;
		else if (nm == "refr" ) bnd.refr  = true;
		else if (nm == "smb"  ) bnd.smb   = true;
		else if (nm == "acc"  ) bnd.acc   = true;
		else if (nm == "lhf"  ) bnd.lhf   = true;
		else if (nm == "shf"  ) bnd.shf   = true;
		else if (nm == "subl" ) bnd.subl  = true;
		else if (nm == "amp"  ) bnd.amp   = true;
		/* unknown names are silently ignored (matches Fortran DEFAULT case) */
	}
}

/* ======================================================================== */
/*  Energy balance                                                           */
/* ======================================================================== */
void semic_energy_balance(SemicState& now, const SemicParam& par,
                          const SemicBnd& bnd, int /*day*/, int /*year*/)
{
	const int nx = par.nx;

	/* 1. Sensible heat flux */
	if (!bnd.shf) {
		for (int i = 0; i < nx; i++)
			now.shf[i] = semic_sensible_heat_flux(
				now.tsurf[i], now.t2m[i], now.wind[i], now.rhoa[i],
				par.csh, cap);
	}

	/* 2. Latent heat flux (sublimation / evaporation) */
	if (!bnd.lhf) {
		for (int i = 0; i < nx; i++) {
			now.subl[i] = 0.0;
			now.evap[i] = 0.0;
			now.lhf[i]  = 0.0;
			semic_latent_heat_flux(
				now.tsurf[i], now.wind[i], now.qq[i], now.sp[i],
				now.rhoa[i], now.mask[i], par.clh,
				now.lhf[i], now.subl[i], now.evap[i]);
		}
	}
	/* sublimation: kg/(s m2) -> m/s */
	for (int i = 0; i < nx; i++)
		now.subl[i] /= rhow;

	/* 3. Upwelling longwave radiation */
	for (int i = 0; i < nx; i++)
		now.lwu[i] = semic_longwave_upward(now.tsurf[i]);

	/* 4. Surface energy balance */
	std::vector<IssmDouble> qsb(nx);
	for (int i = 0; i < nx; i++)
		qsb[i] = (1.0 - now.alb[i]) * now.swd[i]
		        + now.lwd[i] - now.lwu[i]
		        - now.shf[i] - now.lhf[i]
		        - now.qmr_res[i];

	/* 5. Update surface temperature */
	for (int i = 0; i < nx; i++)
		now.qmr[i] = 0.0;

	if (!bnd.tsurf) {
		for (int i = 0; i < nx; i++) {
			now.tsurf[i] += qsb[i] * par.tsticsub / par.ceff;
			/* store residual energy when tsurf > t0 over ice or snow */
			if ((now.mask[i] == 2 || now.hsnow[i] > 0.0) && now.tsurf[i] > t0) {
				now.qmr[i]   = (now.tsurf[i] - t0) * par.ceff / par.tsticsub;
				now.tsurf[i] = t0;
			}
		}
	}

	/* 6. Update 2-m air temperature over ice/snow */
	for (int i = 0; i < nx; i++) {
		if (now.mask[i] == 2 || now.hsnow[i] > 0.0)
			now.t2m[i] += (now.shf[i] + now.lhf[i]) * par.tsticsub / par.ceff;
	}
}

/* ======================================================================== */
/*  Mass balance                                                             */
/* ======================================================================== */
void semic_mass_balance(SemicState& now, const SemicParam& par,
                        const SemicBnd& bnd, int /*day*/, int /*year*/)
{
	const int nx = par.nx;
	const IssmDouble epsil = std::numeric_limits<IssmDouble>::epsilon();

	/* Temporary work arrays */
	std::vector<IssmDouble> qmelt(nx, 0.0), qcold(nx, 0.0);
	std::vector<IssmDouble> above(nx, 0.0), below(nx, 0.0);
	std::vector<IssmDouble> f_rz(nx, 0.0), f_alb(nx, 0.0);
	std::vector<IssmDouble> refrozen_rain(nx, 0.0), refrozen_snow(nx, 0.0);
	std::vector<IssmDouble> snow_to_ice(nx, 0.0);

	/* 1. Diurnal cycle amplitude */
	for (int i = 0; i < nx; i++) {
		if (!bnd.amp)
			now.amp[i] = par.amp;
		semic_diurnal_cycle(
			now.amp[i],
			now.tsurf[i] - t0 + now.qmr[i] / par.ceff * par.tstic,
			above[i], below[i]);
		now.qmr[i] = 0.0;
	}

	/* 2-3. Melt and cold energy where mask >= 1 */
	for (int i = 0; i < nx; i++) {
		if (now.mask[i] >= 1) {
			qmelt[i] = std::max(0.0, above[i] * par.ceff / par.tstic);
			qcold[i] = std::max(0.0, std::fabs(below[i]) * par.ceff / par.tstic);
		}
	}

	/* 4. Potential melt and split into snow/ice melt */
	if (!bnd.melt) {
		for (int i = 0; i < nx; i++)
			now.melt[i] = qmelt[i] / (rhow * clm);
	}
	for (int i = 0; i < nx; i++) {
		now.melted_snow[i] = std::min(now.melt[i], now.hsnow[i] / par.tstic);
		now.melted_ice[i]  = now.melt[i] - now.melted_snow[i];
	}
	if (!bnd.melt) {
		for (int i = 0; i < nx; i++) {
			if (now.mask[i] == 2) {
				now.melt[i] = now.melted_snow[i] + now.melted_ice[i];
			} else {
				now.melt[i]        = now.melted_snow[i];
				now.melted_ice[i]  = 0.0;
			}
		}
	}

	/* 5. Refreezing */
	if (!bnd.refr) {
		for (int i = 0; i < nx; i++) {
			f_rz[i] = par.rcrit;
			IssmDouble pot_refr = qcold[i] / (rhow * clm);
			refrozen_rain[i] = std::min(pot_refr, now.rf[i]);
			refrozen_snow[i] = std::max(pot_refr - refrozen_rain[i], 0.0);
			refrozen_snow[i] = std::min(refrozen_snow[i], now.melted_snow[i]);
			refrozen_rain[i] = f_rz[i] * std::min(now.hsnow[i] / par.tstic, refrozen_rain[i]);
			refrozen_snow[i] = f_rz[i] * std::min(now.hsnow[i] / par.tstic, refrozen_snow[i]);
			now.refr[i]      = refrozen_rain[i] + refrozen_snow[i];
			/* energy released during refreezing that is not used */
			now.qmr[i] -= (1.0 - f_rz[i]) * now.refr[i] * rhow * clm;
		}
	}

	/* 6. Runoff */
	for (int i = 0; i < nx; i++)
		now.runoff[i] = now.melt[i] + now.rf[i] - refrozen_rain[i];

	/* 7. Accumulation */
	if (!bnd.acc) {
		for (int i = 0; i < nx; i++)
			now.acc[i] = now.sf[i] - now.subl[i] + now.refr[i];
	}

	/* 8. SMB of snow */
	if (!bnd.smb) {
		for (int i = 0; i < nx; i++)
			now.smb_snow[i] = now.sf[i] - now.subl[i]
			                 - now.melted_snow[i] + refrozen_snow[i];
	}

	/* 9. Update snow height */
	for (int i = 0; i < nx; i++) {
		if (now.mask[i] == 0)
			now.hsnow[i] = 0.0;
		else
			now.hsnow[i] = std::max(0.0, now.hsnow[i] + now.smb_snow[i] * par.tstic);
	}

	/* 10. Cap snow height (excess becomes ice) */
	for (int i = 0; i < nx; i++) {
		snow_to_ice[i]   = std::max(0.0, now.hsnow[i] - hsmax);
		now.hsnow[i]    -= snow_to_ice[i];
		now.smb_ice[i]   = snow_to_ice[i] / par.tstic
		                 - now.melted_ice[i] + refrozen_rain[i];
		now.hice[i]     += now.smb_ice[i] * par.tstic;
	}

	/* 11. Total SMB */
	if (!bnd.smb) {
		for (int i = 0; i < nx; i++) {
			if (now.mask[i] == 2)
				now.smb[i] = now.smb_snow[i] + now.smb_ice[i]
				           - snow_to_ice[i] / par.tstic;
			else
				now.smb[i] = now.smb_snow[i]
				           + std::max(0.0, now.smb_ice[i] - snow_to_ice[i] / par.tstic);
		}
	}

	/* 12. Snow albedo update */
	for (int i = 0; i < nx; i++)
		f_alb[i] = 1.0 - std::exp(-now.hsnow[i] / (par.hcrit + epsil));

	if (!bnd.alb) {
		switch (par.alb_scheme) {
			case SEMIC_ALB_SLATER:
				for (int i = 0; i < nx; i++)
					now.alb_snow[i] = semic_albedo_slater(
						now.alb_snow[i], now.tsurf[i],
						par.tmin, par.tmax, par.alb_smax, par.alb_smin);
				break;
			case SEMIC_ALB_DENBY:
				for (int i = 0; i < nx; i++)
					now.alb_snow[i] = semic_albedo_denby(
						now.melt[i], par.alb_smax, par.alb_smin, par.mcrit);
				break;
			case SEMIC_ALB_ISBA:
				for (int i = 0; i < nx; i++)
					now.alb_snow[i] = semic_albedo_isba(
						now.alb_snow[i], now.sf[i], now.melt[i],
						par.tstic, par.tau_a, par.tau_f, par.w_crit, par.mcrit,
						par.alb_smin, par.alb_smax);
				break;
			case SEMIC_ALB_NONE:
				/* "none" (lowercase): fix alb_snow to the maximum value each step.
				 * Matches Fortran: trim(par%alb_scheme) .eq. "none" => alb_snow = alb_smax */
				for (int i = 0; i < nx; i++)
					now.alb_snow[i] = par.alb_smax;
				break;
			case SEMIC_ALB_SKIP:
			default:
				/* "None" (capital, annual driver) or unrecognised: leave alb_snow
				 * unchanged.  Matches Fortran run_semic where par%alb_scheme="None"
				 * falls through all branches without updating alb_snow. */
				break;
		}

		/* Grid-averaged albedo by mask */
		for (int i = 0; i < nx; i++) {
			if      (now.mask[i] == 2)
				now.alb[i] = par.albi + f_alb[i] * (now.alb_snow[i] - par.albi);
			else if (now.mask[i] == 1)
				now.alb[i] = par.albl + f_alb[i] * (now.alb_snow[i] - par.albl);
			else
				now.alb[i] = 0.06;
		}

		/* "alex" scheme overrides both alb_snow and alb */
		if (par.alb_scheme == SEMIC_ALB_ALEX) {
			for (int i = 0; i < nx; i++) {
				now.alb_snow[i] = par.alb_smin
				    + (par.alb_smax - par.alb_smin)
				    * (0.5 * std::tanh(par.afac * (now.t2m[i] - par.tmid)) + 0.5);
				now.alb[i] = now.alb_snow[i];
			}
		}
	}

	/* Store residual energy */
	for (int i = 0; i < nx; i++) {
		if (now.mask[i] == 0)
			now.qmr_res[i] = 0.0;
		else
			now.qmr_res[i] = now.qmr[i];
	}
}

/* ======================================================================== */
/*  Combined surface energy-and-mass balance step                           */
/* ======================================================================== */
void semic_surface_step(SemicState& now, const SemicParam& par,
                        const SemicBnd& bnd, int day, int year)
{
	/* quasi-relaxation loop for energy balance (no hsnow update inside) */
	for (int ksub = 0; ksub < par.n_ksub; ksub++)
		semic_energy_balance(now, par, bnd, day, year);
	semic_mass_balance(now, par, bnd, day, year);
}

/* ======================================================================== */
/*  RunSemic – annual driver (replaces Fortran run_semic_)                  */
/* ======================================================================== */
void RunSemic(const IssmDouble* sf_in,   const IssmDouble* rf_in,
              const IssmDouble* swd_in,  const IssmDouble* lwd_in,
              const IssmDouble* wind_in, const IssmDouble* sp_in,
              const IssmDouble* rhoa_in, const IssmDouble* qq_in,
              const IssmDouble* tt_in,
              IssmDouble& tsurf_out, IssmDouble& smb_out,
              IssmDouble& saccu_out, IssmDouble& smelt_out)
{
	const int nloop = 10;
	const int nx    = 1;
	const int ntime = 365;

	SemicParam par;
	par.nx     = nx;
	par.tstic  = 86400.0;
	par.ceff   = 2.0e6;
	par.csh    = 2.0e-3;
	par.clh    = 5.0e-4;
	par.alb_smax  = 0.79;
	par.alb_smin  = 0.6;
	par.albi      = 0.41;
	par.albl      = 0.07;
	par.tmin      = -999.0;
	par.tmax      = 273.15;
	par.hcrit     = 0.028;
	par.rcrit     = 0.85;
	par.amp       = 3.0;
	/* Fortran run_semic sets par%alb_scheme="None" (capital N).
	 * The Fortran mass_balance only matches lowercase "none", so alb_snow is
	 * never updated inside that routine — it keeps its initial value of 0.8.
	 * Use SEMIC_ALB_SKIP to replicate this fall-through behaviour. */
	par.alb_scheme= SEMIC_ALB_SKIP;
	par.tau_a     = 0.008;
	par.tau_f     = 0.24;
	par.w_crit    = 15.0;
	par.mcrit     = 6.0e-8;
	par.n_ksub    = 3;
	par.tsticsub  = par.tstic / par.n_ksub;

	SemicState now;
	semic_alloc(now, nx);

	now.mask [0]      = 2;
	now.hsnow[0]      = 1.0;
	now.hice [0]      = 0.0;
	now.alb  [0]      = 0.8;
	now.tsurf[0]      = 260.0;
	now.alb_snow[0]   = 0.8;
	now.qmr_res[0]    = 0.0;

	SemicBnd bnd;   /* all false */
	std::vector<std::string> bnames; /* empty */
	semic_boundary_define(bnd, bnames);

	tsurf_out = 0.0;
	smb_out   = 0.0;
	saccu_out = 0.0;
	smelt_out = 0.0;

	for (int k = 0; k < nloop; k++) {
		for (int i = 0; i < ntime; i++) {
			now.sf  [0] = sf_in  [i];
			now.rf  [0] = rf_in  [i];
			now.sp  [0] = sp_in  [i];
			now.lwd [0] = lwd_in [i];
			now.swd [0] = swd_in [i];
			now.wind[0] = wind_in[i];
			now.rhoa[0] = rhoa_in[i];
			now.t2m [0] = tt_in  [i];
			now.qq  [0] = qq_in  [i];

			semic_surface_step(now, par, bnd, i + 1, 0);

			if (k == nloop - 1) {
				tsurf_out += now.tsurf[0] / 365.0;
				smb_out   += now.smb[0]   / 365.0;
				saccu_out += now.alb[0]   / 365.0;  /* matches Fortran labelling */
				smelt_out += now.melt[0]  / 365.0;
			}
		}
	}

	semic_dealloc(now);
}

/* ======================================================================== */
/*  RunSemicTransient – transient driver (replaces Fortran                  */
/*                       run_semic_transient_)                              */
/* ======================================================================== */
void RunSemicTransient(int nx, int ntime, int nloop,
                       const IssmDouble* sf_in,    const IssmDouble* rf_in,
                       const IssmDouble* swd_in,   const IssmDouble* lwd_in,
                       const IssmDouble* wind_in,  const IssmDouble* sp_in,
                       const IssmDouble* rhoa_in,  const IssmDouble* qq_in,
                       const IssmDouble* tt_in,
                       const IssmDouble* tsurf_in, const IssmDouble* qmr_in,
                       IssmDouble tstic,
                       IssmDouble hcrit, IssmDouble rcrit,
                       const IssmDouble* mask,   const IssmDouble* hice,
                       const IssmDouble* hsnow,
                       const IssmDouble* albedo, const IssmDouble* albedo_snow,
                       int alb_scheme_int,
                       IssmDouble alb_smax, IssmDouble alb_smin,
                       IssmDouble albi,     IssmDouble albl,
                       const IssmDouble* Tamp,
                       IssmDouble tmin, IssmDouble tmax, IssmDouble tmid,
                       IssmDouble mcrit, IssmDouble wcrit,
                       IssmDouble tau_a, IssmDouble tau_f, IssmDouble afac,
                       bool verbose,
                       IssmDouble* tsurf_out,    IssmDouble* smb_out,
                       IssmDouble* smbi_out,     IssmDouble* smbs_out,
                       IssmDouble* saccu_out,    IssmDouble* smelt_out,
                       IssmDouble* refr_out,     IssmDouble* alb_out,
                       IssmDouble* alb_snow_out,
                       IssmDouble* hsnow_out,    IssmDouble* hice_out,
                       IssmDouble* qmr_out,
                       IssmDouble* runoff_out,   IssmDouble* subl_out)
{
	if (verbose) {
		/* minimal printf – mirrors the Fortran debug prints */
		_printf0_("run_semic_transient: ntime=" << ntime << "  nx=" << nx << "\n");
	}

	SemicParam par;
	par.nx        = nx;
	par.tstic     = tstic;
	par.ceff      = 2.0e6;
	par.csh       = 2.0e-3;
	par.clh       = 5.0e-4;
	par.alb_smax  = alb_smax;
	par.alb_smin  = alb_smin;
	par.albi      = albi;
	par.albl      = albl;
	par.tmin      = tmin;
	par.tmax      = tmax;
	par.hcrit     = hcrit;
	par.rcrit     = rcrit;
	par.tau_a     = tau_a;
	par.tau_f     = tau_f;
	par.w_crit    = wcrit;
	par.mcrit     = mcrit;
	par.afac      = afac;
	par.tmid      = tmid;
	par.n_ksub    = 3;
	par.tsticsub  = par.tstic / par.n_ksub;
	/* Fortran: surface%par%amp = Tamp  (scalar — reads first element only).
	 * bnd%amp is FALSE so semic_mass_balance uses par.amp, not now.amp[i]. */
	par.amp       = Tamp[0];

	switch (alb_scheme_int) {
		case 0: par.alb_scheme = SEMIC_ALB_NONE;   break;
		case 1: par.alb_scheme = SEMIC_ALB_SLATER; break;
		case 2: par.alb_scheme = SEMIC_ALB_DENBY;  break;
		case 3: par.alb_scheme = SEMIC_ALB_ISBA;   break;
		case 4: par.alb_scheme = SEMIC_ALB_ALEX;   break;
		default:
			_error_("RunSemicTransient: unknown albedo scheme " << alb_scheme_int);
	}

	SemicState now;
	semic_alloc(now, nx);

	/* Initialise from inputs */
	for (int i = 0; i < nx; i++) {
		now.mask    [i] = reCast<int,IssmDouble>(mask[i]);
		now.hsnow   [i] = hsnow   [i];
		now.hice    [i] = hice    [i];
		now.tsurf   [i] = tsurf_in[i];
		now.alb     [i] = albedo  [i];
		now.alb_snow[i] = albedo_snow[i];
		now.qmr_res [i] = qmr_in  [i];
		now.amp     [i] = Tamp    [i];
	}

	SemicBnd bnd;
	std::vector<std::string> bnames;
	semic_boundary_define(bnd, bnames);

	for (int k = 0; k < nloop; k++) {
		for (int i = 0; i < ntime; i++) {
			for (int j = 0; j < nx; j++) {
				now.sf  [j] = sf_in  [j];
				now.rf  [j] = rf_in  [j];
				now.sp  [j] = sp_in  [j];
				now.lwd [j] = lwd_in [j];
				now.swd [j] = swd_in [j];
				now.wind[j] = wind_in[j];
				now.rhoa[j] = rhoa_in[j];
				now.t2m [j] = tt_in  [j];
				now.qq  [j] = qq_in  [j];
				now.qmr_res[j] = qmr_in[j];
			}

			semic_surface_step(now, par, bnd, i + 1, 0);

			if (k == nloop - 1) {
				for (int j = 0; j < nx; j++) {
					tsurf_out  [j] = now.tsurf  [j];
					smb_out    [j] = now.smb     [j];
					smbi_out   [j] = now.smb_ice [j];
					smbs_out   [j] = now.smb_snow[j];
					saccu_out  [j] = now.acc     [j];
					smelt_out  [j] = now.melt    [j];
					refr_out   [j] = now.refr    [j];
					alb_out    [j] = now.alb     [j];
					alb_snow_out[j]= now.alb_snow[j];
					hsnow_out  [j] = now.hsnow   [j];
					hice_out   [j] = now.hice    [j];
					qmr_out    [j] = now.qmr_res [j];
					runoff_out [j] = now.runoff  [j];
					subl_out   [j] = now.subl    [j];
				}
			}
		}
	}

	semic_dealloc(now);
}
