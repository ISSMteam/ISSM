/*!\file semic.h
 * \brief C++ port of the SEMIC (Simple Energy/Mass Balance Ice-sheet and
 *        Climate model) surface physics originally written in Fortran by
 *        Mario Krapp.  All physics are translated verbatim from
 *        externalpackages/semic/src/surface_physics.f90 and the two driver
 *        subroutines run_semic.f90 / run_semic_transient.f90 that live in
 *        src/c/modules/SurfaceMassBalancex/.
 *
 *  Public API
 *  ----------
 *  SemicParam     – physical / numerical parameters  (≡ surface_param_class)
 *  SemicState     – prognostic + diagnostic arrays   (≡ surface_state_class)
 *  SemicBnd       – boundary-override flags          (≡ boundary_opt_class)
 *
 *  semic_alloc / semic_dealloc  – allocate / free state arrays
 *  semic_boundary_define        – parse boundary list into SemicBnd flags
 *  semic_energy_balance         – one sub-step of the energy balance
 *  semic_mass_balance           – mass balance (called once per full step)
 *  semic_surface_step           – combined energy+mass (≡ surface_energy_and_mass_balance)
 *
 *  RunSemic          – annual driver  (replaces Fortran run_semic_)
 *  RunSemicTransient – transient driver (replaces Fortran run_semic_transient_)
 */

#ifndef _SEMIC_H_
#define _SEMIC_H_

#include <string>
#include <vector>
#include "../Numerics/types.h"

/* ======================================================================== */
/*  Albedo-scheme integer codes (matches run_semic_transient.f90 convention) */
/* ======================================================================== */
enum SemicAlbedoScheme {
	SEMIC_ALB_SKIP   = -1, /* internal: leave alb_snow unchanged (annual driver "None") */
	SEMIC_ALB_NONE   = 0,  /* "none" lowercase: set alb_snow = alb_smax each step      */
	SEMIC_ALB_SLATER = 1,
	SEMIC_ALB_DENBY  = 2,
	SEMIC_ALB_ISBA   = 3,
	SEMIC_ALB_ALEX   = 4
};

/* ======================================================================== */
/*  Physical constants (≡ module-level parameters in surface_physics.f90)   */
/* ======================================================================== */
namespace semic_const {
	static const IssmDouble pi   = 3.141592653589793238462643;
	static const IssmDouble t0   = 273.15;          /* melting point [K]               */
	static const IssmDouble sigm = 5.67e-8;         /* Stefan-Boltzmann [W/(m2 K4)]    */
	static const IssmDouble eps  = 0.62197;         /* molar-weight ratio water/dry air */
	static const IssmDouble cls  = 2.83e6;          /* latent heat sublimation [J/kg]  */
	static const IssmDouble clm  = 3.30e5;          /* latent heat melting [J/kg]      */
	static const IssmDouble clv  = 2.5e6;           /* latent heat condensation [J/kg] */
	static const IssmDouble cap  = 1000.0;          /* specific heat capacity air [J/(kg K)] */
	static const IssmDouble rhow = 1000.0;          /* density of water [kg/m3]        */
	static const IssmDouble hsmax= 5.0;             /* maximum snow height [m]         */
}

/* ======================================================================== */
/*  Parameter struct                                                         */
/* ======================================================================== */
struct SemicParam {
	int    nx;           /* number of grid points                          */
	int    n_ksub;       /* number of sub-daily time steps                 */
	IssmDouble ceff;         /* surface heat capacity snow/ice [J/(K m2)]      */
	IssmDouble albi;         /* bare-ice albedo                                */
	IssmDouble albl;         /* bare-land albedo                               */
	IssmDouble alb_smax;     /* maximum (fresh) snow albedo                    */
	IssmDouble alb_smin;     /* minimum (old/wet) snow albedo                  */
	IssmDouble hcrit;        /* critical snow height for 50 % snow cover [m]   */
	IssmDouble rcrit;        /* critical snow height for 50 % refreezing [m]   */
	IssmDouble amp;          /* diurnal cycle amplitude [K]                    */
	IssmDouble csh;          /* sensible heat exchange coefficient             */
	IssmDouble clh;          /* latent heat exchange coefficient               */
	IssmDouble tmin;         /* min temperature for Slater albedo decline [K]  */
	IssmDouble tmax;         /* max temperature for Slater albedo decline [K]  */
	IssmDouble tstic;        /* time step [s]                                  */
	IssmDouble tsticsub;     /* sub-time step [s]  (= tstic / n_ksub)          */
	IssmDouble tau_a;        /* dry albedo decline (ISBA) [1/day]              */
	IssmDouble tau_f;        /* wet albedo decline (ISBA) [1/day]              */
	IssmDouble w_crit;       /* critical liquid water (ISBA) [kg/m2]           */
	IssmDouble mcrit;        /* critical melt rate (ISBA / Denby) [m/s]        */
	IssmDouble afac;         /* "alex" albedo param                            */
	IssmDouble tmid;         /* "alex" albedo param [K]                        */
	SemicAlbedoScheme alb_scheme; /* albedo parameterisation to use        */

	SemicParam();        /* constructor sets sensible defaults             */
};

/* ======================================================================== */
/*  Boundary-override flag struct                                            */
/* ======================================================================== */
struct SemicBnd {
	bool t2m;
	bool tsurf;
	bool hsnow;
	bool alb;
	bool melt;
	bool refr;
	bool smb;
	bool acc;
	bool lhf;
	bool shf;
	bool subl;
	bool amp;

	SemicBnd();  /* constructor – all false */
};

/* ======================================================================== */
/*  State struct                                                             */
/* ======================================================================== */
struct SemicState {
	/* prognostic / diagnostic */
	std::vector<IssmDouble> t2m;          /* 2-m air temperature [K]                   */
	std::vector<IssmDouble> tsurf;        /* surface temperature [K]                   */
	std::vector<IssmDouble> hsnow;        /* snow-pack height (water-equiv.) [m]        */
	std::vector<IssmDouble> hice;         /* ice thickness (water-equiv.) [m]          */
	std::vector<IssmDouble> alb;          /* grid-averaged albedo                       */
	std::vector<IssmDouble> alb_snow;     /* snow albedo                                */
	std::vector<IssmDouble> melt;         /* potential surface melt [m/s]              */
	std::vector<IssmDouble> melted_snow;  /* actual melted snow [m/s]                  */
	std::vector<IssmDouble> melted_ice;   /* actual melted ice [m/s]                   */
	std::vector<IssmDouble> refr;         /* refreezing [m/s]                          */
	std::vector<IssmDouble> smb;          /* surface mass balance [m/s]                */
	std::vector<IssmDouble> acc;          /* accumulation [m/s]                        */
	std::vector<IssmDouble> lhf;          /* latent heat flux [W/m2]                   */
	std::vector<IssmDouble> shf;          /* sensible heat flux [W/m2]                 */
	std::vector<IssmDouble> lwu;          /* upwelling longwave radiation [W/m2]       */
	std::vector<IssmDouble> subl;         /* sublimation [m/s]                         */
	std::vector<IssmDouble> evap;         /* evaporation                               */
	std::vector<IssmDouble> smb_snow;     /* SMB of snow [m/s]                         */
	std::vector<IssmDouble> smb_ice;      /* SMB of ice [m/s]                          */
	std::vector<IssmDouble> runoff;       /* surface runoff [m/s]                      */
	std::vector<IssmDouble> qmr;          /* heat flux from melt/refreezing [W/m2]     */
	std::vector<IssmDouble> qmr_res;      /* residual heat flux [W/m2]                 */
	std::vector<IssmDouble> amp;          /* diurnal cycle amplitude [K]               */
	/* forcing */
	std::vector<IssmDouble> sf;           /* snowfall [m/s]                            */
	std::vector<IssmDouble> rf;           /* rainfall [m/s]                            */
	std::vector<IssmDouble> sp;           /* surface pressure [Pa]                     */
	std::vector<IssmDouble> lwd;          /* downwelling longwave radiation [W/m2]     */
	std::vector<IssmDouble> swd;          /* downwelling shortwave radiation [W/m2]    */
	std::vector<IssmDouble> wind;         /* surface wind speed [m/s]                  */
	std::vector<IssmDouble> rhoa;         /* air density [kg/m3]                       */
	std::vector<IssmDouble> qq;           /* air specific humidity [kg/kg]             */
	std::vector<int>    mask;         /* ocean/land/ice mask  0/1/2                */
};

/* ======================================================================== */
/*  Memory management                                                        */
/* ======================================================================== */
void semic_alloc(SemicState& s, int npts);
void semic_dealloc(SemicState& s);

/* ======================================================================== */
/*  Boundary helper                                                          */
/* ======================================================================== */
/*! Parse a list of variable names and set the corresponding override flags.
 *  \param bnd     output boundary flag struct
 *  \param names   list of variable names (up to 30, empty strings ignored) */
void semic_boundary_define(SemicBnd& bnd, const std::vector<std::string>& names);

/* ======================================================================== */
/*  Core physics routines                                                    */
/* ======================================================================== */
void semic_energy_balance(SemicState& now, const SemicParam& par, const SemicBnd& bnd,
                          int day, int year);

void semic_mass_balance  (SemicState& now, const SemicParam& par, const SemicBnd& bnd,
                          int day, int year);

/*! Combined energy+mass balance step (≡ surface_energy_and_mass_balance) */
void semic_surface_step  (SemicState& now, const SemicParam& par, const SemicBnd& bnd,
                          int day, int year);

/* ======================================================================== */
/*  Elemental physics helpers (scalar, inlined equivalents of Fortran        */
/*  ELEMENTAL subroutines — also useful for unit tests)                      */
/* ======================================================================== */
IssmDouble semic_sensible_heat_flux(IssmDouble ts, IssmDouble ta, IssmDouble wind, IssmDouble rhoa,
                                IssmDouble csh, IssmDouble cap_air);

void   semic_latent_heat_flux  (IssmDouble ts, IssmDouble wind, IssmDouble shum, IssmDouble sp,
                                IssmDouble rhoatm, int mask, IssmDouble clh,
                                IssmDouble& lhf, IssmDouble& subl, IssmDouble& evap);

IssmDouble semic_longwave_upward   (IssmDouble ts);

void   semic_diurnal_cycle     (IssmDouble amp, IssmDouble tmean, IssmDouble& above, IssmDouble& below);

IssmDouble semic_albedo_slater(IssmDouble alb_snow, IssmDouble tsurf, IssmDouble tmin, IssmDouble tmax,
                           IssmDouble alb_smax, IssmDouble alb_smin);
IssmDouble semic_albedo_denby (IssmDouble melt, IssmDouble alb_smax, IssmDouble alb_smin, IssmDouble mcrit);
IssmDouble semic_albedo_isba  (IssmDouble alb, IssmDouble sf, IssmDouble melt, IssmDouble tstic,
                           IssmDouble tau_a, IssmDouble tau_f, IssmDouble w_crit, IssmDouble mcrit,
                           IssmDouble alb_smin, IssmDouble alb_smax);

/* ======================================================================== */
/*  High-level drivers  (replace the two Fortran entry points)               */
/* ======================================================================== */

/*! Annual driver: loops 365 days × nloop spinup cycles for a single vertex.
 *  Input arrays have dimension [365].
 *  Used by Element::SmbSemic() (method 0).
 *
 *  \param sf_in     snowfall rate [m/s], 365 entries
 *  \param rf_in     rainfall rate [m/s], 365 entries
 *  \param swd_in    downwelling SW radiation [W/m2], 365 entries
 *  \param lwd_in    downwelling LW radiation [W/m2], 365 entries
 *  \param wind_in   surface wind speed [m/s], 365 entries
 *  \param sp_in     surface pressure [Pa], 365 entries
 *  \param rhoa_in   air density [kg/m3], 365 entries
 *  \param qq_in     specific humidity [kg/kg], 365 entries
 *  \param tt_in     2-m temperature [K], 365 entries
 *  \param tsurf_out mean annual surface temperature [K]        (output)
 *  \param smb_out   mean annual SMB [m/s]                      (output)
 *  \param saccu_out mean annual albedo (historically labelled) (output)
 *  \param smelt_out mean annual melt [m/s]                     (output)
 */
void RunSemic(const IssmDouble* sf_in, const IssmDouble* rf_in,
              const IssmDouble* swd_in, const IssmDouble* lwd_in,
              const IssmDouble* wind_in, const IssmDouble* sp_in,
              const IssmDouble* rhoa_in, const IssmDouble* qq_in,
              const IssmDouble* tt_in,
              IssmDouble& tsurf_out, IssmDouble& smb_out,
              IssmDouble& saccu_out, IssmDouble& smelt_out);

/*! Transient driver: advances one time step for nx grid points.
 *  Used by Element::SmbSemicTransient() (method 1).
 *
 *  Input/output arrays all have dimension [nx].
 *
 *  \param nx          number of grid points
 *  \param ntime       number of sub-time steps within a call (usually 1)
 *  \param nloop       number of spinup repetitions (usually 1)
 *  \param sf_in       snowfall   [m/s]
 *  \param rf_in       rainfall   [m/s]
 *  \param swd_in      downwelling SW [W/m2]
 *  \param lwd_in      downwelling LW [W/m2]
 *  \param wind_in     wind speed [m/s]
 *  \param sp_in       surface pressure [Pa]
 *  \param rhoa_in     air density [kg/m3]
 *  \param qq_in       specific humidity [kg/kg]
 *  \param tt_in       2-m temperature [K]
 *  \param tsurf_in    previous surface temperature [K]
 *  \param qmr_in      residual heat flux from previous step [W/m2]
 *  \param tstic       time step [s]
 *  \param hcrit       critical snow height for snow cover [m]
 *  \param rcrit       critical snow height for refreezing [m]
 *  \param mask        ocean/land/ice mask 0/1/2
 *  \param hice        ice thickness (water-equiv.) [m]   (in/out via out arrays)
 *  \param hsnow       snow height (water-equiv.) [m]     (in/out via out arrays)
 *  \param albedo      grid-averaged albedo
 *  \param albedo_snow snow albedo
 *  \param alb_scheme  albedo scheme code (SemicAlbedoScheme)
 *  \param alb_smax    maximum snow albedo
 *  \param alb_smin    minimum snow albedo
 *  \param albi        bare-ice albedo
 *  \param albl        bare-land albedo
 *  \param Tamp        diurnal cycle amplitude [K]
 *  \param tmin        min temperature for Slater scheme [K]
 *  \param tmax        max temperature for Slater scheme [K]
 *  \param tmid        "alex" parameter [K]
 *  \param mcrit       critical melt rate [m/s]
 *  \param wcrit       critical liquid water (ISBA) [kg/m2]
 *  \param tau_a       dry albedo decline (ISBA)
 *  \param tau_f       wet albedo decline (ISBA)
 *  \param afac        "alex" parameter
 *  \param verbose     print debug info
 *  \param tsurf_out … output arrays (same units as inputs)
 */
void RunSemicTransient(int nx, int ntime, int nloop,
                       const IssmDouble* sf_in,    const IssmDouble* rf_in,
                       const IssmDouble* swd_in,   const IssmDouble* lwd_in,
                       const IssmDouble* wind_in,  const IssmDouble* sp_in,
                       const IssmDouble* rhoa_in,  const IssmDouble* qq_in,
                       const IssmDouble* tt_in,
                       const IssmDouble* tsurf_in, const IssmDouble* qmr_in,
                       IssmDouble tstic,
                       IssmDouble hcrit, IssmDouble rcrit,
                       const IssmDouble* mask, const IssmDouble* hice, const IssmDouble* hsnow,
                       const IssmDouble* albedo, const IssmDouble* albedo_snow,
                       int alb_scheme,
                       IssmDouble alb_smax, IssmDouble alb_smin, IssmDouble albi, IssmDouble albl,
                       const IssmDouble* Tamp,
                       IssmDouble tmin, IssmDouble tmax, IssmDouble tmid,
                       IssmDouble mcrit, IssmDouble wcrit,
                       IssmDouble tau_a, IssmDouble tau_f, IssmDouble afac,
                       bool verbose,
                       IssmDouble* tsurf_out,  IssmDouble* smb_out,
                       IssmDouble* smbi_out,   IssmDouble* smbs_out,
                       IssmDouble* saccu_out,  IssmDouble* smelt_out,
                       IssmDouble* refr_out,   IssmDouble* alb_out,
                       IssmDouble* alb_snow_out,
                       IssmDouble* hsnow_out,  IssmDouble* hice_out,
                       IssmDouble* qmr_out,
                       IssmDouble* runoff_out, IssmDouble* subl_out);

#endif /* _SEMIC_H_ */
