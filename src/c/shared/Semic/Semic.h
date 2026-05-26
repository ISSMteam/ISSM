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
	static const double pi   = 3.141592653589793238462643;
	static const double t0   = 273.15;          /* melting point [K]               */
	static const double sigm = 5.67e-8;         /* Stefan-Boltzmann [W/(m2 K4)]    */
	static const double eps  = 0.62197;         /* molar-weight ratio water/dry air */
	static const double cls  = 2.83e6;          /* latent heat sublimation [J/kg]  */
	static const double clm  = 3.30e5;          /* latent heat melting [J/kg]      */
	static const double clv  = 2.5e6;           /* latent heat condensation [J/kg] */
	static const double cap  = 1000.0;          /* specific heat capacity air [J/(kg K)] */
	static const double rhow = 1000.0;          /* density of water [kg/m3]        */
	static const double hsmax= 5.0;             /* maximum snow height [m]         */
}

/* ======================================================================== */
/*  Parameter struct                                                         */
/* ======================================================================== */
struct SemicParam {
	int    nx;           /* number of grid points                          */
	int    n_ksub;       /* number of sub-daily time steps                 */
	double ceff;         /* surface heat capacity snow/ice [J/(K m2)]      */
	double albi;         /* bare-ice albedo                                */
	double albl;         /* bare-land albedo                               */
	double alb_smax;     /* maximum (fresh) snow albedo                    */
	double alb_smin;     /* minimum (old/wet) snow albedo                  */
	double hcrit;        /* critical snow height for 50 % snow cover [m]   */
	double rcrit;        /* critical snow height for 50 % refreezing [m]   */
	double amp;          /* diurnal cycle amplitude [K]                    */
	double csh;          /* sensible heat exchange coefficient             */
	double clh;          /* latent heat exchange coefficient               */
	double tmin;         /* min temperature for Slater albedo decline [K]  */
	double tmax;         /* max temperature for Slater albedo decline [K]  */
	double tstic;        /* time step [s]                                  */
	double tsticsub;     /* sub-time step [s]  (= tstic / n_ksub)          */
	double tau_a;        /* dry albedo decline (ISBA) [1/day]              */
	double tau_f;        /* wet albedo decline (ISBA) [1/day]              */
	double w_crit;       /* critical liquid water (ISBA) [kg/m2]           */
	double mcrit;        /* critical melt rate (ISBA / Denby) [m/s]        */
	double afac;         /* "alex" albedo param                            */
	double tmid;         /* "alex" albedo param [K]                        */
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
	std::vector<double> t2m;          /* 2-m air temperature [K]                   */
	std::vector<double> tsurf;        /* surface temperature [K]                   */
	std::vector<double> hsnow;        /* snow-pack height (water-equiv.) [m]        */
	std::vector<double> hice;         /* ice thickness (water-equiv.) [m]          */
	std::vector<double> alb;          /* grid-averaged albedo                       */
	std::vector<double> alb_snow;     /* snow albedo                                */
	std::vector<double> melt;         /* potential surface melt [m/s]              */
	std::vector<double> melted_snow;  /* actual melted snow [m/s]                  */
	std::vector<double> melted_ice;   /* actual melted ice [m/s]                   */
	std::vector<double> refr;         /* refreezing [m/s]                          */
	std::vector<double> smb;          /* surface mass balance [m/s]                */
	std::vector<double> acc;          /* accumulation [m/s]                        */
	std::vector<double> lhf;          /* latent heat flux [W/m2]                   */
	std::vector<double> shf;          /* sensible heat flux [W/m2]                 */
	std::vector<double> lwu;          /* upwelling longwave radiation [W/m2]       */
	std::vector<double> subl;         /* sublimation [m/s]                         */
	std::vector<double> evap;         /* evaporation                               */
	std::vector<double> smb_snow;     /* SMB of snow [m/s]                         */
	std::vector<double> smb_ice;      /* SMB of ice [m/s]                          */
	std::vector<double> runoff;       /* surface runoff [m/s]                      */
	std::vector<double> qmr;          /* heat flux from melt/refreezing [W/m2]     */
	std::vector<double> qmr_res;      /* residual heat flux [W/m2]                 */
	std::vector<double> amp;          /* diurnal cycle amplitude [K]               */
	/* forcing */
	std::vector<double> sf;           /* snowfall [m/s]                            */
	std::vector<double> rf;           /* rainfall [m/s]                            */
	std::vector<double> sp;           /* surface pressure [Pa]                     */
	std::vector<double> lwd;          /* downwelling longwave radiation [W/m2]     */
	std::vector<double> swd;          /* downwelling shortwave radiation [W/m2]    */
	std::vector<double> wind;         /* surface wind speed [m/s]                  */
	std::vector<double> rhoa;         /* air density [kg/m3]                       */
	std::vector<double> qq;           /* air specific humidity [kg/kg]             */
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
double semic_sensible_heat_flux(double ts, double ta, double wind, double rhoa,
                                double csh, double cap_air);

void   semic_latent_heat_flux  (double ts, double wind, double shum, double sp,
                                double rhoatm, int mask, double clh,
                                double& lhf, double& subl, double& evap);

double semic_longwave_upward   (double ts);

void   semic_diurnal_cycle     (double amp, double tmean, double& above, double& below);

double semic_albedo_slater(double alb_snow, double tsurf, double tmin, double tmax,
                           double alb_smax, double alb_smin);
double semic_albedo_denby (double melt, double alb_smax, double alb_smin, double mcrit);
double semic_albedo_isba  (double alb, double sf, double melt, double tstic,
                           double tau_a, double tau_f, double w_crit, double mcrit,
                           double alb_smin, double alb_smax);

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
void RunSemic(const double* sf_in, const double* rf_in,
              const double* swd_in, const double* lwd_in,
              const double* wind_in, const double* sp_in,
              const double* rhoa_in, const double* qq_in,
              const double* tt_in,
              double& tsurf_out, double& smb_out,
              double& saccu_out, double& smelt_out);

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
                       const double* sf_in,    const double* rf_in,
                       const double* swd_in,   const double* lwd_in,
                       const double* wind_in,  const double* sp_in,
                       const double* rhoa_in,  const double* qq_in,
                       const double* tt_in,
                       const double* tsurf_in, const double* qmr_in,
                       double tstic,
                       double hcrit, double rcrit,
                       const double* mask, const double* hice, const double* hsnow,
                       const double* albedo, const double* albedo_snow,
                       int alb_scheme,
                       double alb_smax, double alb_smin, double albi, double albl,
                       const double* Tamp,
                       double tmin, double tmax, double tmid,
                       double mcrit, double wcrit,
                       double tau_a, double tau_f, double afac,
                       bool verbose,
                       double* tsurf_out,  double* smb_out,
                       double* smbi_out,   double* smbs_out,
                       double* saccu_out,  double* smelt_out,
                       double* refr_out,   double* alb_out,
                       double* alb_snow_out,
                       double* hsnow_out,  double* hice_out,
                       double* qmr_out,
                       double* runoff_out, double* subl_out);

#endif /* _SEMIC_H_ */
