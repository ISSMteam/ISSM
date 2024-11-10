subroutine run_semic_transient(nx, ntime, nloop, sf_in, rf_in, swd_in, lwd_in, wind_in, &
      sp_in, rhoa_in, qq_in, tt_in, tsurf_in, qmr_in, &
      tstic, &
      hcrit, rcrit, &
      mask, hice, hsnow, &
      albedo, albedo_snow, &
      alb_scheme, alb_smax, alb_smin, albi, albl, &
      Tamp, &
      tmin, tmax, tmid, mcrit, w_crit, tau_a, tau_f, afac, verbose, &
      tsurf_out, smb_out, smbi_out, smbs_out, saccu_out, smelt_out,  refr_out, alb_out, & 
      alb_snow_out,hsnow_out,hice_out,qmr_out, runoff_out, subl_out) !{{{

   use utils
   use surface_physics
   implicit none

   ! declare surface physics class
   type(surface_physics_class) :: surface
   ! declare forcing class
   !type(forc_class) :: forc
   ! declare validation class
   !type(vali_class) :: vali	! validation not needed here

   integer, parameter:: dp=kind(0.d0)  !< define precision (machine specific)
   integer :: i, k, n
   integer :: nnx, nny
   integer :: year=0
   integer :: day =1 !< not used value.
   integer, intent(in) :: nx, ntime, nloop   ! number of grid / number of times
   logical, intent(in) :: verbose            ! verbosity
   logical :: debug=.false.

   ! forcing data    
   ! input argument format array size (nx, ntime)...
   double precision, intent(in), dimension(nx):: sf_in    ! snow fall rate [m/s]
   double precision, intent(in), dimension(nx):: rf_in    ! rain fall rate [m/s]
   double precision, intent(in), dimension(nx):: swd_in   ! downwelling shortwave radiation [W/m2]
   double precision, intent(in), dimension(nx):: lwd_in   ! downwelling longwave radiation [W/m2]
   double precision, intent(in), dimension(nx):: wind_in  ! surface wind speed [m/s]
   double precision, intent(in), dimension(nx):: sp_in    ! surface pressure [Pa]
   double precision, intent(in), dimension(nx):: rhoa_in  ! air density [kg/m3]
   double precision, intent(in), dimension(nx):: qq_in    ! air specific humidity [kg/kg]
   double precision, intent(in), dimension(nx):: tt_in    ! air temperature [K]

   ! input data
   double precision, intent(in) :: tstic  ! time step from ISSM [sec].

   ! output data
   ! Ice surface Temperature [K]
   double precision, intent(out), dimension(nx):: tsurf_out    
   ! surface mass balance=(Accu-Melt) [m/s]
   double precision, intent(out), dimension(nx):: smb_out     
   double precision, intent(out), dimension(nx):: smbi_out     ! SMB ice  [water equivalent m/s]
   double precision, intent(out), dimension(nx):: smbs_out     ! SMB snow [water equivalent m/s]
   double precision, intent(out), dimension(nx):: saccu_out    ! accumulation [m/s]
   double precision, intent(out), dimension(nx):: smelt_out    ! ablation [m/s]
   double precision, intent(out), dimension(nx):: refr_out     ! freezing [m/s]
   double precision, intent(out), dimension(nx):: alb_out      ! grid-averaged albedo [no unit] 
   double precision, intent(out), dimension(nx):: alb_snow_out 
   double precision, intent(out), dimension(nx):: hice_out    
   double precision, intent(out), dimension(nx):: hsnow_out   
   double precision, intent(out), dimension(nx):: qmr_out     
   double precision, intent(out), dimension(nx):: runoff_out
   double precision, intent(out), dimension(nx):: subl_out

   double precision :: total_time, start, finish

   ! set parameters
   !character (len=256) :: name         ! not used(?)
   !character (len=256) :: boundary(30) ! not used(?)
   !character (len=256), intent(in) :: alb_scheme  !< name of albedo scheme
   integer, intent(in)          :: alb_scheme
   !integer :: n_ksub    
   !double precision             :: ceff         !< surface heat heat capacity of snow/ice [J/W m2]
   double precision, intent(in), dimension(nx):: albedo
   double precision, intent(in), dimension(nx):: albedo_snow !< spatial..
   double precision, intent(in), dimension(nx):: hsnow
   double precision, intent(in), dimension(nx):: hice
   double precision, intent(in), dimension(nx):: tsurf_in    !< input temperature [K]
   double precision, intent(in), dimension(nx):: qmr_in 
   double precision, intent(in), dimension(nx):: mask

   double precision, intent(in) :: albi
   double precision, intent(in) :: albl
   double precision, intent(in) :: alb_smax
   double precision, intent(in) :: alb_smin
   double precision, intent(in) :: hcrit !< critical snow height for which grid cell is 50% snow covered [m]
   double precision, intent(in) :: rcrit !< critical snow height fro which refreezing fraction is 50% [m]
   double precision, intent(in) :: Tamp
   !double precision    :: csh
   !double precision    :: clh
   double precision, intent(in) :: tmin
   double precision, intent(in) :: tmax
   !double precision    :: tsticsub
   ! parameters for isba albedo scheme.
   double precision, intent(in) :: tau_a  !< critical liquide water concent for "isba" albedo scheme [kg/m2]
   double precision, intent(in) :: tau_f
   double precision, intent(in) :: w_crit
   double precision, intent(in) :: mcrit
   double precision, intent(in) :: afac !< param
   double precision, intent(in) :: tmid !< param for "alex" albedo parameterization [K]

   if (debug) then
      print*,'   ntime: ', ntime
      print*,'   nx   : ', nx
   end if

   ! set vector length
   surface%par%nx = nx

   ! FIXME should be user input
   !boundary = "" "" ""
   if (debug) then
      print*, "run_semic_transient: initialize parameters."
   end if
   surface%par%tstic = tstic      !< time step [s]
   surface%par%ceff= 2.0e6_dp     !< surface heat capacity of snow/ice [J/K/m2]
   surface%par%csh = 2.0e-3_dp    !< turbulent heat exchange coefficient 
   surface%par%clh = 5.0e-4_dp    !< turbulent heat exchange coefficient [no unit]
   surface%par%alb_smax = alb_smax !0.79_dp !< max snow albedo
   surface%par%alb_smin = alb_smin !0.6_dp  !< min snow albedo
   surface%par%albi = albi ! 0.41_dp     !< albedo for ice
   surface%par%albl = albl ! 0.07_dp     !< albedo for land
   surface%par%tmin = tmin ! -999_dp
   surface%par%tmax = tmax ! 273.15_dp
   surface%par%hcrit = hcrit !0.028_dp   !< critical snow height for which grid cell is 50 % snow colvered [m]
   surface%par%rcrit = rcrit !0.85_dp    !< refreezing fraction is 50% [m]
   surface%par%amp   = Tamp !3.0_dp   !< amplitude of diurnal cycle [K]
   if (alb_scheme == 0) then
      surface%par%alb_scheme="none"
   else if (alb_scheme == 1) then
      surface%par%alb_scheme = "slater"
   else if (alb_scheme == 2) then
      surface%par%alb_scheme = "denby"
   else if (alb_scheme == 3) then
      surface%par%alb_scheme = "isba"
   else
      print*, "ERROR: current albedo scheme is not available."
      call exit(1)
   end if 
   surface%par%tau_a  = tau_a  !0.008_dp
   surface%par%tau_f  = tau_f  !0.24_dp
   surface%par%w_crit = w_crit !15.0_dp ! default value
   surface%par%mcrit  = mcrit  !6.0e-8_dp
   surface%par%n_ksub = 3      ! sub ...
   ! snow albedo of alex
   surface%par%afac   = afac
   surface%par%tmid   = tmid
   
   ! initialize sub-daily time step tsticsub
   surface%par%tsticsub = surface%par%tstic / dble(surface%par%n_ksub)

   ! allocate necessary arrays for surface_physics module
   call surface_alloc(surface%now,surface%par%nx)

   ! initialise prognostic variables
   if (debug) then
      print*,"run_semic_transient: initialize variables."
   end if
   ! these values will be updated through "surface_energy_and_mass_balance" function.
   surface%now%mask    (:) = mask       (:) ! 2.0_dp  !loi_mask(:nx)
   if (debug) then
      print*,"run_semic_transient: initialize variables: mask"
   end if
   surface%now%hsnow   (:) = hsnow      (:) ! initial snow height...
   surface%now%hice    (:) = hice       (:) ! initial ice height..
   surface%now%tsurf   (:) = tsurf_in   (:) !< initial ice surface temperature
   surface%now%alb     (:) = albedo     (:) !< initial albedo for energy balance.
   surface%now%alb_snow(:) = albedo_snow(:) !< initial albedo for ISBA albedo method.
   if (debug) then
      print*,"run_semic_transient: initialize variables. DONE."
   end if

   if (debug) then
      !print*, "====== global variable =========="
      !print*, "nloop          :", nloop
      !print*, "nx             :", surface%par%nx
      !print*, "======  parameters ======"
      !print*, "csh            :", surface%par%csh
      !print*, "clh            :", surface%par%clh
      !print*, "albeo scheme   :", surface%par%alb_scheme
      !print*, "albeo ice      :", surface%par%albi
      !print*, "tstic          :", surface%par%tstic
      !print*, "tsticsub       :", surface%par%tsticsub
      !print*, "n_ksub         :", surface%par%n_ksub
      print*, "====== inputs ========="
      print*, "hsnow          :", hsnow
      print*, "======  state variables ======"
      print*, "hsnow          :", surface%now%hsnow
      print*, "hice           :", surface%now%hice
      print*, "albeo          :", surface%now%alb
      print*, "albeo snow     :", surface%now%alb_snow
      print*, "mask           :", surface%now%mask
      print*, "tsurf          :", surface%now%tsurf
      print*, "sf             :", sf_in
   end if

   ! define boundary conditions (not used, here!)
   call surface_boundary_define(surface%bnd,surface%par%boundary)
   !call print_boundary_opt(surface%bnd)

   ! input with single value
   do k =1,nloop
      do i =1,ntime
         if (debug) then
            print*,"run_semic_transient: forcing data: ntime = ", i
         end if
         surface%now%sf   = sf_in  !(:,i)
         surface%now%rf   = rf_in  !(:,i)
         surface%now%sp   = sp_in  !(:,i)
         surface%now%lwd  = lwd_in !(:,i)
         surface%now%swd  = swd_in !(:,i)
         surface%now%wind = wind_in!(:,i)
         surface%now%rhoa = rhoa_in!(:,i)
         surface%now%t2m  = tt_in  !(:,i)
         surface%now%qq   = qq_in  !(:,i)
         ! qmr_res is used to "energy_balance" in semic.
         surface%now%qmr_res = qmr_in

         ! calculate prognostic and diagnsotic variables
         call surface_energy_and_mass_balance(surface%now,surface%par,surface%bnd,i,year)
         
         if (debug) then
            print*,"done..."
         end if
         if (k == nloop) then
            tsurf_out         = surface%now%tsurf
            ! melt - potential surface melt [m/s]
            ! smb = SMB_ice + SMB_snow
            ! smbi  - SMB_ice  (water equivalent m/sec)
            ! smbs  - SMB_snow (water equivalent m/sec)
            smb_out           =surface%now%smb      ! smb = smb_snow + smb_ice
            smbi_out          =surface%now%smb_ice  ! Csi (snow>ice) - melted_ice + refrezon_rain.
            smbs_out          =surface%now%smb_snow ! smb_snow = snowfall - sublimiation - melted_snow + refrozen_snow
            saccu_out         =surface%now%acc      ! acc      = snowfall - sublimiation - refreezing 
            smelt_out         =surface%now%melt     ! potential surface melt = melt_ice + melt_snow
            refr_out          =surface%now%refr     ! refreezing values. [m/sec]
            alb_out           =surface%now%alb
            alb_snow_out      =surface%now%alb_snow
            hsnow_out         =surface%now%hsnow
            hice_out          =surface%now%hice
            qmr_out           =surface%now%qmr_res
            runoff_out        =surface%now%runoff   ! runoff values. [m/sec]
            subl_out          =surface%now%subl     ! subl values. [m/sec]
         end if
      end do
   end do

   ! de-allocate surface_physics arrays
   call surface_dealloc(surface%now)

end subroutine run_semic_transient ! }}}
