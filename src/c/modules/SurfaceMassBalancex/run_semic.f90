subroutine run_semic(sf_in, rf_in, swd_in, lwd_in, wind_in, &
      sp_in, rhoa_in, qq_in, tt_in, tsurf_out, &
      smb_out, saccu_out, smelt_out)

   use utils
   use surface_physics

   implicit none

   ! declare surface physics class
   type(surface_physics_class) :: surface
   ! declare forcing class
   type(forc_class) :: forc
   ! declare validation class
   !type(vali_class) :: vali	! validation not needed here

   integer, parameter:: dp=kind(0.d0)  !< define precision (machine specific)
   integer :: i, k, nx, nloop, ntime, year, day

   ! forcing data    
   double precision, intent(in),dimension(1,365) :: sf_in  ! snow fall rate [m/s]
   double precision, intent(in),dimension(1,365) :: rf_in  ! rain fall rate [m/s]
   double precision, intent(in),dimension(1,365) :: swd_in ! downwelling shortwave radiation [W/m2]
   double precision, intent(in),dimension(1,365) :: lwd_in ! downwelling longwave radiation [W/m2]
   double precision, intent(in),dimension(1,365) :: wind_in! surface wind speed [m/s]
   double precision, intent(in),dimension(1,365) :: sp_in  ! surface pressure [Pa]
   double precision, intent(in),dimension(1,365) :: rhoa_in! air density [kg/m3]
   double precision, intent(in),dimension(1,365) :: qq_in  ! air specific humidity [kg/kg]
   double precision, intent(in),dimension(1,365) :: tt_in  ! air temperature [K]

   ! output data
   double precision :: tsurf_out  ! Ice surface Temperature [K]
   double precision :: smb_out  ! surface mass balance=(Accu-Melt) [m/s]
   double precision :: saccu_out  ! accumulation [m/s]
   double precision :: smelt_out  ! ablation [m/s]

   double precision :: total_time, start, finish

   ! set parameters
   character (len=256) :: name         ! not used(?)
   character (len=256) :: boundary(30) ! not used(?)
   character (len=256) :: alb_scheme   
   integer :: n_ksub    
   double precision    :: ceff
   double precision    :: albi
   double precision    :: albl
   double precision    :: alb_smax
   double precision    :: alb_smin
   double precision    :: hcrit
   double precision    :: rcrit
   double precision    :: amp
   double precision    :: csh
   double precision    :: clh
   double precision    :: tmin
   double precision    :: tmax
   double precision    :: tstic
   double precision    :: tsticsub
   double precision    :: tau_a
   double precision    :: tau_f
   double precision    :: w_crit
   double precision    :: mcrit
   double precision    :: afac
   double precision    :: tmid

   nloop = 10
   nx = 1
   ntime = 365
   year = 0

   ! set vector length
   surface%par%nx = nx

   ! set input (forcing data)
   allocate(forc%sf(nx,ntime))
   allocate(forc%rf(nx,ntime))
   allocate(forc%swd(nx,ntime))
   allocate(forc%lwd(nx,ntime))
   allocate(forc%wind(nx,ntime))
   allocate(forc%sp(nx,ntime))
   allocate(forc%rhoa(nx,ntime))
   allocate(forc%tt(nx,ntime))
   allocate(forc%qq(nx,ntime))

!	write(*,*) sf_in
   forc%sf(:,:) = sf_in
   forc%rf(:,:) = rf_in
   forc%swd(:,:) = swd_in
   forc%lwd(:,:) = lwd_in
   forc%wind(:,:) = wind_in
   forc%sp(:,:) = sp_in
   forc%rhoa(:,:) = rhoa_in
   forc%qq(:,:) = qq_in
   forc%tt(:,:) = tt_in

   ! FIXME should be user input
   !boundary = "" "" ""
   surface%par%tstic = 86400.0_dp
   surface%par%ceff= 2.0e6_dp
   surface%par%csh = 2.0e-3_dp
   surface%par%clh = 5.0e-4_dp
   surface%par%alb_smax = 0.79_dp
   surface%par%alb_smin = 0.6_dp
   surface%par%albi = 0.41_dp
   surface%par%albl = 0.07_dp
   surface%par%tmin = -999_dp
   surface%par%tmax = 273.15_dp
   surface%par%hcrit = 0.028_dp
   surface%par%rcrit = 0.85_dp
   surface%par%amp = 3.0_dp
   surface%par%alb_scheme = "None"
   surface%par%tau_a = 0.008_dp
   surface%par%tau_f = 0.24_dp
   surface%par%w_crit = 15.0_dp
   surface%par%mcrit = 6.0e-8_dp
   surface%par%n_ksub = 3.0_dp
   
   ! initialize sub-daily time step tsticsub
   surface%par%tsticsub = surface%par%tstic / dble(surface%par%n_ksub)

   ! allocate necessary arrays for surface_physics module
   call surface_alloc(surface%now,surface%par%nx)

   ! initialise prognostic variables
   surface%now%mask(:) = 2.0_dp !loi_mask(:nx)
   surface%now%hsnow(:) = 1.0_dp
   surface%now%hice(:)  = 0.0_dp
   surface%now%alb(:) = 0.8_dp
   surface%now%tsurf(:) = 260.0_dp
   surface%now%alb_snow(:) = 0.8_dp
   !surface%now%acc(:) = 0.0_dp
   !surface%now%smb(:) = 0.0_dp
   !surface%now%melt(:) = 0.0_dp
   surface%now%qmr_res(:) = 0.0_dp

   tsurf_out = 0.0_dp
   smb_out = 0.0_dp
   saccu_out = 0.0_dp
   smelt_out = 0.0_dp

   ! define boundary conditions (not used, here!)
   call surface_boundary_define(surface%bnd,surface%par%boundary)
   !call print_boundary_opt(surface%bnd)

   do k=1,nloop ! re-iterate 'nloop' times

   day = 1

   do i=1,ntime ! loop over one year


   ! read input for i-th day of year
   surface%now%sf = forc%sf(1,day)
   surface%now%rf = forc%rf(1,day)
   surface%now%sp = forc%sp(1,day)
   surface%now%lwd = forc%lwd(1,day)
   surface%now%swd = forc%swd(1,day)
   surface%now%wind = forc%wind(1,day)
   surface%now%rhoa = forc%rhoa(1,day)
   surface%now%t2m = forc%tt(1,day)
   surface%now%qq = forc%qq(1,day)

   ! calculate prognostic and diagnsotic variables
   call cpu_time(start)
   call surface_energy_and_mass_balance(surface%now,surface%par,surface%bnd,day,-1)
   call cpu_time(finish)
   total_time = total_time + (finish - start)

   if (k==nloop) then 
      tsurf_out=tsurf_out+surface%now%tsurf(1)*1.0_dp/365.0_dp
      smb_out=smb_out+surface%now%smb(1)*1.0_dp/365.0_dp
      saccu_out=saccu_out+surface%now%alb(1)*1.0_dp/365.0_dp
      smelt_out=smelt_out+surface%now%melt(1)*1.0_dp/365.0_dp
   endif
   day = day + 1

   end do

   end do

   ! de-allocate surface_physics arrays
   call surface_dealloc(surface%now)

   !write(*,*) 'total time for surface_physics:', nloop, total_time

end subroutine run_semic
