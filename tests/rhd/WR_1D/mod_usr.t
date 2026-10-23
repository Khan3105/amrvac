!>  Based on WR_1D from amrvac 2.2.4 - adapted by Kian Tedeskino 2026
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_eos, only: eos
  use mod_fld
  use mod_physics, only: phys_to_primitive, phys_to_conserved
  use mod_multigrid_coupling


  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: StefBoltz

  double precision :: cak_Q, cak_a, cak_base, cak_x0, cak_x1
  integer :: it_start_cak
  double precision :: rho_bound, v_inf, Mdot
  double precision :: R_star, M_star

  integer :: i_v1, i_v2, i_p
  integer :: i_Trad, i_Tgas, i_Mdot, i_Opal, i_CAK, i_lambda, i_fld_R
  integer :: i_Gamma, i_Lum, i_F1, i_F2, i_grE

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out
  logical :: fixed_lum, Cak_in_D, read_cak_table

  logical          :: use_beta_law = .false.
  double precision :: beta_vel = 0.5d0, beta_b = 1.d0, T_out_ic = 5.d4


contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()
    ! Check that the EoS keeps the mean molecular weight constant
    call check_eos_constant_mu()

    ! Choose coordinate system as 2D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_1D")
    
    ! Set independent units BEFORE hd_activate derives the EOS/FLD constants.
    ! eos_init has already read eos_list; params_read reads wind_list in cgs.
    call params_read(par_files)
    unit_length = R_star
    unit_density = rho_bound
    unit_velocity = v_inf

    ! Initialize units
    usr_set_parameters => initglobaldata_usr

    ! A routine for initial conditions is always required
    usr_init_one_grid => initial_conditions

    ! Specify other user routines, for a list see mod_usr_methods.t
    ! Boundary conditions
    usr_special_bc => boundary_conditions
    usr_special_mg_bc => mg_boundary_conditions

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => OPAL_and_CAK

    !> Additional variables
    usr_process_grid => update_extravars

    ! Refine mesh near base
    usr_refine_grid => refine_base

    ! Active the physics module
    call hd_activate()

    i_v1 = var_set_extravar("v1", "v1")
    i_p = var_set_extravar("p","p")
    i_Trad = var_set_extravar("Trad", "Trad")
    i_Tgas = var_set_extravar("Tgas", "Tgas")
    i_Mdot = var_set_extravar("Mdot", "Mdot")
    i_Opal = var_set_extravar("OPAL", "OPAL")
    i_CAK = var_set_extravar("CAK", "CAK")
    i_lambda = var_set_extravar("lambda", "lambda")
    i_fld_R = var_set_extravar("fld_R", "fld_R")
    i_Gamma = var_set_extravar("Gamma", "Gamma")
    i_Lum = var_set_extravar("Lum", "Lum")
    i_F1 = var_set_extravar("F1", "F1")
    i_grE = var_set_extravar("grE", "grE")

  end subroutine usr_init

  !> Stop at start-up unless the EoS keeps the mean molecular weight constant
  subroutine check_eos_constant_mu()
    use mod_global_parameters, only: mype
    if (trim(eos%eos_type) /= 'FI') then
      if (mype == 0) write(*,*) "EoS type is ", trim(eos%eos_type), " but 'FI' is required"
      call mpistop("Equations assume constant mean molecular weight, set eos_type to 'FI'")                  ! short reason
    endif
  end subroutine check_eos_constant_mu


  subroutine initglobaldata_usr
    use mod_global_parameters
    use mod_opal_opacity, only: init_opal_table
    use mod_cak_opacity, only: init_cak_table

    use mod_fld

    integer :: i

    if (SI_unit) call mpistop('WR_1D opacity tables and input parameters use cgs')
    if (.not.hd_energy .or. .not.hd_radiation_fld .or. .not.hd_gravity) &
         call mpistop('WR_1D requires HD energy, FLD and gravity')
    if (ndim /= 1 .or. .not.slab) call mpistop('WR_1D requires Cartesian_1D')
    if (cak_Q <= 0.d0 .or. min(cak_a,cak_base) <= 0.d0 .or. &
         max(cak_a,cak_base) >= 1.d0 .or. cak_x1 <= cak_x0) &
         call mpistop('Invalid CAK parameters')

    !> Initialise Opal
    call init_opal_table(fld_opal_table)

    !> Initialise CAK tables
    call init_cak_table(fld_opal_table)

    R_star = R_star/unit_length
    M_star = M_star/unit_mass
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    Mdot = Mdot*unit_time/unit_mass
    kappa_e = kappa_e/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)
    StefBoltz = const_sigmaSB*unit_temperature**4/unit_radflux

    !> Very bad initial guess for gradE using kappa_e
    gradE = -3.d0*kappa_e*rho_bound*F_bound/c_norm
    E_out = arad_norm*(1.5d4/unit_temperature)**4

    if(mype==0)then
    print*, '==== WR_1D parameters ==============================='
    print*, 'R_star', R_star
    print*, 'M_star (Msun)', M_star*unit_mass/M_sun
    print*, 'L_bound (Lsun)', L_bound*(unit_radflux*unit_length**2)/L_sun
    print*, 'log10(L_bound)', log10(L_bound*(unit_radflux*unit_length**2)/L_sun)
    print*, 'Mdot_ini (Msun/yr)', Mdot*(unit_mass/unit_time)/M_sun*year
    print*, 'kappa_e', kappa_e
    print*, 'F_bound (ergs/cm^2/s)', F_bound*unit_radflux
    print*, 'rho_bound (g/cm^3)', rho_bound*unit_density
    print*, 'v_inf (km/s)', v_inf*unit_velocity/1.d5
    print*, 'Gamma_e_bound', Gamma_e_bound
    print*, '==========================================================='
    endif

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound, kappa_e, &
    R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak, fixed_lum, Cak_in_D, read_cak_table, &
    use_beta_law, beta_vel, beta_b, T_out_ic

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/kappa_e
    F_bound = L_bound/(4*dpi*R_star**2)
    Mdot = Mdot*M_sun/year

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_constants
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer :: ii

    if (use_beta_law) then
      call beta_law_profile(ixI^L, ixO^L, w, x)
    else
      do ii = ixOmin1,ixOmax1
        w(ii,rho_) = read_initial_conditions(x(ii,1),2)
        w(ii,mom(1)) = read_initial_conditions(x(ii,1),3)
        w(ii,e_) = read_initial_conditions(x(ii,1),4)
        w(ii,r_e) = read_initial_conditions(x(ii,1),5)
        w(ii,i_diff_mg) = read_initial_conditions(x(ii,1),6)
      enddo
    endif

    call hd_to_conserved(ixI^L,ixO^L,w,x)

  end subroutine initial_conditions


  function read_initial_conditions(r_in,index) result(var)
    integer, intent(in) :: index
    double precision, intent(in) :: r_in
    double precision :: var

    double precision :: w(1:6), w_mo(1:6), w_po(1:6)
    integer :: ll

    w(:) = 0.d0

    open(unit=1, file='1D_stable.blk')
    read(1,*) !>header
    read(1,*) !>header
    read(1,*) !>header
    read(1,*) w !> first line of data
    do ll = 1,1024
      w_mo = w
      read(1,*) w
        if (w(1) .gt. min(9.d0,r_in)) then
          w_po = w
          goto 8765
        endif
      w_po = w
    enddo

8765 CLOSE(1)
    var = w_mo(index) + (w_po(index) - w_mo(index))/(w_po(1) - w_mo(1))*(r_in - w_mo(1))

  end function read_initial_conditions

  !> beta-law initial condition (Moens+2022 Eqs. 51-54), primitive variables in code units
  subroutine beta_law_profile(ixI^L, ixO^L, w, x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    integer, parameter :: n_int = 2000
    integer          :: ii, k
    double precision :: r, rho, E, T, p, E_out_ic, dr, ra, rb

    E_out_ic = const_rad_a * T_out_ic**4 / unit_pressure   ! a_r T_out^4 / unit_pressure

    do ii = ixOmin1, ixOmax1
      r   = x(ii,1)
      rho = Mdot / (4.d0*dpi*r**2 * beta_vel_of(r))  ! Eq. 51, using beta_vel_of(r)

      E  = E_out_ic                     ! Eq. 53: trapezoid from r to xprobmax1
      dr = (xprobmax1 - r)/n_int
      do k = 1, n_int
        ra = r + (k-1)*dr
        rb = r + k*dr
        E  = E + 0.5d0*(integrand_E(ra)+integrand_E(rb))*dr
      enddo

      T = (E/arad_norm)**0.25d0                           ! from E and arad_norm
      p = rho * T
      w(ii,rho_)   = rho
      w(ii,mom(1)) = beta_vel_of(r)
      w(ii,e_)     = p
      w(ii,r_e)    = E
    enddo
  end subroutine beta_law_profile

  !> beta-law velocity in code units
  double precision function beta_vel_of(r)
    double precision, intent(in) :: r
    beta_vel_of = v_inf * (1 - beta_b*R_star/r)**beta_vel                   ! Eq. 52
  end function beta_vel_of

  !> Energy density integrand helper
  double precision function integrand_E(r)
    double precision, intent(in) :: r
    integrand_E = 3*kappa_e*L_bound*Mdot/(16*dpi**2*c_norm*r**4*beta_vel_of(r))
  end function integrand_E


  subroutine boundary_conditions(qdt,qt,ixI^L,ixB^L,iB,w,x)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)             :: ixI^L, ixB^L, iB
    double precision, intent(in)    :: qdt, qt, x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixI^S), Temp(ixI^S)
    double precision :: Temp0, rho0, T_out, n
    double precision :: Local_gradE(ixI^S), F_adv
    double precision :: Local_tauout(ixB^S)
    double precision :: Local_Tout(ixB^S)

    double precision :: kappa_out, W_geo(ixB^S)

    integer :: ix^D, ixP^L, ixC^L
    double precision :: wprim(ixI^S,1:nw), kappa_D(ixI^S), lambda_D(ixI^S), R_D(ixI^S)
    double precision :: D_base

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,rho_) = dexp(2*dlog(w(ix1+1,rho_)) - dlog(w(ix1+2,rho_)))
      enddo

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,mom(1)) = w(ix1+1,mom(1))
      enddo

      where(w(ixB^S,mom(1)) .lt. 0.d0)
       w(ixB^S,mom(1)) = 0.d0
      endwhere

      call get_kappa_OPAL(ixI^L,ixI^L,w,x,kappa)
      do ix1 = ixBmin1,ixBmax1
        kappa(ix1) = kappa(ixBmax1+1)
      enddo

      F_adv = 4.d0/3.d0*(w(ixBmax1,mom(1))/w(ixBmax1,rho_))*w(ixBmax1,r_e)
                 !> * 4*dpi*xprobmin1**2 this turned F_out into a luminosity
      if (F_adv .ne. F_adv) F_adv = 0.d0

      !> In 3.3 w(:,i_diff_mg) holds D*dt during the implicit step,
      !> so compute D in the first interior cell here
      !ixPmin1 = ixBmax1; ixPmax1 = ixBmax1+2
      !ixCmin1 = ixBmax1+1; ixCmax1 = ixBmax1+1
      !wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
      !call phys_to_primitive(ixI^L,ixP^L,wprim,x)
      !call fld_get_opacity_prim(wprim, x, ixI^L, ixC^L, kappa_D, fld_fl)
      !call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixC^L, lambda_D, R_D, 1, fld_fl)
      !D_base = c_norm*lambda_D(ixCmin1)/(kappa_D(ixCmin1)*wprim(ixCmin1,rho_))
      !> Updated mod_fld.t such that D is not multiplied by qdt
      Local_gradE(ixI^S) = -(F_bound-F_adv)/w(nghostcells+1,i_diff_mg)
      gradE = Local_gradE(nghostcells)
      !Run was using same bad guess for gradE for entire run, update:
      mg%bc(1, mg_iphi)%bc_value = gradE 

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,r_e) = w(ix1+2,r_e) &
        + (x(ix1,1)-x(ix1+2,1))*Local_gradE(ix1+1)
      enddo

      temp(ixB^S) = (w(ixB^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
      w(ixB^S,e_) = w(ixB^S,rho_)*temp(ixB^S)/(eos%gamma-1.d0) + half*w(ixB^S,mom(1))**2/w(ixB^S,rho_)

    case(2)

      !> Compute mean kappa in outer blocks
      call get_kappa_OPAL(ixI^L,ixI^L,w,x,kappa)

      kappa_out = kappa(ixImax1-nghostcells)
      ! kappa_out = kappa_e

      if (kappa_out .ne. kappa_out) kappa_out = kappa_e
      kappa_out = max(kappa_out,kappa_e)
      kappa_out = min(kappa_out,20*kappa_e)

      Local_tauout(ixB^S) = kappa_out*w(ixB^S,rho_)*R_star**2/(3*x(ixB^S,1))
      W_geo(ixB^S) = 0.5d0*(1.d0 - sqrt(1.d0 - (R_star/x(ixB^S,1))**2.d0)) !> Geometric dilution factor
      Local_Tout(ixB^S) = (F_bound/StefBoltz*(W_geo(ixB^S) + 3.d0/4.d0*Local_tauout(ixB^S)))**0.25d0

      T_out = Local_Tout(ixBmin1)

      T_out = max(1.5d4/unit_temperature, T_out)
      ! T_out = max(3.5d4/unit_temperature, T_out)
      E_out = const_rad_a*(T_out*unit_temperature)**4.d0/unit_pressure
      !Same as above, update:
      mg%bc(2, mg_iphi)%bc_value = E_out

      do ix1 = ixBmin1,ixBmax1
        w(ix1,r_e) = 2*w(ix1-1,r_e) - w(ix1-2,r_e)
      enddo

      !> Old Code: Stefan-Boltzmann law for outer BC, but this crashes with new Lucy E_out  
      !> w(ixB^S,r_e) = const_rad_a*(Local_Tout(ixB^S)*unit_temperature)**4.d0/unit_pressure
      !> Use 1/r^2 law of freely streaming radiation 
      do ix1 = ixBmin1,ixBmax1
        w(ix1,r_e) = w(ixBmin1-1,r_e)*(x(ixBmin1-1,1)/x(ix1,1))**2
      enddo

    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling
    
    integer, intent(in) :: iB
    integer :: iigrid,igrid,mpierr
    double precision :: sendbuf(2),recvbuf(2),value
    logical :: owns_boundary

    owns_boundary=.false.
    if (mg%is_allocated) then
       do iigrid=1,igridstail
          igrid=igrids(iigrid)
          if (ps(igrid)%is_physical_boundary(iB)) owns_boundary=.true.
       enddo
    endif
    sendbuf=0.d0
    if (owns_boundary) then
       sendbuf(2)=1.d0
       if (iB == 1) sendbuf(1)=gradE
       if (iB == 2) sendbuf(1)=E_out
    endif
    call MPI_ALLREDUCE(sendbuf,recvbuf,2,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,mpierr)
    if (recvbuf(2)>0.d0) then
       if (nint(recvbuf(2))/=1) call mpistop('WR_1D: expected one owner per 1D boundary')
       value=recvbuf(1)
    else
       if (iB == 1) value=gradE
       if (iB == 2) value=E_out
    endif
    if (.not.ieee_is_finite(value)) call mpistop('WR_1D: non-finite MG boundary value')
    select case(iB)
    case(1)
       mg%bc(iB,mg_iphi)%bc_type=mg_bc_neumann
    case(2)
       mg%bc(iB,mg_iphi)%bc_type=mg_bc_dirichlet
    case default
       call mpistop('WR_1D: unknown MG boundary')
    end select
    mg%bc(iB,mg_iphi)%bc_value=value
  end subroutine mg_boundary_conditions


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,ndim)

    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixI^S) = x(ixI^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,1) = -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length)

  end subroutine set_gravitation_field

  !> Calculate w(iw)=w(iw)+qdt*SOURCE[wCT,qtC,x] within ixO for all indices
  !> iw=iwmin...iwmax.  wCT is at time qCT
  subroutine PseudoPlanar(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in)    :: qdt, qtC, qt
    double precision, intent(in)    :: wCT(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: ppsource(ixO^S,1:nw)

    double precision :: k_cak(ixO^S), rad_flux(ixI^S,1:ndim)

    call PseudoPlanarSource(ixI^L,ixO^L,wCT,x,ppsource)
    w(ixO^S,rho_) = w(ixO^S,rho_) + qdt*ppsource(ixO^S,rho_) !> OK
    w(ixO^S,mom(1)) = w(ixO^S,mom(1)) + qdt*ppsource(ixO^S,mom(1)) !> OK
    w(ixO^S,e_) = w(ixO^S,e_) + qdt*ppsource(ixO^S,e_) !> OK
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      call get_kappa_CAK(ixI^L,ixO^L,wCT,x,k_cak)

      if (fixed_lum) then
        !> Fixed L = L_bound
        w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
          + qdt*wCT(ixO^S,rho_)*L_bound/(4*dpi*x(ixO^S,1)**2)/const_c*k_cak(ixO^S)*unit_velocity
        w(ixO^S,e_) = w(ixO^S,e_) &
            + qdt*wCT(ixO^S,mom(1))*L_bound/(4*dpi*x(ixO^S,1)**2)/const_c*k_cak(ixO^S)*unit_velocity
      else
        !> Local flux
        call fld_get_radflux(wCT, x, ixI^L, ixO^L, rad_flux, fld_fl)
        w(ixO^S,mom(1)) = w(ixO^S,mom(1)) &
          + qdt*wCT(ixO^S,rho_)*rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity
        w(ixO^S,e_) = w(ixO^S,e_) &
            + qdt*wCT(ixO^S,mom(1))*rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity
      endif
    endif

  end subroutine PseudoPlanar

  subroutine get_dt_cak(w,ixI^L,ixO^L,dtnew,dx^D,x)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: dx^D, x(ixI^S,1:ndim)
    double precision, intent(in)    :: w(ixI^S,1:nw)
    double precision, intent(inout) :: dtnew

    double precision :: radius(ixI^S)
    double precision :: mass

    double precision :: dt_cak
    double precision :: k_cak(ixO^S), rad_flux(ixI^S,1:ndim)

    call get_kappa_CAK(ixI^L,ixO^L,w,x,k_cak)

    if (fixed_lum) then
      !> Fixed L = L_bound
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/(L_bound/(4*dpi*x(ixO^S,1)**2)/const_c*k_cak(ixO^S)*unit_velocity &
      -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    else
      !> Local flux
      call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, fld_fl)
      dt_cak = courantpar*minval(dsqrt(dxlevel(1)/abs(rad_flux(ixO^S,1)/const_c*k_cak(ixO^S)*unit_velocity &
      -const_G*mass/radius(ixI^S)**2*(unit_time**2/unit_length))))
    endif

    dtnew = min(dt_cak, dtnew)

  end subroutine get_dt_cak

  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)

    double precision :: rad_flux(ixI^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,1:ndim)
    double precision :: radius(ixO^S),  pert(ixO^S)

    double precision :: edd(ixI^S,1:ndim,1:ndim), lambda(ixI^S), fld_R(ixI^S)
    double precision :: wprim(ixI^S,1:nw)

    integer :: rdir

    source(ixO^S,1:nw) = zero

    rdir = 1

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S)

    !> de/dt = -2 (e+p) v_r/r
     source(ixO^S,e_) = -two*(w(ixO^S,e_)+pth(ixO^S))*v(ixO^S,rdir)/radius(ixO^S)

    !> dEr/dt = -2 (E v_r + F_r)/r
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, fld_fl)
    source(ixO^S,r_e) = source(ixO^S,r_e) - two*rad_flux(ixO^S,rdir)/radius(ixO^S)

    source(ixO^S,r_e) = source(ixO^S,r_e) - two*w(ixO^S,r_e)*v(ixO^S,rdir)/radius(ixO^S)

    wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_eddington(wprim, x, ixI^L, ixO^L, edd, lambda, fld_R, nghostcells, fld_fl)
    ! Missing transverse radiation work in spherical geometry, checked this analytically:
    ! -v/r*(P_theta_theta+P_phi_phi) = -(1-f_rr)*E*v/r.
    source(ixO^S,r_e) = source(ixO^S,r_e) &
         -(one-edd(ixO^S,1,1))*v(ixO^S,rdir)*w(ixO^S,r_e)/radius(ixO^S)
  end subroutine PseudoPlanarSource


  subroutine OPAL_and_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixI^S)
    double precision :: wcons(ixI^S,1:nw)
    
    double precision :: OPAL(ixO^S), CAK(ixO^S)

    !> Conserved from primitive vars, needed for get_kappa_OPAL
    wcons(ixI^S,1:nw) = w(ixI^S,1:nw)
    call phys_to_conserved(ixI^L,ixI^L,wcons,x)

    !> Get OPAL opacities by reading from table
    call get_kappa_OPAL(ixI^L,ixO^L,wcons,x,OPAL)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    if (Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixI^L,ixO^L,wcons,x,CAK)
      else
        call get_kappa_CAK(ixI^L,ixO^L,wcons,x,CAK)
      endif
    else
      CAK(ixO^S) = 0.d0
    endif

    !> Add OPAL and CAK for total opacity
    kappa(ixO^S) = OPAL(ixO^S) + CAK(ixO^S)

  end subroutine OPAL_and_CAK


  subroutine get_kappa_OPAL(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    integer :: ix^D
    double precision :: Temp(ixI^S)
    double precision :: kappaval, rho0, Temp0

    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Temp)

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,kappaval)
        kappa(ix^D) = kappaval/unit_opacity
    {enddo\ }

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)
    double precision :: xx(ixO^S), alpha(ixO^S)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)

    call gradient(vel,ixI^L,ixO^L,1,gradvI)
    gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))

    xx(ixO^S) = 1.d0-xprobmin1/x(ixO^S,1)

    alpha(ixO^S) = cak_a

    where (xx(ixO^S) .le. cak_x0)
      alpha(ixO^S) = cak_base
    elsewhere (xx(ixO^S) .le. cak_x1)
      alpha(ixO^S) = cak_base + (cak_a - cak_base)&
      *(xx(ixO^S) - cak_x0)/(cak_x1 - cak_x0)
    endwhere

    kappa(ixO^S) = kappa_e*cak_Q/(1-alpha(ixO^S)) &
    *(gradv(ixO^S)*unit_velocity/(w(ixO^S,rho_)*const_c*cak_Q*kappa_e))**alpha(ixO^S)

    if (x(ixImax1,1) .ge. xprobmax1) then
      kappa(ixOmax1) = kappa(ixOmax1-1)
    endif

  end subroutine get_kappa_CAK

  subroutine get_kappa_CAK2(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_cak_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    double precision :: Temp(ixI^S), rho0, temp0, gradv0, kap0
    integer :: ix^D

    double precision :: alpha, Qbar, Q0, kappa_e_t
    double precision :: tau, M_t
    double precision :: vel(ixI^S), gradv(ixO^S), gradvI(ixI^S)

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)

    call gradient(vel,ixI^L,ixO^L,1,gradvI)
    gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))

    !> Get CAK opacities by reading from table
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Temp)

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix^D)*(unit_velocity/unit_length)
        call set_cak_opacity(rho0,Temp0,alpha,Qbar,Q0,kappa_e_t)

        tau = (kappa_e*unit_opacity)*rho0*const_c/gradv0
        M_t = Qbar/(1-alpha)*((1+Q0*tau)**(1-alpha) - 1)/(Q0*tau)
        kap0 = (kappa_e*unit_opacity)*M_t

        kappa(ix^D) = kap0/unit_opacity

        kappa(ix^D) = min(50*kappa_e,kappa(ix^D))
    {enddo\ }


    if (x(ixImax1,1) .ge. xprobmax1) then
      kappa(ixOmax1) = kappa(ixOmax1-1)
    endif

  end subroutine get_kappa_CAK2


  subroutine refine_base(igrid,level,ixG^L,ix^L,qt,w,x,refine,coarsen)
    ! Enforce additional refinement or coarsening
    ! One can use the coordinate info in x and/or time qt=t_n and w(t_n) values w.
    ! you must set consistent values for integers refine/coarsen:
    ! refine = -1 enforce to not refine
    ! refine =  0 doesn't enforce anything
    ! refine =  1 enforce refinement
    ! coarsen = -1 enforce to not coarsen
    ! coarsen =  0 doesn't enforce anything
    ! coarsen =  1 enforce coarsen
    use mod_global_parameters

    integer, intent(in) :: igrid, level, ixG^L, ix^L
    double precision, intent(in) :: qt, w(ixG^S,1:nw), x(ixG^S,1:ndim)
    integer, intent(inout) :: refine, coarsen

    refine=0
    coarsen=0

    !> Refine close to base
    refine = -1

    if (qt .gt. 1.d0) then
      if (any(x(ixG^S,1) < 1.d0)) refine=1
    endif

    if (qt .gt. 2.d0) then
      if (any(x(ixG^S,1) < 1.d0)) refine=1
    endif

    if (qt .gt. 4.d0) then
      if (any(x(ixG^S,1) < 1.d0)) refine=1
    endif

  end subroutine refine_base


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
    double precision                   :: g_grav(ixI^S)
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: kappa(ixI^S), OPAL(ixO^S), CAK(ixO^S)
    double precision                   :: vel(ixI^S), gradv(ixI^S), gradE(ixI^S)
    double precision                   :: rad_flux(ixI^S,1:ndim), Lum(ixO^S)
    double precision                   :: pth(ixI^S), lambda(ixI^S), fld_R(ixI^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass
    double precision :: wprim(ixI^S,1:nw)

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_opacity_prim(wprim, x, ixI^L, ixO^L, kappa, fld_fl)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, fld_fl)

    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Tgas)
    call hd_get_trad(w, x, ixI^L, ixO^L, Trad)

    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)
    call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)

    g_rad(ixO^S) = (OPAL(ixO^S)+CAK(ixO^S))*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv)
    call gradient(w(ixI^S,r_e),ixI^L,ixO^L,1,gradE)

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R, nghostcells,fld_fl)

    Lum(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun

    w(ixO^S,i_v1) = w(ixO^S,mom(1))/w(ixO^S,rho_)
    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    w(ixO^S,i_p) = pth(ixO^S)

    w(ixO^S,i_Trad) = Trad(ixO^S)*unit_temperature
    w(ixO^S,i_Tgas) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,i_Mdot) = 4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
    *unit_density*unit_velocity/M_sun*year
    w(ixO^S,i_Opal) = OPAL(ixO^S)/kappa_e
    w(ixO^S,i_CAK) = CAK(ixO^S)/kappa_e
    w(ixO^S,i_lambda) = lambda(ixO^S)
    w(ixO^S,i_fld_R) = fld_R(ixO^S)
    w(ixO^S,i_Gamma) = big_gamma(ixO^S)
    w(ixO^S,i_Lum) = Lum(ixO^S)
    w(ixO^S,i_F1) = rad_flux(ixO^S,1)/F_bound
    w(ixO^S,i_grE) = gradE(ixO^S)

  end subroutine update_extravars
end module mod_usr
