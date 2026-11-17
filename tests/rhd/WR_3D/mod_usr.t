!> This is a template for a new user problem
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_eos, only: eos
  use mod_fld
  use mod_physics, only: phys_to_primitive, phys_to_conserved
  use mod_multigrid_coupling

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

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

  double precision :: kappa_e, L_bound, Gamma_e_bound, F_bound, gradE, E_out
  logical :: fixed_lum, Cak_in_D
  logical :: read_cak_table = .false.
  character(len=64) :: cak_table = ''
  logical :: CAK_zero = .false.

  integer :: i_v1, i_v2, i_v3, i_p
  integer :: i_Trad, i_Tgas, i_Mdot, i_Opal, i_CAK, i_CAK2, i_lambda, i_edd
  integer :: i_Gamma, i_Lum, i_F1, i_F2, i_F3, i_gradE

  logical          :: use_beta_law = .false.
  double precision :: beta_vel = 0.5d0, beta_b = 1.d0, T_out_ic = 5.d4

  double precision :: sum_time
  double precision, allocatable :: vr_sumt(:), rho_sumt(:), rho2_sumt(:), vr2_sumt(:), sumt(:)
  double precision, allocatable :: rhovr_sumt(:), rho2vr_sumt(:)

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Check that the EoS keeps the mean molecular weight constant
    call check_eos_constant_mu()

    ! Choose coordinate system as 3D Cartesian with three components for vectors
    call set_coordinate_system("Cartesian_3D")

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

    !> lasy fix for inexplicable pressure
    usr_internal_bc => Fix_pressure

    ! PseudoPlanar correction
    usr_source => PseudoPlanar

    ! Graviatational field
    usr_gravity => set_gravitation_field

    ! Special Opacity
    usr_special_opacity => OPAL_and_CAK

    !> Set maximum value for diffusion coefficient
    usr_special_diffcoef => ceil_diffcoef

    ! Refine mesh near base
    usr_refine_grid => refine_base

    ! Write out energy levels and temperature
    usr_write_analysis => collapse_to_1D

    !> Additional variables
    usr_process_grid => update_extravars

    ! Output routines
    ! usr_aux_output    => specialvar_output
    ! usr_add_aux_names => specialvarnames_output

    ! Timestep for PseudoPlanar
    ! usr_get_dt => get_dt_cak

    ! Activate the physics module
    call hd_activate()

    i_v1 = var_set_extravar("v1", "v1")
    i_v2 = var_set_extravar("v2", "v2")
    i_v3 = var_set_extravar("v3", "v3")
    i_p = var_set_extravar("p","p")
    i_Trad = var_set_extravar("Trad", "Trad")
    i_Tgas = var_set_extravar("Tgas", "Tgas")
    i_Mdot = var_set_extravar("Mdot", "Mdot")
    i_Opal = var_set_extravar("OPAL", "OPAL")
    i_CAK = var_set_extravar("CAK", "CAK")
    i_CAK2 = var_set_extravar("CAK2", "CAK2")
    i_lambda = var_set_extravar("lambda", "lambda")
    i_Edd = var_set_extravar("Edd", "Edd")
    i_Gamma = var_set_extravar("Gamma", "Gamma")
    i_Lum = var_set_extravar("Lum", "Lum")
    i_F1 = var_set_extravar("F1", "F1")
    i_F2 = var_set_extravar("F2", "F2")
    i_F3 = var_set_extravar("F3", "F3")
    i_gradE = var_set_extravar("gradE", "gradE")

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

    !> Initialise Opal
    call init_opal_table(fld_opal_table)

    !> Initialise CAK tables
    if (len_trim(cak_table) == 0) cak_table = fld_opal_table
    call init_cak_table(trim(cak_table))

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
    print*, '==== WR_3D parameters ==============================='
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

    sum_time = 0.d0

    allocate(vr_sumt(domain_nx1))
    vr_sumt = 0.d0
    allocate(rho_sumt(domain_nx1))
    rho_sumt = 0.d0
    allocate(rho2_sumt(domain_nx1))
    rho2_sumt = 0.d0
    allocate(vr2_sumt(domain_nx1))
    vr2_sumt = 0.d0
    allocate(rhovr_sumt(domain_nx1))
    rhovr_sumt = 0.d0
    allocate(rho2vr_sumt(domain_nx1))
    rho2vr_sumt = 0.d0

    !> Some checks
    if (SI_unit) call mpistop('WR_3D opacity tables and input parameters use cgs')
    if (.not.hd_energy .or. .not.hd_radiation_fld .or. .not.hd_gravity) &
          call mpistop('WR_3D requires HD energy, FLD and gravity')
    if (ndim /= 3 .or. .not.slab) call mpistop('WR_3D requires Cartesian_3D')

  end subroutine initglobaldata_usr

  !> Read parameters from a file
  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound, kappa_e, &
    R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak, fixed_lum, Cak_in_D, &
    read_cak_table, CAK_zero, use_beta_law, beta_vel, beta_b, T_out_ic, cak_table

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

    ! print*, Gamma_e_bound, dpi, const_G, M_star, const_c

  end subroutine params_read

  !> A routine for specifying initial conditions
  subroutine initial_conditions(ixI^L, ixO^L, w, x)
    use mod_constants
    use mod_physics

    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision :: kappa(ixO^S), lambda(ixO^S), fld_R(ixO^S)
    double precision :: vel(ixI^S)
    double precision :: T_out, E_out, E_gauge
    double precision :: T_in, E_in, rr(ixI^S), bb, temp(ixI^S)

    integer :: ii, ix^D

!    w(ixI^S,rho_)   = rho_bound
!    w(ixI^S,mom(:)) = 0.d0
!
!    where(x(ixI^S,1) .gt. 1.d0)
!      vel(ixI^S) =  v_inf*(1 - 0.999d0*( 1.d0/x(ixI^S,1)))**0.5d0
!      w(ixI^S,rho_) = Mdot/(4.d0*dpi*vel(ixI^S)*x(ixI^S,1)**2)
!    endwhere
!
!    w(ixI^S,mom(1)) = vel(ixI^S)*w(ixI^S,rho_)
!
!    !> Outer/Inner temperature
!    T_out = 28445.836732569689/unit_temperature
!    E_out = const_rad_a*(T_out*unit_temperature)**4/unit_pressure
!    T_in = T_bound
!    E_in = const_rad_a*(T_in*unit_temperature)**4/unit_pressure
!
!    !>Very bad initial profile using constant gradE
!    ! w(ixI^S,r_e) = E_out + gradE*(x(ixI^S,1)-xprobmax1)
!    ! w(ixI^S,r_e) = E_in + gradE*(x(ixI^S,1)-xprobmin1)
!    rr(ixI^S) = dsqrt(x(ixI^S,1)-xprobmin1)*(16*x(ixI^S,1)**2 + 8*x(ixI^S,1)+ 6) &
!                /(15*x(ixI^S,1)**(5.d0/2))
!    bb = -kappa_e*L_bound*Mdot*unit_velocity*3/(16*dpi**2*const_c*v_inf)
!
!    ! bb = bb*2
!
!    ! rr(ixI^S) = 2*dsqrt(x(ixI^S,1)-xprobmin1)/dsqrt(x(ixI^S,1))
!    ! bb = kappa_e*F_bound*Mdot*unit_velocity*3/(4*dpi*const_c*v_inf)
!    w(ixI^S,r_e) = E_in + bb*rr(ixI^S)
!
!    E_gauge = E_in + bb*dsqrt(xprobmax1-xprobmin1)&
!              *(16*xprobmax1**2 + 8*xprobmax1+ 6)/(15*xprobmax1**(5.d0/2))
!
!    w(ixI^S,r_e) = w(ixI^S,r_e) - E_gauge + E_out
!
!      temp(ixI^S) = (w(ixI^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
!      w(ixI^S,e_) = w(ixI^S,rho_)*temp(ixI^S)/(eos%gamma-1.d0) + half*w(ixI^S,mom(1))**2/w(ixI^S,rho_)
!
!
!    call fld_get_opacity(w, x, ixI^L, ixO^L, kappa)
!    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R, nth_for_gradient,fld_fl)
!
!    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)*lambda(ixO^S)/(kappa(ixO^S)*w(ixO^S,rho_))
!    w(ixO^S,i_diff_mg) = (const_c/unit_velocity)/(3.d0*kappa(ixO^S)*w(ixO^S,rho_))
!
    if (use_beta_law) then
      call beta_law_profile(ixI^L, ixO^L, w, x)
    else
      do ii = ixOmin1,ixOmax1
        w(ii,:,:,rho_) = read_initial_conditions(x(ii,3,3,1),2)
        w(ii,:,:,mom(1)) = read_initial_conditions(x(ii,3,3,1),3)
        w(ii,:,:,e_) = read_initial_conditions(x(ii,3,3,1),4)
        w(ii,:,:,r_e) = read_initial_conditions(x(ii,3,3,1),5)
        w(ii,:,:,i_diff_mg) = read_initial_conditions(x(ii,3,3,1),6)
      enddo
    endif

    !> Adding sinusoidal perturbation in y and z
    do ix2 = ixOmin2,ixOmax2
      do ix3 = ixOmin3,ixOmax3
        w(ixOmin1:ixOmax1,ix2,ix3,mom(2)) = 1.d-1*dsin(x(ixOmin1:ixOmax1,ix2,ix3,1)-1.d0) &
             *dcos(6*x(ixOmin1:ixOmax1,ix2,ix3,2))*w(ixOmin1:ixOmax1,ix2,ix3,rho_)
        w(ixOmin1:ixOmax1,ix2,ix3,mom(3)) = 1.d-1*dsin(x(ixOmin1:ixOmax1,ix2,ix3,1)-1.d0) &
             *dcos(6*x(ixOmin1:ixOmax1,ix2,ix3,3))*w(ixOmin1:ixOmax1,ix2,ix3,rho_)
      enddo
    enddo

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
    read(1,*) !> header
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
      r   = x(ii,ixOmin2,ixOmin3,1)
      rho = Mdot / (4.d0*dpi*r**2 * beta_vel_of(r))  ! Eq. 51, using beta_vel_of(r)

      E  = E_out_ic                     ! Eq. 53: trapezoid from r to xprobmax1
      dr = (xprobmax1 - r)/n_int
      do k = 1, n_int
        ra = r + (k-1)*dr
        rb = r + k*dr
        E  = E + 0.5d0*(integrand_E(ra)+integrand_E(rb))*dr
      enddo

      T = (E/arad_norm)**0.25d0
      p = rho * T

      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,rho_)   = rho
      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(1)) = beta_vel_of(r)
      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(2)) = 0.d0
      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,mom(3)) = 0.d0
      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,e_)     = p
      w(ii,ixOmin2:ixOmax2,ixOmin3:ixOmax3,r_e)    = E
      
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

    double precision :: Temp(ixI^S)
    double precision :: Local_gradE(ixI^S)

    integer :: ix^D

    select case (iB)

    case(1)
      w(ixB^S,rho_) = rho_bound
      do ix1 = ixBmax1-1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_) = dexp(2*dlog(w(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_)) - dlog(w(ix1+2,ixBmin2:ixBmax2,ixBmin3:ixBmax3,rho_)))
      enddo

      !w(ixB^S,mom(2)) = 0.d0

      do ix1 = ixBmax1,ixBmin1,-1
        w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1)) = w(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(1))
        w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2)) = w(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(2))
        w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3)) = w(ix1+1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,mom(3))
      enddo

      !> Floor value for negative inflow
      where(w(ixB^S,mom(1)) .lt. -0.05d0*rho_bound)
       w(ixB^S,mom(1)) = -0.05d0*rho_bound
      endwhere

      !> Ceil value for positive outflow
      where(w(ixB^S,mom(1)) .gt. 0.1d0*rho_bound)
       w(ixB^S,mom(1)) = 0.1d0*rho_bound
      endwhere

      !> per-column dE/dr from the FLD closure (D from the first interior cell)
      call inner_gradE(ixI^L, w, x, Local_gradE)

      !> Extrapolate using gradE, but impose some rational ceil value for Erad near boundary
      do ix1 = ixBmax1,ixBmin1,-1
        do ix2 = ixBmin2,ixBmax2
          do ix3 = ixBmin3,ixBmax3
            w(ix1,ix2,ix3,r_e) = min(2.d0,w(ix1+2,ix2,ix3,r_e)) &
            + (x(ix1,ix2,ix3,1)-x(ix1+2,ix2,ix3,1))*Local_gradE(ix1+1,ix2,ix3)
          enddo
        enddo
      enddo

        temp(ixB^S) = (w(ixB^S,r_e)*unit_pressure/const_rad_a)**0.25d0/unit_temperature
        w(ixB^S,e_) = w(ixB^S,rho_)*temp(ixB^S)/(eos%gamma-1.d0) + half*(w(ixB^S,mom(1))**2+w(ixB^S,mom(2))**2+w(ixB^S,mom(3))**2)/w(ixB^S,rho_)


    case(2)

      !> Ghost cells only. The Lucy value E_out for the MG solve is
      !> computed in outer_Eout, called from mg_boundary_conditions.
      !> Old Code: Stefan-Boltzmann law for outer BC, but this crashes with new Lucy E_out  
      !> w(ixB^S,r_e) = const_rad_a*(Local_Tout(ixB^S)*unit_temperature)**4.d0/unit_pressure
      !> Use 1/r^2 law of freely streaming radiation 
      do ix1 = ixBmin1,ixBmax1
        w(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,r_e) = w(ixBmin1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,r_e)*(x(ixBmin1-1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1)/x(ix1,ixBmin2:ixBmax2,ixBmin3:ixBmax3,1))**2
      enddo

      ! print*, 'E_out', E_out
    case default
      call mpistop('boundary not known')
    end select
  end subroutine boundary_conditions

  !> MG boundary values for E, the same on every MPI rank.
  !> octree-mg takes ONE value per face, so take the area-weighted mean
  !> over all leaf blocks that touch the face, on all ranks.
  subroutine mg_boundary_conditions(iB)
    use mod_global_parameters
    use mod_multigrid_coupling

    integer, intent(in) :: iB
    integer          :: iigrid, igrid
    double precision :: col(ixG^T), sendbuf(2), recvbuf(2), dA, value

    sendbuf = 0.d0
    if (mg%is_allocated) then
      do iigrid = 1, igridstail; igrid = igrids(iigrid)
        if (.not. ps(igrid)%is_physical_boundary(iB)) cycle
        select case (iB)
        case (1)
          call inner_gradE(ixG^LL, ps(igrid)%w, ps(igrid)%x, col)
        case (2)
          call outer_Eout(ixG^LL, ps(igrid)%w, ps(igrid)%x, col)
        end select
        ! every column has the same value in x1; use the interior lateral cells
        dA = rnode(rpdx2_,igrid)*rnode(rpdx3_,igrid)
        sendbuf(1) = sendbuf(1) + dA*sum(col(ixMlo1,ixMlo2:ixMhi2,ixMlo3:ixMhi3))
        sendbuf(2) = sendbuf(2) + dA*dble(block_nx2*block_nx3)
      enddo
    endif
    ! collective: fld_set_mg_bounds calls this on every rank for iB = 1, 2
    call MPI_ALLREDUCE(sendbuf, recvbuf, 2, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierrmpi)

    if (recvbuf(2) > 0.d0) then
      value = recvbuf(1)/recvbuf(2)
    else
      ! start-up (no mesh / mg yet): keep the guess from initglobaldata_usr
      value = merge(gradE, E_out, iB == 1)
    endif
    if (.not. ieee_is_finite(value)) call mpistop('WR_3D: non-finite MG boundary value')

    select case (iB)
    case (1)
      gradE = value
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_neumann
    case (2)
      E_out = value
      mg%bc(iB, mg_iphi)%bc_type = mg_bc_dirichlet
    case default
      call mpistop('WR_3D: unknown MG boundary')
    end select
    mg%bc(iB, mg_iphi)%bc_value = value
  end subroutine mg_boundary_conditions


  !> dE/dr at the inner boundary from the FLD closure (Moens+2022 Eq. 55),
  !> per lateral column, with D from the first interior cell (1D item 14)
  !> and F_adv = 4/3 v E from the same cell. Result is copied along x1.
  subroutine inner_gradE(ixI^L, w, x, gradE_col)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: gradE_col(ixI^S)

    integer          :: ixC^L, ix2, ix3
    double precision :: wprim(ixI^S,1:nw), kap(ixI^S), lam(ixI^S), fR(ixI^S)
    double precision :: D_col, F_adv

    ! first interior cell in x1; the flux limiter needs one neighbour
    ! laterally, so skip the outermost lateral ghost column on each side
    ixCmin1 = ixImin1+nghostcells; ixCmax1 = ixCmin1
    ixCmin2 = ixImin2+1;           ixCmax2 = ixImax2-1
    ixCmin3 = ixImin3+1;           ixCmax3 = ixImax3-1

    wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L, ixI^L, wprim, x)
    call fld_get_opacity_prim(wprim, x, ixI^L, ixC^L, kap, fld_fl)
    call fld_get_fluxlimiter_prim(wprim, x, ixI^L, ixC^L, lam, fR, 1, fld_fl)

    gradE_col(ixI^S) = 0.d0
    do ix3 = ixCmin3, ixCmax3
      do ix2 = ixCmin2, ixCmax2
        D_col = c_norm*lam(ixCmin1,ix2,ix3)/(kap(ixCmin1,ix2,ix3)*wprim(ixCmin1,ix2,ix3,rho_))
        F_adv = max(0.d0, 4.d0/3.d0*wprim(ixCmin1,ix2,ix3,mom(1))*wprim(ixCmin1,ix2,ix3,r_e))
        gradE_col(ixImin1:ixImax1,ix2,ix3) = max(-(F_bound-F_adv)/D_col, -200.d0)
      enddo
    enddo
    gradE_col(ixImin1:ixImax1,ixImin2,:) = gradE_col(ixImin1:ixImax1,ixImin2+1,:)
    gradE_col(ixImin1:ixImax1,ixImax2,:) = gradE_col(ixImin1:ixImax1,ixImax2-1,:)
    gradE_col(ixImin1:ixImax1,:,ixImin3) = gradE_col(ixImin1:ixImax1,:,ixImin3+1)
    gradE_col(ixImin1:ixImax1,:,ixImax3) = gradE_col(ixImin1:ixImax1,:,ixImax3-1)
  end subroutine inner_gradE

  !> Lucy (1971) E at the outer boundary (Poniatowski+2021 Eq. 7, 1D item 17),
  !> per lateral column, at the first ghost cell. Result is copied along x1.
  subroutine outer_Eout(ixI^L, w, x, Eout_col)
    use mod_global_parameters
    integer, intent(in)           :: ixI^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: Eout_col(ixI^S)

    integer          :: ixL, ixG1, ix2, ix3
    double precision :: kappa(ixI^S), kap, tau, W_geo, r

    ixL  = ixImax1-nghostcells   ! last interior cell
    ixG1 = ixL+1                 ! first ghost cell
    call get_kappa_OPAL(ixI^L, ixI^L, w, x, kappa)

    do ix3 = ixImin3, ixImax3
      do ix2 = ixImin2, ixImax2
        kap = kappa(ixL,ix2,ix3)
        if (kap .ne. kap) kap = kappa_e
        kap = min(max(kap, kappa_e), 20*kappa_e)
        r = x(ixG1,ix2,ix3,1)
        tau   = kap*w(ixG1,ix2,ix3,rho_)*R_star**2/(3*r)
        W_geo = 0.5d0*(1.d0 - sqrt(1.d0 - (R_star/r)**2))
        Eout_col(ixImin1:ixImax1,ix2,ix3) = arad_norm*F_bound/StefBoltz*(W_geo + 0.75d0*tau)
      enddo
    enddo
  end subroutine outer_Eout

  subroutine Fix_pressure(level,qt,qdt,ixI^L,ixO^L,w,x)
    use mod_global_parameters

    integer, intent(in)             :: ixI^L,ixO^L,level
    double precision, intent(in)    :: qt, qdt
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in)    :: x(ixI^S,1:ndim)

    double precision :: pth(ixI^S), mean_p

    !> fix density as well
    !where ((w(ixO^S,rho_) .lt. 1.d-5) .and. (x(ixO^S,1) .lt. 4.d0))
    !  w(ixO^S,rho_) = 1.d-5
    !endwhere

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)

    !if (any(press(ixO^S) .lt. 0.d0)) then
      mean_p = max(sum(pth(ixO^S))/(block_nx1*block_nx2*block_nx3),small_pressure)
      where (pth(ixO^S) .le. small_pressure)
        w(ixO^S,e_) = mean_p/(eos%gamma - 1) +  0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1)/w(ixO^S,rho_)
      end where
    !endif

    !> Temperature ceil, Tmax 1d6
    where ((pth(ixO^S)/w(ixO^S,rho_)*unit_temperature) .gt. 1.d6)
      w(ixO^S,e_) = 1.d6/unit_temperature*w(ixO^S,rho_)/(eos%gamma - 1) &
          + 0.5d0 * sum(w(ixO^S, mom(:))**2, dim=ndim+1)/w(ixO^S,rho_)
    end where

    where (w(ixO^S,r_e) .gt. 0.5d0)
        w(ixO^S,r_e) = 0.5d0
    end where


  end subroutine Fix_pressure


  !> Calculate gravitational acceleration in each dimension
  subroutine set_gravitation_field(ixI^L,ixO^L,wCT,x,gravity_field)
    use mod_global_parameters
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(in)    :: wCT(ixI^S,1:nw)
    double precision, intent(out)   :: gravity_field(ixI^S,1:ndim)

    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixI^S) = x(ixI^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    gravity_field(ixI^S,:) = 0.d0
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
    w(ixO^S,mom(2)) = w(ixO^S,mom(2)) + qdt*ppsource(ixO^S,mom(2)) !> OK
    w(ixO^S,mom(3)) = w(ixO^S,mom(3)) + qdt*ppsource(ixO^S,mom(3)) !> OK
    w(ixO^S,e_) = w(ixO^S,e_) + qdt*ppsource(ixO^S,e_) !> OK
    w(ixO^S,r_e) = w(ixO^S,r_e) + qdt*ppsource(ixO^S,r_e) !> TROUBLEMAKER

    if (.not. Cak_in_D) then
      if (read_cak_table) then
        call get_kappa_CAK2(ixI^L,ixO^L,wCT,x,k_cak)
      else
        call get_kappa_CAK(ixI^L,ixO^L,wCT,x,k_cak)
      endif

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

    if (read_cak_table) then
      call get_kappa_CAK2(ixI^L,ixO^L,w,x,k_cak)
    else
      call get_kappa_CAK(ixI^L,ixO^L,w,x,k_cak)
    endif

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
  !

  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)

    double precision :: rad_flux(ixI^S,1:ndim)
    double precision :: pth(ixI^S), v(ixO^S,1:ndim)
    double precision :: radius(ixO^S), pert(ixO^S)
    double precision :: edd(ixI^S,1:ndim,1:ndim), lambda(ixI^S), fld_R(ixI^S)
    double precision :: wprim(ixI^S,1:nw)
    integer :: rdir, pdir, tdir

    source(ixO^S,1:nw) = zero

    rdir = 1
    pdir = 2
    tdir = 3

    v(ixO^S,rdir) = w(ixO^S,mom(rdir))/w(ixO^S,rho_)
    v(ixO^S,pdir) = w(ixO^S,mom(pdir))/w(ixO^S,rho_)
    v(ixO^S,tdir) = w(ixO^S,mom(tdir))/w(ixO^S,rho_)

    radius(ixO^S) = x(ixO^S,rdir) ! + half*block%dx(ixO^S,rdir)

    !> Correction for spherical fluxes:
    !> drho/dt = -2 rho v_r/r
    source(ixO^S,rho_) = -two*w(ixO^S,rho_)*v(ixO^S,rdir)/radius(ixO^S)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)

    !> dm_r/dt = +(rho*v_p**2 + 2pth)/r -2 (rho*v_r**2 + pth)/r
    !> dm_phi/dt = - 3*rho*v_p m_r/r
    source(ixO^S,mom(rdir)) = - 2*w(ixO^S,rho_)*v(ixO^S,rdir)**two/radius(ixO^S) &
                              + w(ixO^S,rho_)*v(ixO^S,pdir)**two/radius(ixO^S) &
                              + w(ixO^S,rho_)*v(ixO^S,tdir)**two/radius(ixO^S)

    source(ixO^S,mom(pdir)) = - 3*v(ixO^S,rdir)*v(ixO^S,pdir)*w(ixO^S,rho_)/radius(ixO^S)
    source(ixO^S,mom(tdir)) = - 3*v(ixO^S,rdir)*v(ixO^S,tdir)*w(ixO^S,rho_)/radius(ixO^S)

    !> de/dt = -2 (e + p)v_r/r
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

    double precision :: OPAL(ixO^S), CAK(ixO^S)
    double precision :: wcons(ixI^S,1:nw)

    !> 3.3 passes PRIMITIVE w here; the helpers below expect conserved w
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

    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

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
    double precision :: n, rho0, Temp0

    !> Get OPAL opacities by reading from table
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Temp)

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        call set_opal_opacity(rho0,Temp0,n)
        kappa(ix^D) = n/unit_opacity
    {enddo\ }

    where(kappa(ixO^S) .ne. kappa(ixO^S))
      kappa(ixO^S) = kappa_e
    endwhere

    !> Lower limit is electron scattering:
    where(kappa(ixO^S) .lt. kappa_e)
      kappa(ixO^S) = kappa_e
    endwhere

    !> test without opal
    ! kappa(ixO^S) = kappa_e

  end subroutine get_kappa_OPAL

  subroutine get_kappa_CAK(ixI^L,ixO^L,w,x,kappa)
    use mod_global_parameters
    use mod_opal_opacity
    use mod_fld

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out):: kappa(ixO^S)

    integer :: ix^D
    double precision :: vel(ixI^S), gradv(ixI^S), gradvI(ixI^S)
    double precision :: xx(ixO^S), alpha(ixO^S)

    if (.not. CAK_zero) then

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

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

    if (x(ixImax1,nghostcells, nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,ixOmin2:ixOmax2,ixOmin3:ixOmax3) = kappa(ixOmax1-1,ixOmin2:ixOmax2,ixOmin3:ixOmax3)
    endif

    ! if (it .le. it_start_cak) then
    !   kappa(ixO^S) = kappa(ixO^S)*dexp(-w(ixO^S,rho_)*kappa_e)
    ! endif

    !{do ix^D=ixOmin^D,ixOmax^D\ }
    !    if (xx(ix^D) .lt. cak_x0) then
    !      kappa(ix^D) = min(2*kappa_e,kappa(ix^D))
    !    else if (xx(ix^D) .lt. cak_x1) then
    !      kappa(ix^D) = min( (2 + (20.d0 - 2.d0)/(cak_x1 - cak_x0)*(xx(ix^D)-cak_x0) )*kappa_e, kappa(ix^D))
    !    else
    !      kappa(ix^D) = min(20*kappa_e,kappa(ix^D))
    !    endif
    !
    !  ! if (kappa(ix^D) .gt. 3*kappa_e) &
    !  !   kappa(ix^D) = 3*kappa_e + half*(kappa(ix^D)-3*kappa_e)
    !{enddo\ }

    {do ix^D=ixOmin^D,ixOmax^D\ }
      kappa(ix^D) = min(50*kappa_e,kappa(ix^D))
    {enddo\ }

    else

    !> test with no cak
    kappa(ixO^S) = 0.d0

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
    double precision :: vel(ixI^S), gradv(ixI^S), gradvI(ixI^S)

    if (.not. CAK_zero) then

    !> Get CAK opacities from gradient in v_r (This is maybe a weird approximation)
    !> Need diffusion coefficient depending on direction?
    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv,nghostcells)

    ! call gradient(vel,ixI^L,ixO^L,1,gradvI)
    ! gradv(ixO^S) = gradvI(ixO^S)

    !> Absolute value of gradient:
    gradv(ixO^S) = abs(gradv(ixO^S))

    !> Get CAK opacities by reading from table
      call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,Temp)

    {do ix^D=ixOmin^D,ixOmax^D\ }
        rho0 = w(ix^D,rho_)*unit_density
        Temp0 = Temp(ix^D)*unit_temperature
        Temp0 = max(Temp0,1.d4)
        gradv0 = gradv(ix^D)*(unit_velocity/unit_length)
        call set_cak_opacity(rho0,Temp0,alpha, Qbar, Q0, kappa_e_t)

        tau = (kappa_e*unit_opacity)*rho0*const_c/gradv0
        M_t = Qbar/(1-alpha)*((1+Q0*tau)**(1-alpha) - 1)/(Q0*tau)
        kap0 = (kappa_e*unit_opacity)*M_t

        kappa(ix^D) = kap0/unit_opacity

        if (kappa(ix^D) .ne. kappa(ix^D)) kappa(ix^D) = 0.d0

        kappa(ix^D) = min(50*kappa_e,kappa(ix^D))
    {enddo\ }


    if (x(ixImax1,nghostcells,nghostcells,1) .ge. xprobmax1) then
      kappa(ixOmax1,:,:) = kappa(ixOmax1-1,:,:)
    endif

    else

    kappa(ixO^S) = 0.d0

    endif

  end subroutine get_kappa_CAK2

  subroutine ceil_diffcoef(w, wCT, x, ixI^L, ixO^L)
    use mod_global_parameters
    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(inout) :: w(ixI^S, 1:nw)
    double precision, intent(in) :: wCT(ixI^S, 1:nw)
    double precision, intent(in) :: x(ixI^S, 1:ndim)


    where (w(ixO^S,i_diff_mg) .gt. 1.5d3) &
      w(ixO^S,i_diff_mg) = 1.5d3

    where (w(ixO^S,i_diff_mg) .lt. 1.d-2) &
      w(ixO^S,i_diff_mg) = 1.d-2


  end subroutine ceil_diffcoef

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

    double precision :: lim_1, lim_2, lim_3, lim_4

    lim_3 = 1.5d0
    lim_2 = 2.5d0
    lim_1 = 4.d0

    !refine= -1
    !coarsen= -1

    if (all(x(ixG^S,1) < lim_3)) then
      if (level > 4) coarsen=1
      if (level < 4) refine=1
    elseif (all(x(ixG^S,1) < lim_2)) then
      if (level > 3) coarsen=1
      if (level < 3) refine=1
    elseif (all(x(ixG^S,1) < lim_1)) then
      if (level > 2) coarsen=1
      if (level < 2) refine=1
    endif

  end subroutine refine_base


  subroutine collapse_to_1D()
    use mod_global_parameters

    integer          :: iigrid, igrid, jj_blk,ii, jj, nbx, i, il, ih, dn_mdot
    integer          :: lvl, ibx
    integer          :: np_mdot, nc
    double precision :: ratio, sf_mdot, dx_l1

    double precision :: radflux(ixG^T,ndim)

    double precision, allocatable :: rp_mdot(:), p_mdot(:)
    double precision, allocatable :: rp_lum(:),p_lum(:)
    double precision, allocatable :: mdot_S(:), mdot_R(:)
    double precision, allocatable :: lum_S(:), lum_R(:)
    integer, allocatable :: jp_mdot(:)

    double precision :: rr(1:domain_nx1), rr_S(1:domain_nx1), rr_R(1:domain_nx1)
    double precision :: vr(1:domain_nx1), vr_S(1:domain_nx1), vr_R(1:domain_nx1)
    double precision :: mdot(1:domain_nx1), lum(1:domain_nx1)
    double precision :: rho(1:domain_nx1), rho_S(1:domain_nx1), rho_R(1:domain_nx1)
    double precision :: rho2(1:domain_nx1), rho2_S(1:domain_nx1), rho2_R(1:domain_nx1)
    double precision :: vr2(1:domain_nx1), vr2_S(1:domain_nx1), vr2_R(1:domain_nx1)
    double precision :: rhovr(1:domain_nx1), rhovr_S(1:domain_nx1), rhovr_R(1:domain_nx1)
    double precision :: rho2vr(1:domain_nx1), rho2vr_S(1:domain_nx1), rho2vr_R(1:domain_nx1)

    integer :: lvl_h(1:domain_nx1), lvl_h_S(1:domain_nx1), lvl_h_R(1:domain_nx1)
    integer :: lvl_l(1:domain_nx1), lvl_l_S(1:domain_nx1), lvl_l_R(1:domain_nx1)

    ! if (refine_max_level .ne. 1)
    ! call mpistop("collapse_to_1D doesnt work YET with mpi")

    !> #R_star -1 in simulation
    np_mdot = floor((xprobmax1-xprobmin1)/R_star)

    allocate(rp_mdot(1:np_mdot))
    allocate(p_mdot(1:np_mdot))
    allocate(p_lum(1:np_mdot))
    allocate(jp_mdot(1:np_mdot))
    allocate(mdot_S(1:np_mdot))
    allocate(mdot_R(1:np_mdot))
    allocate(lum_S(1:np_mdot))
    allocate(lum_R(1:np_mdot))

    rr = 0.d0
    rp_mdot = 0.d0
    vr = 0.0d0
    rho = 0.d0
    rho2 = 0.d0
    vr2 = 0.d0
    mdot = 0.d0
    p_mdot = 0.d0
    lum = 0.d0
    p_lum = 0.d0
    lvl_h = 0.d0
    lvl_l = 2* refine_max_level
    rhovr = 0.d0
    rho2vr = 0.d0

    !> Reconstruct radius at level 1
    do jj = 1,domain_nx1
      rr(jj) = xprobmin1 + (jj-0.5d0)*(xprobmax1-xprobmin1)/domain_nx1
    enddo

    !> Choose radii at which to save mdot
    !> This is done at every stellar radii
    do jj = 1,np_mdot
      ! rp_mdot(jj) = xprobmin1 + (jj-1)* R_star
      dx_l1 = (xprobmax1-xprobmin1)/domain_nx1
      dn_mdot = floor(1.d0*domain_nx1/np_mdot+smalldouble)
      rp_mdot(jj) = xprobmin1 + dx_l1/2 + (jj-1)*dn_mdot*dx_l1
    enddo

    !> Find cells that are closest to radii where we want to track mdot
    do jj = 1,domain_nx1
      do ii = 1,np_mdot
        if (rr(jj) - rp_mdot(ii) .le. smalldouble) &
          jp_mdot(ii) = jj
          ! rp_mdot(ii) = rr(jj)
      enddo
    enddo

    !> loop over all grids for this proc
    do iigrid=1,igridstail; igrid=igrids(iigrid);
      block=>ps(igrid)

      !> widht in cells per cell of lvl 1: On lvl 1: 1, lvl 2: 2, lvl 3: 4
      nc = 2**(node(plevel_,igrid)-1)

      !> Map the block index to the global index
      !> nbx = nr of blocks between i=0 and current block
      jj_blk = (node(pig1_,igrid)-1)*block_nx1/nc

      if (nc > block_nx1) &
        call mpistop("collapse_to_1D doesnt work, reduce amr")

      !> Calculate radflux in block for luminosity
      call fld_get_radflux(block%w, block%x, ixG^LL, ixM^LL, radflux, fld_fl)

      !> For all cells in the current block, average velocity over lateral direction.
      !> Take into account possible amr
      do ii = 1,block_nx1/nc
        !> relate global index jj to local index ii on grid
        il = nghostcells + 1 + (ii-1)*nc
        ih = nghostcells + 1 + (ii)*nc -1
        jj = jj_blk + ii

        !> Velocity
        vr(jj) = vr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,mom(1))/block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_))/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1

        !> Density
        rho(jj) = rho(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_))/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1

        !> Density squared
        rho2(jj) = rho2(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_)**2)/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1

        !> radial velocity squared
        vr2(jj) = vr2(jj) &
        + sum((block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,mom(1))/block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_))**2)/(domain_nx2*domain_nx3)/nc**3

        !> Density weighted velocity
        rhovr(jj) = rhovr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,mom(1)))/(domain_nx2*domain_nx3)/nc**3

        !> Density squared weighted velocity
        rho2vr(jj) = rho2vr(jj) &
        + sum(block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,mom(1))*block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,rho_))/(domain_nx2*domain_nx3)/nc**3

        !> Mass loss rate
        mdot(jj) = mdot(jj) &
        + 4*dpi*sum(block%x(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)**2*block%w(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,mom(1)))/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1

        !> Luminosity
        lum(jj) = lum(jj) &
        + 4*dpi*sum(block%x(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1)**2*radflux(il:ih,ixMlo2:ixMhi2,ixMlo3:ixMhi3,1))/(domain_nx2*domain_nx3)/nc**3 !> average value on lvl 1

        !> Highest/Lowest amr level
        lvl_h(jj) = node(plevel_,igrid)
        lvl_l(jj) = node(plevel_,igrid)
      enddo

    enddo

    !> communicate velocity
    vr_S=vr
    call mpi_allreduce(vr_S,vr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    vr=vr_R

    !> communicate density
    rho_S=rho
    call mpi_allreduce(rho_S,rho_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho=rho_R

    !> communicate density squared
    rho2_S=rho2
    call mpi_allreduce(rho2_S,rho2_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho2=rho2_R

    !> communicate vr squared
    vr2_S=vr2
    call mpi_allreduce(vr2_S,vr2_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    vr2=vr2_R

    !> communicate rho*vr
    rhovr_S=rhovr
    call mpi_allreduce(rhovr_S,rhovr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rhovr=rhovr_R

    !> communicate rho*vr
    rho2vr_S=rho2vr
    call mpi_allreduce(rho2vr_S,rho2vr_R,domain_nx1,&
    MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierrmpi)
    rho2vr=rho2vr_R

    !> Only keep mdot in interested radii
    do ii = 1,np_mdot
      jj = jp_mdot(ii)
      p_mdot(ii) = mdot(jj)
      p_lum(ii) = lum(jj)
    enddo

    !> communicate mdot array
    mdot_S=p_mdot
    call mpi_reduce(mdot_S,mdot_R,np_mdot,&
    MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    p_mdot=mdot_R*unit_density*unit_velocity*unit_length**2/M_sun*year

    !> communicate lum array
    lum_S=p_lum
    call mpi_reduce(lum_S,lum_R,np_mdot,&
    MPI_DOUBLE_PRECISION,MPI_SUM,0,icomm,ierrmpi)
    p_lum=lum_R*unit_radflux*unit_length**2/L_sun

    !> communicate highest amr level
    lvl_h_S=lvl_h
    call mpi_allreduce(lvl_h_S,lvl_h_R,domain_nx1,&
    MPI_INTEGER,MPI_MAX,icomm,ierrmpi)
    lvl_h = lvl_h_R
    lvl_l_S=lvl_l
    call mpi_allreduce(lvl_l_S,lvl_l_R,domain_nx1,&
    MPI_INTEGER,MPI_MIN,icomm,ierrmpi)
    lvl_l = lvl_l_R

    !> integrate over time
    sum_time = sum_time + dt
    vr_sumt(1:domain_nx1) = vr_sumt(1:domain_nx1) + vr(1:domain_nx1)*dt
    rho_sumt(1:domain_nx1) = rho_sumt(1:domain_nx1) + rho(1:domain_nx1)*dt
    rho2_sumt(1:domain_nx1) = rho2_sumt(1:domain_nx1) + rho2(1:domain_nx1)*dt
    vr2_sumt(1:domain_nx1) = vr2_sumt(1:domain_nx1) + vr2(1:domain_nx1)*dt
    rhovr_sumt(1:domain_nx1) = rhovr_sumt(1:domain_nx1) + rhovr(1:domain_nx1)*dt
    rho2vr_sumt(1:domain_nx1) = rho2vr_sumt(1:domain_nx1) + rho2vr(1:domain_nx1)*dt

    !> Write out average velocity profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_vr',status='replace')
      write(unitanalysis,*) 'r | vr_int_theta | vr_int_theta_int_t | lvl_h lvl_l'
      do i=1,domain_nx1
        write(unitanalysis,'(8f11.7,2i4)') rr(i), vr(i), vr_sumt(i)/sum_time, vr2(i), vr2_sumt(i)/sum_time,dsqrt( abs(vr2_sumt(i)/sum_time - (vr_sumt(i)/sum_time)**2 )), rhovr_sumt(i)/rho_sumt(i), rho2vr_sumt(i)/rho2_sumt(i) ,lvl_h(i), lvl_l(i)
enddo
    close(unitanalysis)
    endif

    !> Write out average density profile
    if (mype==0) then
      !> Always update file to give the last mean snapshot
      open(unit=unitanalysis,file=trim(base_filename)//'_rho',status='replace')
      write(unitanalysis,*) 'r | rho_int_theta | rho_int_theta_int_t'
      do i=1,domain_nx1
        write(unitanalysis,'(6f11.7)') rr(i), rho(i), rho_sumt(i)/sum_time, rho2(i), rho2_sumt(i)/sum_time, sum_time*rho2_sumt(i)/rho_sumt(i)**2
      enddo
    close(unitanalysis)
    endif

    !> Write out average mdot at different radii
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'Mdot at t | r= ', rp_mdot
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_mdot',access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_mdot
       endif
       close(unitanalysis)
     endif

    !> Write out average lum at different radii
    if (mype==0) then
       if (global_time<smalldouble) then ! if very 1st iteration
         open(unit=unitanalysis,file=trim(base_filename)//'_lum',status='replace')
         write(unitanalysis,'(a16,*(f7.2))') 'Lum at t | r= ', rp_mdot
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_lum
       else
         open(unit=unitanalysis,file=trim(base_filename)//'_lum',access='append')
         write(unitanalysis,'(f10.3,*(E15.6))') global_time, p_lum
       endif
       close(unitanalysis)
     endif

  end subroutine collapse_to_1D


  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
    double precision                   :: g_grav(ixI^S), pth(ixI^S)
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: wprim(ixI^S,1:nw), OPAL(ixO^S), CAK(ixO^S), CAK2(ixO^S)
    double precision                   :: vel(ixI^S), gradv(ixI^S)
    double precision                   :: rad_flux(ixI^S,1:ndim), Lum(ixO^S)
    double precision                   :: pp_rf(ixO^S), lambda(ixI^S), fld_R(ixI^S), gradOE(ixI^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)

    wprim(ixI^S,1:nw) = w(ixI^S,1:nw)
    call phys_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_radflux(w, x, ixI^L, ixO^L, rad_flux, fld_fl)

    call eos%get_thermal_pressure(w,x,ixI^L,ixO^L,pth)
    call eos%get_temperature_from_etot(w, x, ixI^L, ixO^L, Tgas)
    call hd_get_trad(w, x, ixI^L, ixO^L, Trad)

    call get_kappa_OPAL(ixI^L,ixO^L,w,x,OPAL)
    call get_kappa_CAK(ixI^L,ixO^L,w,x,CAK)
    call get_kappa_CAK2(ixI^L,ixO^L,w,x,CAK2)

    if (read_cak_table) then
      g_rad(ixO^S) = (OPAL(ixO^S)+CAK2(ixO^S))*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    else  
      g_rad(ixO^S) = (OPAL(ixO^S)+CAK(ixO^S))*rad_flux(ixO^S,1)/(const_c/unit_velocity)
    endif

    g_grav(ixO^S) = const_G*mass/radius(ixO^S)**2*(unit_time**2/unit_length)
    big_gamma(ixO^S) = g_rad(ixO^S)/g_grav(ixO^S)

    vel(ixI^S) = w(ixI^S,mom(1))/w(ixI^S,rho_)
    call gradient(vel,ixI^L,ixO^L,1,gradv)

    pp_rf(ixO^S) = two*rad_flux(ixO^S,1)/x(ixO^S,1)*dt

    call fld_get_fluxlimiter(w, x, ixI^L, ixO^L, lambda, fld_R, 1, fld_fl)

    Lum(ixO^S) = 4*dpi*rad_flux(ixO^S,1)*(x(ixO^S,1)*unit_length)**2*unit_radflux/L_sun

    call gradient(w(ixI^S,r_e),ixI^L,ixO^L,1,gradOE,nghostcells)

    w(ixO^S,i_v1) = w(ixO^S,mom(1))/w(ixO^S,rho_)
    w(ixO^S,i_v2) = w(ixO^S,mom(2))/w(ixO^S,rho_)
    w(ixO^S,i_v3) = w(ixO^S,mom(3))/w(ixO^S,rho_)
    w(ixO^S,i_p) = pth(ixO^S)

    w(ixO^S,i_Trad) = Trad(ixO^S)*unit_temperature
    w(ixO^S,i_Tgas) = Tgas(ixO^S)*unit_temperature
    w(ixO^S,i_Mdot) = 4*dpi*w(ixO^S,mom(1))*radius(ixO^S)**2 &
    *unit_density*unit_velocity/M_sun*year
    w(ixO^S,i_Opal) = OPAL(ixO^S)/kappa_e
    w(ixO^S,i_CAK) = CAK(ixO^S)/kappa_e
    w(ixO^S,i_CAK2) = CAK2(ixO^S)/kappa_e
    w(ixO^S,i_lambda) = lambda(ixO^S)
    w(ixO^S,i_Edd) = lambda(ixO^S) + lambda(ixO^S)**2 * fld_R(ixO^S)**2
    w(ixO^S,i_Gamma) = big_gamma(ixO^S)
    w(ixO^S,i_Lum) = Lum(ixO^S)
    w(ixO^S,i_F1) = rad_flux(ixO^S,1)/F_bound
    w(ixO^S,i_F2) = rad_flux(ixO^S,2)/F_bound
    w(ixO^S,i_F3) = rad_flux(ixO^S,3)/F_bound
    w(ixO^S,i_gradE) = gradOE(ixO^S)

  end subroutine update_extravars

end module mod_usr
