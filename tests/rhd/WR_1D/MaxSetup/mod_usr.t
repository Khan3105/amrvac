! WR_1D for amrvac3.3, commit 76aed67aa5607ef39410615d64ffb5686d9b274e.
! Apply rhd-src.patch and retain the repository 1D_stable.blk.
module mod_usr

  ! Include a physics module
  use mod_hd
  use mod_eos, only: eos
  use mod_fld

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  double precision, parameter :: M_sun = 1.9891000d33
  double precision, parameter :: R_sun = 6.9599000d10
  double precision, parameter :: L_sun = 3.8268000d33
  double precision, parameter :: year = 365.25*24*60*60

  double precision :: StefBoltz

  double precision :: cak_Q=2300.d0, cak_a=0.5d0, cak_base=0.66d0
  double precision :: cak_x0=0.4d0, cak_x1=0.65d0
  integer :: it_start_cak=0 ! accepted for compatibility; the original hook was inactive
  double precision :: rho_bound=2.5d-8, v_inf=1.d8, Mdot=1.d-5
  double precision :: T_bound=3.2d5, R_star=1.d0, M_star=10.d0

  integer :: i_v1, i_v2, i_p
  integer :: i_Trad, i_Tgas, i_Mdot, i_Opal, i_CAK, i_lambda, i_fld_R
  integer :: i_Gamma, i_Lum, i_F1, i_F2, i_grE

  double precision :: kappa_e=0.2d0, Gamma_e_bound=0.4d0
  double precision :: L_bound, F_bound, gradE, E_out
  logical :: fixed_lum=.false., Cak_in_D=.true., read_cak_table=.false.
  double precision, allocatable :: initial_profile(:,:)
  integer :: n_profile

contains

  !> This routine should set user methods, and activate the physics module
  subroutine usr_init()

    ! Use a Cartesian operator plus the spherical divergence corrections below.
    call set_coordinate_system("Cartesian_1D")

    ! Set independent units BEFORE hd_activate derives the EOS/FLD constants.
    ! eos_init has already read eos_list; params_read reads wind_list in cgs.
    call params_read(par_files)
    unit_length = R_star
    unit_density = rho_bound
    unit_velocity = v_inf

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

    ! Receives the freshly computed physical D before the collective MG BC hook.
    usr_special_diffcoef => radiation_diffcoef

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

  !> Set units-dependent wind parameters after EOS finalisation.
  subroutine initglobaldata_usr
    use mod_opal_opacity, only: init_opal_table
    use mod_cak_opacity, only: init_cak_table

    if (eos%eos_type /= 'FI') &
         call mpistop('WR_1D requires eos_type=FI: FLD exchange assumes constant gamma')
    if (SI_unit) call mpistop('WR_1D opacity tables and input parameters use cgs')
    if (.not.hd_energy .or. .not.hd_radiation_fld .or. .not.hd_gravity) &
         call mpistop('WR_1D requires HD energy, FLD and gravity')
    if (ndim /= 1 .or. .not.slab) call mpistop('WR_1D requires Cartesian_1D')
    if (cak_Q <= 0.d0 .or. min(cak_a,cak_base) <= 0.d0 .or. &
         max(cak_a,cak_base) >= 1.d0 .or. cak_x1 <= cak_x0) &
         call mpistop('Invalid CAK parameters')

    call init_opal_table(fld_opal_table)
    if (read_cak_table) call init_cak_table(fld_opal_table)

    R_star = R_star/unit_length
    M_star = M_star/unit_mass
    T_bound = T_bound/unit_temperature
    rho_bound = rho_bound/unit_density
    v_inf = v_inf/unit_velocity
    Mdot = Mdot*unit_time/unit_mass
    kappa_e = kappa_e/unit_opacity
    F_bound = F_bound/unit_radflux
    L_bound = L_bound/(unit_radflux*unit_length**2)
    StefBoltz = const_sigmaSB*unit_temperature**4/unit_radflux

    ! Defined on EVERY rank before the startup multigrid callback.
    gradE = -3.d0*kappa_e*rho_bound*F_bound/c_norm
    E_out = arad_norm*(1.5d4/unit_temperature)**4
    call load_initial_profile()

    if (mype == 0) then
       print *, 'WR_1D: unit length, density, velocity:', &
            unit_length,unit_density,unit_velocity
       print *, 'WR_1D: unit temperature, c_norm, arad_norm:', &
            unit_temperature,c_norm,arad_norm
       print *, 'WR_1D: L/Lsun, F_bound:', &
            L_bound*unit_radflux*unit_length**2/L_sun,F_bound
    endif
  end subroutine initglobaldata_usr

  subroutine params_read(files)
    use mod_global_parameters, only: unitpar
    character(len=*), intent(in) :: files(:)
    integer                      :: n

    namelist /wind_list/ cak_Q, cak_a, cak_base, cak_x0, cak_x1, rho_bound, kappa_e, &
    T_bound, R_star, M_star, v_inf, Mdot, Gamma_e_bound, it_start_cak, fixed_lum, Cak_in_D, read_cak_table

    do n = 1, size(files)
       open(unitpar, file=trim(files(n)), status="old")
       rewind(unitpar)
       read(unitpar, wind_list, end=113)
113    close(unitpar)
    end do

    if (min(R_star,M_star,rho_bound,v_inf,kappa_e,T_bound) <= 0.d0 .or. &
         Gamma_e_bound <= 0.d0) call mpistop('Invalid WR wind normalization')

    R_star = R_star*R_sun
    M_star = M_star*M_sun
    L_bound = Gamma_e_bound * 4.0 * dpi * const_G * M_star * const_c/kappa_e
    F_bound = L_bound/(4*dpi*R_star**2)
    Mdot = Mdot*M_sun/year

  end subroutine params_read

  !> Read the existing test profile once per rank, using its actual row count.
  subroutine load_initial_profile()
    integer :: iu, ios, j, n_cells
    double precision :: profile_time
    character(len=256) :: header

    open(newunit=iu,file='1D_stable.blk',status='old',action='read',iostat=ios)
    if (ios /= 0) call mpistop('WR_1D: keep the repository 1D_stable.blk beside amrvac.par')
    read(iu,'(A)',iostat=ios) header
    if (ios /= 0) call mpistop('WR_1D: invalid profile header')
    read(iu,*,iostat=ios) n_profile,n_cells
    if (ios /= 0 .or. n_profile < 2) call mpistop('WR_1D: invalid profile row count')
    read(iu,*,iostat=ios) profile_time
    if (ios /= 0) call mpistop('WR_1D: invalid profile time')
    allocate(initial_profile(6,n_profile))
    do j=1,n_profile
       read(iu,*,iostat=ios) initial_profile(:,j)
       if (ios /= 0) call mpistop('WR_1D: incomplete initial profile')
    enddo
    close(iu)
    if (any(.not.ieee_is_finite(initial_profile))) call mpistop('WR_1D: non-finite profile')
    if (any(initial_profile(2,:) <= 0.d0) .or. &
         any(initial_profile(4:6,:) <= 0.d0)) call mpistop('WR_1D: non-positive profile')
    if (any(initial_profile(1,2:) <= initial_profile(1,:n_profile-1))) &
         call mpistop('WR_1D: profile coordinates must increase')
  end subroutine load_initial_profile

  subroutine initial_conditions(ixI^L,ixO^L,w,x)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: row(6), fraction
    integer :: j, lo, hi, mid

    w(ixO^S,1:nw) = 0.d0
    do j=ixOmin1,ixOmax1
       lo=1; hi=n_profile
       do while (hi-lo > 1)
          mid=(lo+hi)/2
          if (initial_profile(1,mid) <= x(j,1)) then
             lo=mid
          else
             hi=mid
          endif
       enddo
       fraction=(x(j,1)-initial_profile(1,lo))/ &
            (initial_profile(1,hi)-initial_profile(1,lo))
       ! Clamp at the two endpoints; do not extrapolate the profile from r=9.
       fraction=max(0.d0,min(1.d0,fraction))
       row=initial_profile(:,lo)+fraction*(initial_profile(:,hi)-initial_profile(:,lo))
       w(j,rho_)=row(2)
       w(j,mom(1))=row(3) ! primitive velocity in the input file
       w(j,e_)=row(4)     ! primitive gas pressure in the input file
       w(j,r_e)=row(5)
       w(j,i_diff_mg)=row(6) ! startup guess only; FLD recomputes it each step
    enddo
    call hd_to_conserved(ixI^L,ixO^L,w,x)
  end subroutine initial_conditions

  !> Corrected weighted-optical-depth outer condition, evaluated at the face.
  subroutine outer_radiation_value(ixI^L,j,w,x,value)
    use mod_opal_opacity, only: set_opal_opacity
    integer, intent(in) :: ixI^L,j
    double precision, intent(in) :: w(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision, intent(out) :: value
    double precision :: temp(ixI^S),kappa_cgs,kappa_out,rho_face,tau_sc
    integer :: ixO^L

    ixOmin1=j; ixOmax1=j
    call eos%get_temperature_from_etot(w,x,ixI^L,ixO^L,temp)
    call set_opal_opacity(w(j,rho_)*unit_density, &
         max(1.d4,temp(j)*unit_temperature),kappa_cgs)
    kappa_out=min(20.d0*kappa_e,max(kappa_e,kappa_cgs/unit_opacity))
    rho_face=w(j,rho_)*(x(j,1)/xprobmax1)**2
    tau_sc=kappa_out*rho_face*R_star**2/(3.d0*xprobmax1)
    ! T^4 = (F_bound/sigma_SB) * (3*tau_sc/4); E=a*T^4.
    value=max(arad_norm*(1.5d4/unit_temperature)**4, &
         arad_norm*(F_bound/StefBoltz)*(3.d0*tau_sc/4.d0))
    if (.not.ieee_is_finite(value) .or. value <= 0.d0) &
         call mpistop('WR_1D: invalid outer radiation boundary')
  end subroutine outer_radiation_value

  subroutine boundary_conditions(qdt,qt,ixI^L,ixB^L,iB,w,x)
    integer, intent(in) :: ixI^L,ixB^L,iB
    double precision, intent(in) :: qdt,qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision :: temp(ixI^S),Rfactor(ixI^S),D_adj,F_adv
    integer :: j,first,mirror

    select case(iB)
    case(1)
       first=ixBmax1+1
       w(ixBmax1,rho_)=rho_bound
       do j=ixBmax1-1,ixBmin1,-1
          w(j,rho_)=exp(2.d0*log(w(j+1,rho_))-log(w(j+2,rho_)))
       enddo
       w(ixB^S,mom(1))=max(0.d0,w(first,mom(1)))
       ! F_adv and F_bound are both FLUXES, with no 4*pi*r^2 here.
       F_adv=(4.d0/3.d0)*w(ixBmax1,mom(1))/w(ixBmax1,rho_)*w(first,r_e)
       D_adj=w(first,i_diff_mg)
       if (.not.ieee_is_finite(D_adj) .or. D_adj <= 0.d0) &
            call mpistop('WR_1D: invalid physical D at inner boundary; apply the FLD patch')
       gradE=-(F_bound-F_adv)/D_adj
       do j=ixBmax1,ixBmin1,-1
          w(j,r_e)=w(first,r_e)+(x(j,1)-x(first,1))*gradE
       enddo
       if (any(.not.ieee_is_finite(w(ixB^S,r_e))) .or. any(w(ixB^S,r_e)<=0.d0)) &
            call mpistop('WR_1D: non-positive inner radiation energy')
       temp(ixB^S)=(w(ixB^S,r_e)/arad_norm)**0.25d0
       call eos%get_Rfactor(w,x,ixI^L,ixB^L,Rfactor)
       w(ixB^S,e_)=w(ixB^S,rho_)*Rfactor(ixB^S)*temp(ixB^S)*eos%inv_gamma_minus_1 &
            +0.5d0*w(ixB^S,mom(1))**2/w(ixB^S,rho_)
       w(ixB^S,i_diff_mg)=D_adj
    case(2)
       first=ixBmin1-1
       call outer_radiation_value(ixI^L,first,w,x,E_out)
       do j=ixBmin1,ixBmax1
          mirror=2*ixBmin1-j-1
          w(j,r_e)=max(small_r_e,2.d0*E_out-w(mirror,r_e))
       enddo
       w(ixB^S,i_diff_mg)=w(first,i_diff_mg)
    case default
       call mpistop('WR_1D: unknown boundary')
    end select
  end subroutine boundary_conditions

  !> Called per block after computing physical D; no MPI collectives here.
  subroutine radiation_diffcoef(w,wprim,x,ixI^L,ixO^L)
    integer, intent(in) :: ixI^L,ixO^L
    double precision, intent(inout) :: w(ixI^S,1:nw)
    double precision, intent(in) :: wprim(ixI^S,1:nw),x(ixI^S,1:ndim)
    double precision :: F_adv
    integer :: j

    if (any(.not.ieee_is_finite(w(ixO^S,i_diff_mg))) .or. &
         any(w(ixO^S,i_diff_mg)<=0.d0)) call mpistop('WR_1D: invalid diffusion coefficient')
    if (block%is_physical_boundary(1)) then
       j=ixMlo1
       F_adv=(4.d0/3.d0)*w(j-1,mom(1))/w(j-1,rho_)*w(j,r_e)
       gradE=-(F_bound-F_adv)/w(j,i_diff_mg)
    endif
    if (block%is_physical_boundary(2)) &
         call outer_radiation_value(ixI^L,ixMhi1,w,x,E_out)
  end subroutine radiation_diffcoef

  !> Invoked collectively at startup and, with the FLD patch, before EACH solve.
  subroutine mg_boundary_conditions(iB)
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


  subroutine PseudoPlanarSource(ixI^L,ixO^L,w,x,source)
    use mod_global_parameters

    integer, intent(in)           :: ixI^L, ixO^L
    double precision, intent(in)  :: w(ixI^S,1:nw), x(ixI^S,1:ndim)
    double precision, intent(out) :: source(ixO^S,1:nw)

    double precision :: rad_flux(ixI^S,1:ndir)
    double precision :: pth(ixI^S),v(ixO^S,1:ndim)
    double precision :: radius(ixO^S),  pert(ixO^S)

    double precision :: edd(ixI^S,1:ndim,1:ndim)
    double precision :: lambda(ixI^S), fld_R(ixI^S)
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

    wprim = w
    call hd_to_primitive(ixI^L,ixI^L,wprim,x)
    call fld_get_eddington(wprim,x,ixI^L,ixO^L,edd,lambda,fld_R,2,fld_fl)
    ! Missing transverse radiation work in spherical geometry:
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

    wcons = w
    call hd_to_conserved(ixI^L,ixI^L,wcons,x)

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




  subroutine update_extravars(igrid,level,ixI^L,ixO^L,qt,w,x)
    use mod_global_parameters
    integer, intent(in)             :: igrid,level,ixI^L,ixO^L
    double precision, intent(in)    :: qt,x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                   :: g_rad(ixI^S), big_gamma(ixI^S)
    double precision                   :: g_grav(ixI^S)
    double precision                   :: Tgas(ixI^S),Trad(ixI^S)
    double precision                   :: kappa(ixO^S), OPAL(ixO^S), CAK(ixO^S)
    double precision                   :: vel(ixI^S), gradv(ixI^S), gradE(ixI^S)
    double precision                   :: rad_flux(ixI^S,1:ndim), Lum(ixO^S)
    double precision                   :: pth(ixI^S), lambda(ixI^S), fld_R(ixI^S)
    integer                            :: idim
    double precision :: radius(ixI^S)
    double precision :: mass

    radius(ixO^S) = x(ixO^S,1)*unit_length
    mass = M_star*(unit_density*unit_length**3.d0)
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
