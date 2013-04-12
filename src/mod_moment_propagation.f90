MODULE MOMENT_PROPAGATION
    use precision_kinds
    use constants, only: x, y, z
    use system, only: tmax, tmom, pbc, supercell
    use mod_lbmodel, only: lbm
    implicit none
    private
    public :: init, propagate, deallocate_propagated_quantity!, not_yet_converged
    real(dp), allocatable, dimension(:,:,:,:,:) :: Propagated_Quantity
    real(dp), allocatable, dimension(:,:,:,:,:) :: Propagated_Quantity_Adsorbed
    integer(i2b), parameter :: now=0, next=1, past=-1
    real(dp), dimension(x:z, past:next) :: vacf
    real(dp) :: lambda, lambda_s ! lambda bulk and surface
    type type_tracer
        real(dp) :: ka, kd, K, z, Db, Ds !K=ka/kd, z=tracer charge
    end type
    type(type_tracer) :: tracer
    contains

! ==============================================================================

SUBROUTINE INIT
  use system, only: phi, fluid, solid,&
                     lx, ly, lz, rho, n
  use look_at_supercell, only: is_interfacial
  use input, only: input_dp
  real(dp) :: boltz_weight, Pstat, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi
  integer(i2b) :: i, j, k, l, l_inv, ip, jp, kp

  tracer%ka = input_dp('tracer_ka')
  tracer%kd = input_dp('tracer_kd')
  if( .not. test(tracer%ka) ) stop 'problem in tracer%ka in module moment_propagation'
  if( .not. test(tracer%kd) ) stop 'problem in tracer%kd in module moment_propagation'
  if(tracer%kd==0.0_dp) then
    tracer%K = 0.0_dp
  else
    tracer%K = (tracer%ka)/(tracer%kd)
  end if
  tracer%z = input_dp('tracer_z')
  tracer%Db = input_dp('tracer_Db') ! bulk diffusion coefficient of tracer, i.e. the molecular diffusion coefficient
  tracer%Ds = input_dp('tracer_Ds') ! surface diffusion coefficient of tracer
  if (tracer%Db <= 0.0_dp ) stop 'tracer_Db as readen in input is invalid'
  if (tracer%Ds < 0.0_dp ) stop 'tracer_Ds as readen in input file is invalid'

  lambda = calc_lambda()
  lambda_s = calc_lambda_s()

  vacf = 0.0_dp

  call test_and_allocate_what_is_needed_for_moment_propagation

  ! the sum of all boltzman weights is the sum over all exp(-z*phi) where supercell%node%nature == fluid. Note that is_interfacial == fluid + at interface
  Pstat = sum(exp(-tracer%z*phi), mask=(supercell%node%nature==fluid))&
         +sum(exp(-tracer%z*phi)*tracer%K, mask=(is_interfacial))

  do concurrent (i=1:lx, j=1:ly, k=1:lz, supercell%node(i,j,k)%nature==fluid )
    boltz_weight = exp(-tracer%z*phi(i,j,k))/Pstat ! boltz_weight=1/Pstat if tracer%z=0
    do concurrent (l=lbm%lmin+1:lbm%lmax)
      ip = pbc (i+lbm%vel(l)%coo(x) ,x)
      jp = pbc (j+lbm%vel(l)%coo(y) ,y)
      kp = pbc (k+lbm%vel(l)%coo(z) ,z)
      if (supercell%node(ip,jp,kp)%nature==solid) cycle
      exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
      exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
      fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
      scattprop = calc_scattprop( n(i,j,k,l), rho(i,j,k), lbm%vel(l)%a0, lambda, fermi)
      vacf(:,past) = vacf(:,past) + boltz_weight * scattprop * lbm%vel(l)%coo(:)**2
      l_inv = lbm%vel(l)%inv
      scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv), rho(ip,jp,kp), lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)
      Propagated_Quantity(:,i,j,k,now) = Propagated_Quantity(:,i,j,k,now) &
                 + exp_min_dphi * scattprop_p * lbm%vel(l_inv)%coo(:) * boltz_weight
    end do
    if(is_interfacial(i,j,k)) Propagated_Quantity_Adsorbed(:,i,j,k,now) = 0.0_dp
  end do

  print*, 0, vacf(x,past), vacf(y,past), vacf(z,past)
  open(unit=99, file='output/vacf.dat')
  write(99,*)'# time t, VACF_x(t), VACF_y(t), VACF_z(t)'
  write(99,*) 0, vacf(x,past), vacf(y,past), vacf(z,past)
  close(99)
END SUBROUTINE INIT

! ==============================================================================

SUBROUTINE PROPAGATE(it, is_converged)
  use system, only: lx, ly, lz, supercell, fluid, solid, n, rho
  use look_at_supercell, only: is_interfacial
  integer(kind=i2b), intent(in) :: it
  real(kind=dp) :: fermi, restpart, scattprop, scattprop_p, exp_dphi
  integer(kind=i2b), parameter :: now=0, next=1, past=-1
  real(dp), dimension(3) :: u_star
  integer(kind=i2b) :: i, j, k, l, l_inv, ip, jp, kp
  logical :: error
  logical, intent(out) :: is_converged

  lambda = calc_lambda()
  lambda_s = calc_lambda_s()

  error=.false.

  fluidnodes: do concurrent (i=1:lx, j=1:ly, k=1:lz, supercell%node(i,j,k)%nature==fluid )
    u_star = 0.0_dp ! the average velocity at r
    restpart = 1.0_dp ! fraction of particles staying at r; decreases in the loop over neighbours

    do l = lbm%lmin+1, lbm%lmax
      ip = pbc(i+lbm%vel(l)%coo(x) ,x)
      jp = pbc(j+lbm%vel(l)%coo(y) ,y)
      kp = pbc(k+lbm%vel(l)%coo(z) ,z)
      if ( supercell%node(ip,jp,kp)%nature /= fluid ) cycle
      fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
      scattprop = calc_scattprop( n(i,j,k,l), rho(i,j,k), lbm%vel(l)%a0, lambda, fermi)
      restpart = restpart - scattprop
      scattprop_p = calc_scattprop( n(ip,jp,kp,lbm%vel(l)%inv), rho(ip,jp,kp), lbm%vel((lbm%vel(l)%inv))%a0, lambda, 1.0_dp-fermi)
      Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity(:,i,j,k,next) + Propagated_Quantity(:,ip,jp,kp,now)*scattprop_p
      u_star = u_star + scattprop * lbm%vel(l)%coo(:)
    end do

    vacf(:,now) = vacf(:,now) + Propagated_Quantity(:,i,j,k,now)*u_star(:)

    if (is_interfacial(i,j,k)) then
      restpart = restpart - tracer%ka
      if (restpart<0.0_dp) error=.true.

      Propagated_Quantity(:,i,j,k,next) = &
          Propagated_Quantity (:,i,j,k,next) &
        + restpart * Propagated_Quantity (:,i,j,k,now) &
        + Propagated_Quantity_Adsorbed (:,i,j,k,now) * tracer%kd

      Propagated_Quantity_Adsorbed(:,i,j,k,next) = &
          Propagated_Quantity_Adsorbed(:,i,j,k,now) * (1.0_dp - tracer%kd) &
        + Propagated_Quantity(:,i,j,k,now)*tracer%ka

      vel: do concurrent (l=lbm%lmin+1:lbm%lmax)
          ip = pbc(i+lbm%vel(l)%coo(x) ,x)
          jp = pbc(j+lbm%vel(l)%coo(y) ,y)
          kp = pbc(k+lbm%vel(l)%coo(z) ,z)
        if (.not. is_interfacial (ip,jp,kp)) cycle ! is_interfacial is fluid AND interface
        fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp)) ! 1/2 when tracer has no charge
        scattprop = calc_scattprop( n(i,j,k,l), rho(i,j,k), lbm%vel(l)%a0, lambda_s, fermi)
        restpart = restpart - scattprop
        l_inv = lbm%vel(l)%inv
        scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv), rho(ip,jp,kp), lbm%vel(l_inv)%a0, lambda_s, 1.0_dp-fermi)
        Propagated_Quantity_adsorbed (:,i,j,k,next) = &
            Propagated_Quantity_adsorbed (:,i,j,k,next) &
          + Propagated_Quantity_adsorbed (:,ip,jp,kp,now)*scattprop_p
        u_star = u_star + scattprop* lbm%vel(l)%coo(:)
      end do vel



    else if( .not. is_interfacial(i,j,k)) then
      Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity(:,i,j,k,next) + restpart*Propagated_Quantity(:,i,j,k,now)
    end if

  end do fluidnodes

  if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

  if(modulo(it,(tmax-tmom)/10)==0) print*,it,vacf(x,now),vacf(y,now),vacf(z,now) ! print to user every 1/100 steps. X should be read in input file. To be implemented.

  ! back to the futur: the futur is now, and reinit futur
  Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
  Propagated_Quantity(:,:,:,:,next) = 0.0_dp
  Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
  Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp

  open(unit=99, file='output/vacf.dat', access='append')
  write(99,*) it, vacf(x,now), vacf(y,now), vacf(z,now)

  vacf(:,past) = vacf(:,now)
  vacf(:,now) = 0.0_dp

  if( it>2 .and. all(abs(vacf)<1._dp/(2._dp*lx*ly*lz/tracer%Db)) .and. all(abs(vacf)<1.e-10) ) then
    is_converged = .true.
  else
    is_converged = .false.
  end if


END SUBROUTINE PROPAGATE



! ==============================================================================

REAL(DP) PURE FUNCTION CALC_EXP_DPHI( i, j, k, ip, jp, kp)
  integer(i2b), intent(in) :: i, j, k, ip, jp, kp
  if( tracer%z == 0.0_dp ) then
    calc_exp_dphi = 1.0_dp
  else
    calc_exp_dphi = exp( tracer%z * dphi(i,j,k,ip,jp,kp) )
  end if
END FUNCTION CALC_EXP_DPHI

! ==============================================================================

REAL(DP) PURE FUNCTION DPHI(i,j,k,ip,jp,kp)
  use system, only: phi, elec_slope, lx, ly, lz
  integer(i2b), intent(in) :: i, j, k, ip, jp, kp
  dphi = phi(ip,jp,kp) - phi(i,j,k)
  if      (i==lx .and. ip==1) then
    dphi = dphi + elec_slope(x)*(lx+1)
  else if(i==1 .and. ip==lx) then
    dphi = dphi - elec_slope(x)*(lx+1)
  else if(j==ly .and. jp==1) then
    dphi = dphi + elec_slope(y)*(ly+1)
  else if(j==1 .and. jp==ly) then
    dphi = dphi - elec_slope(y)*(ly+1)
  else if(k==lz .and. kp==1) then
    dphi = dphi + elec_slope(z)*(lz+1)
  else if(k==1 .and. kp==lz) then
    dphi = dphi - elec_slope(z)*(lz+1)
  end if
END FUNCTION DPHI

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_LAMBDA()
  use system, only: kBT
  calc_lambda = 4.0_dp*tracer%Db/kBT
END FUNCTION CALC_LAMBDA

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_LAMBDA_S()
  use system, only: kBT
  calc_lambda_s = 4.0_dp*tracer%Ds/kBT
END FUNCTION CALC_LAMBDA_S

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_SCATTPROP(n,rho,w,lambda,fermi)
  real(dp), intent(in) :: n, rho, w, lambda, fermi
  calc_scattprop = n/rho - w + lambda*w*fermi
END FUNCTION CALC_SCATTPROP

! ==============================================================================

SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY
  if (allocated(Propagated_Quantity)) deallocate(Propagated_Quantity)
  if (allocated(Propagated_Quantity_Adsorbed)) deallocate(Propagated_Quantity_Adsorbed)
END SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY

! ==============================================================================

SUBROUTINE TEST_AND_ALLOCATE_WHAT_IS_NEEDED_FOR_MOMENT_PROPAGATION
  use system, only: lx, ly, lz
  ! Propagated_Quantity is the probability vector of arriving at r at time t
  if(allocated(Propagated_Quantity)) stop 'Propagated quantity should not be allocated in init_moment_propagation'
  allocate(Propagated_Quantity(x:z,lx,ly,lz,now:next), source=0.0_dp)
  if(allocated(Propagated_Quantity_Adsorbed)) stop 'Propagated_Quantity_Adsorbed should not be allocated yet'
  allocate(Propagated_Quantity_Adsorbed(x:z,lx,ly,lz,now:next), source=0.0_dp)
END SUBROUTINE TEST_AND_ALLOCATE_WHAT_IS_NEEDED_FOR_MOMENT_PROPAGATION

! ==============================================================================

LOGICAL PURE FUNCTION TEST(ka)
  real(dp), intent(in) :: ka
  test = .true.
  if( ka < 0.0_dp ) test = .false.
END FUNCTION TEST

! ==============================================================================

!LOGICAL PURE FUNCTION NOT_YET_CONVERGED(t)
!  integer(i2b), intent(in) :: t
!  if(t < (lbound(vacf,2)+2)) then
!    not_yet_converged = .true.
!  else if( vacf(x,lbound(vacf,2)+t) /= vacf(x,lbound(vacf,2)+t-1) &
!      .or. vacf(y,lbound(vacf,2)+t) /= vacf(y,lbound(vacf,2)+t-1) &
!      .or. vacf(z,lbound(vacf,2)+t) /= vacf(z,lbound(vacf,2)+t-1)  ) then
!    not_yet_converged = .true.
!  else
!    not_yet_converged = .false. ! ie is converged
!  end if
!END FUNCTION NOT_YET_CONVERGED

END MODULE MOMENT_PROPAGATION
