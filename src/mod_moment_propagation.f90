MODULE MOMENT_PROPAGATION
  use precision_kinds, only: dp, i2b
  use constants, only: x, y, z
  use system, only: tmax, tmom
  implicit none
  private
  public :: init,&
            propagate,&
            deallocate_propagated_quantity!, not_yet_converged
  real(dp), allocatable, dimension(:,:,:,:,:) :: Propagated_Quantity
  real(dp), allocatable, dimension(:,:,:,:,:) :: Propagated_Quantity_Adsorbed
  integer(i2b), parameter :: now=0, next=1, past=-1
  real(dp), dimension(x:z, past:next) :: vacf
  real(dp), dimension(x:z), save :: int_vacf
  real(dp) :: lambda
  TYPE TYPE_TRACER
    real(dp) :: ka, kd, K, z, D !K=ka/kd, z=tracer charge
  END TYPE TYPE_TRACER
  type(type_tracer) :: tracer
  contains

! ==============================================================================

SUBROUTINE INIT
  use system, only: phi, z_tracer, NbVel, inside, fluid, solid,&
                     lx, ly, lz, a0, c, vel_inv, plusx, plusy, plusz, rho, n
  use supercell, only: is_interfacial
  use input, only: input_dp
  implicit none
  real(dp) :: boltz_weight, Pstat, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi
  integer(i2b) :: i, j, k, l, l_inv, ip, jp, kp

  tracer%ka = input_dp('tracer_ka')
  tracer%kd = input_dp('tracer_kd')
  if( .not. test(tracer%ka) ) stop 'test(ka) in mod_moment_propagation has detected a problem with tracer%ka'
  if( .not. test(tracer%kd) ) stop 'test(kd) in mod_moment_propagation has detected a problem with tracer%kd'
  if(tracer%kd==0.0_dp) then
    tracer%K = 0.0_dp
  else
    tracer%K = (tracer%ka)/(tracer%kd)
  end if
  tracer%z = input_dp('z_tracer') ! charge of tracer
  if( tracer%z /= z_tracer ) stop 'problem in tracer charge GRRRR'

  lambda = calc_lambda()

  vacf = 0.0_dp

  call test_and_allocate_what_is_needed_for_moment_propagation

  ! the sum of all boltzman weights is the sum over all exp(-z*phi) where inside == fluid. Note that is_interfacial == fluid + at interface
  Pstat = sum(exp(-tracer%z*phi), mask=(inside==fluid))&
         +sum(exp(-tracer%z*phi)*tracer%K, mask=(is_interfacial))

  do concurrent (i=1:lx, j=1:ly, k=1:lz, inside(i,j,k)==fluid )
    boltz_weight = exp(-tracer%z*phi(i,j,k))/Pstat ! boltz_weight=1 if tracer%z=0
    do concurrent (l=2:NbVel, inside(plusx(i+c(x,l)), plusy(j+c(y,l)), plusz(k+c(z,l)))==fluid)
      ip = plusx(i+c(x,l)) ; jp = plusy(j+c(y,l)) ; kp = plusz(k+c(z,l))
      exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
      exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
      fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
      scattprop = calc_scattprop( n(i,j,k,l), rho(i,j,k), a0(l), lambda, fermi)
      vacf(:,past) = vacf(:,past) + boltz_weight * scattprop * c(:,l)**2
      l_inv = vel_inv(l)
      scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv), rho(ip,jp,kp), a0(l_inv), lambda, 1.0_dp-fermi)
      Propagated_Quantity(:,i,j,k,now) = Propagated_Quantity(:,i,j,k,now) &
                 + exp_min_dphi * scattprop_p * c(:,l_inv) * boltz_weight
    end do
    if(is_interfacial(i,j,k)) Propagated_Quantity_Adsorbed(:,i,j,k,now) = 0.0_dp
  end do

  print*, 0, vacf(x,past), vacf(y,past), vacf(z,past)
  open(unit=99, file='output/vacf.dat')
  write(99,*)'# time t, VACF_x(t), VACF_y(t), VACF_z(t), Integrate[vacf_x,t], Integrate[vacf_y,t], Integrate[vacf_z,t]'
  write(99,*) 0, vacf(x,past), vacf(y,past), vacf(z,past), 0.0_dp, 0.0_dp, 0.0_dp ! time, vacf_x, vacf_y, vacf_z, int_vacf_x, int_vacf_y, int_vacf_z
  close(99)
END SUBROUTINE INIT

! ==============================================================================

SUBROUTINE PROPAGATE(it)
  use system, only: lx, ly, lz, inside, fluid, solid,&
                     n, rho, a0, c, vel_inv,&
                     NbVel, plusx, plusy, plusz
  use supercell, only: is_interfacial
!  use parallel, only: Nparallel_threads => Nthread
!$ use OMP_LIB
  implicit none
  integer(kind=i2b), intent(in) :: it
  real(kind=dp) :: lambda, fermi, restpart, scattprop, scattprop_p
  integer(kind=i2b), parameter :: now=0, next=1, past=-1
  real(dp), dimension(3) :: u_star
  integer(kind=i2b) :: i, j, k, l, ip, jp, kp
  logical :: error
  real(dp), allocatable, dimension(:,:,:,:,:) :: prop_quant, prop_quant_ad
  integer(i2b) :: thread, Nthread, kl, ku


  lambda = calc_lambda()

  error=.false.

!!$ print*,OMP_GET_NUM_THREADS(), OMP_GET_NUM_PROCS()


  Nthread = lz!min(lz/4,16)
!$ Nthread = OMP_GET_NUM_PROCS()

!$ CALL OMP_SET_NUM_THREADS(NTHREAD)

!$OMP PARALLEL DEFAULT(FIRSTPRIVATE) &
!$OMP SHARED(Nthread,lx,ly,lz,propagated_quantity,propagated_quantity_adsorbed,n,rho,inside,is_interfacial,error) &
!$OMP PRIVATE(thread,kl,ku,prop_quant,prop_quant_ad,u_star,restpart,scattprop,scattprop_p,ip,jp,kp,i,j,k,l,fermi) &
!$OMP REDUCTION(+:vacf)
!$OMP DO
  do thread = 1, Nthread
    if ( thread == 1 ) then
      kl = 1
      ku = lz/Nthread
    else if ( thread == Nthread) then
      kl = (thread-1)*lz/Nthread + 1
      ku = lz
    else
      kl = (thread-1)*lz/Nthread + 1
      ku = thread*lz/Nthread
    end if
  
    allocate ( prop_quant( x:z, 1:lx, 1:ly, kl-1:ku+1, now:next) )
    allocate ( prop_quant_ad( x:z, 1:lx, 1:ly, kl-1:ku+1, now:next) )
  
    if ( thread == 1 .and. Nthread /= 1) then
      prop_quant(:,:,:, kl-1, :) = propagated_quantity(:,:,:, lz, :) 
      prop_quant(:,:,:, kl:ku+1, :) = propagated_quantity(:,:,:, kl:ku+1, :)
      prop_quant_ad(:,:,:, kl-1, :) = propagated_quantity_adsorbed(:,:,:, lz, :) 
      prop_quant_ad(:,:,:, kl:ku+1, :) = propagated_quantity_adsorbed(:,:,:, kl:ku+1, :)
    elseif ( thread == 1 .and. Nthread == 1) then
      prop_quant(:,:,:, kl-1, :) = propagated_quantity(:,:,:, lz, :) 
      prop_quant(:,:,:, kl:ku,:) = propagated_quantity(:,:,:, kl:ku, :)
      prop_quant(:,:,:, ku+1, :) = propagated_quantity(:,:,:, 1 , :)
      prop_quant_ad(:,:,:, kl-1, :) = propagated_quantity_adsorbed(:,:,:, lz, :) 
      prop_quant_ad(:,:,:, kl:ku, :) = propagated_quantity_adsorbed(:,:,:, kl:ku, :)
      prop_quant_ad(:,:,:, ku+1, :) = propagated_quantity_adsorbed(:,:,:, 1, :)
    else if ( thread == Nthread) then
      prop_quant(:,:,:, kl-1:ku, :) = propagated_quantity(:,:,:, kl-1:ku, :)
      prop_quant(:,:,:, ku+1, :) = propagated_quantity(:,:,:, 1, :)
      prop_quant_ad(:,:,:, kl-1:ku, :) = propagated_quantity_adsorbed(:,:,:, kl-1:ku, :)
      prop_quant_ad(:,:,:, ku+1, :) = propagated_quantity_adsorbed(:,:,:, 1, :)
    else
      prop_quant(:,:,:, kl-1:ku+1, :) = propagated_quantity(:,:,:, kl-1:ku+1, :)
      prop_quant_ad(:,:,:, kl-1:ku+1, :) = propagated_quantity_adsorbed(:,:,:, kl-1:ku+1, :)
    end if
  
    do concurrent ( i = 1:lx, j = 1:ly, k = kl:ku, inside(i,j,k) == fluid )
      u_star = 0.0_dp
      restpart = 1.0_dp
  
      do l = 2, nbvel
        ip = plusx( i + c(x,l) )
        jp = plusy( j + c(y,l) )
        kp = plusz( k + c(z,l) )
  
        if ( inside(ip,jp,kp) /= fluid ) cycle
        fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
        scattprop = calc_scattprop( n(i,j,k,l), rho(i,j,k), a0(l), lambda, fermi)
        restpart = restpart - scattprop
        scattprop_p = calc_scattprop( n(ip,jp,kp,vel_inv(l)), rho(ip,jp,kp), a0(vel_inv(l)), lambda, 1.0_dp-fermi)
        if( thread == 1 .and. plusz( k + c(z,l)) == lx ) then
          kp = 0
        else if ( thread == Nthread .and. plusz( k + c(z,l)) == 1 ) then
          kp = lx+1
        end if
        prop_quant(:,i,j,k,next) = prop_quant(:,i,j,k,next) + prop_quant(:,ip,jp,kp,now)*scattprop_p
        u_star = u_star + scattprop*c(:,l)
      end do
  
      if ( is_interfacial(i,j,k) ) then
        restpart = restpart - tracer%ka
        prop_quant(:,i,j,k,next) = prop_quant(:,i,j,k,next) + restpart*prop_quant(:,i,j,k,now) &
                                            + prop_quant_ad(:,i,j,k,now)*tracer%kd
        prop_quant_ad(:,i,j,k,next) = prop_quant_ad(:,i,j,k,now)*(1.0_dp-tracer%kd) &
                                                    + prop_quant(:,i,j,k,now)*tracer%ka
      else if( .not. is_interfacial(i,j,k)) then
        prop_quant(:,i,j,k,next) = prop_quant(:,i,j,k,next) + restpart*prop_quant(:,i,j,k,now)
      end if
  
      if (restpart<0.0_dp) error=.true.
      vacf(:,now) = vacf(:,now) + prop_quant(:,i,j,k,now)*u_star(:)
    end do

  propagated_quantity( :, :, :, kl:ku, :) = prop_quant(:,:,:,kl:ku,:)
  propagated_quantity_adsorbed( :, :, :, kl:ku, :) = prop_quant_ad(:,:,:,kl:ku,:)
  deallocate ( prop_quant )
  deallocate ( prop_quant_ad )

end do
!$OMP END DO
!$OMP END PARALLEL


  if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

  if(modulo(it,(tmax-tmom)/10)==0) print*,it,vacf(x,now),vacf(y,now),vacf(z,now) ! print to user every 1/100 steps. X should be read in input file. To be implemented.

  ! back to the futur: the futur is now, and reinit futur

!$OMP PARALLEL WORKSHARE
  Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
  Propagated_Quantity(:,:,:,:,next) = 0.0_dp
  Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
  Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp
!$OMP END PARALLEL WORKSHARE

  int_vacf(:) = int_vacf(:)+ 0.5_dp*(vacf(:,now)+vacf(:,past)) ! trapeze integration

  open(unit=99, file='output/vacf.dat', access='append')
  write(99,*) it, vacf(x,now), vacf(y,now), vacf(z,now), int_vacf(x), int_vacf(y), int_vacf(z)

  vacf(:,past) = vacf(:,now)
  vacf(:,now) = 0.0_dp

END SUBROUTINE PROPAGATE

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_EXP_DPHI( i, j, k, ip, jp, kp)
  implicit none
  integer(i2b), intent(in) :: i, j, k, ip, jp, kp
  if( tracer%z == 0.0_dp ) then
    calc_exp_dphi = 1.0_dp
  else if( dphi(i,j,k,ip,jp,kp) == 0.0_dp ) then
    calc_exp_dphi = 1.0_dp
  else
    calc_exp_dphi = exp( tracer%z * dphi(i,j,k,ip,jp,kp) )
  end if
END FUNCTION CALC_EXP_DPHI

! ==============================================================================

REAL(DP) PURE FUNCTION DPHI(i,j,k,ip,jp,kp)
  use system, only: phi, elec_slope_x, elec_slope_y, elec_slope_z, lx, ly, lz
  implicit none
  integer(i2b), intent(in) :: i, j, k, ip, jp, kp

  dphi = phi(ip,jp,kp) - phi(i,j,k)
  if      (i==lx .and. ip==1) then
    dphi = dphi + elec_slope_x*(lx+1)
  else if(i==1 .and. ip==lx) then
    dphi = dphi - elec_slope_x*(lx+1)
  else if(j==ly .and. jp==1) then
    dphi = dphi + elec_slope_y*(ly+1)
  else if(j==1 .and. jp==ly) then
    dphi = dphi - elec_slope_y*(ly+1)
  else if(k==lz .and. kp==1) then
    dphi = dphi + elec_slope_z*(lz+1)
  else if(k==1 .and. kp==lz) then
    dphi = dphi - elec_slope_z*(lz+1)
  end if
END FUNCTION DPHI

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_LAMBDA()
  use system, only: D_tracer, kBT
  implicit none
  calc_lambda = 4.0_dp*D_tracer/kBT
END FUNCTION CALC_LAMBDA

! ==============================================================================

REAL(DP) PURE FUNCTION CALC_SCATTPROP(n,rho,w,lambda,fermi)
  implicit none
  real(dp), intent(in) :: n, rho, w, lambda, fermi
  calc_scattprop = n/rho - w + lambda*w*fermi
END FUNCTION CALC_SCATTPROP

! ==============================================================================

SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY
  implicit none
  if (allocated(Propagated_Quantity)) deallocate(Propagated_Quantity)
  if (allocated(Propagated_Quantity_Adsorbed)) deallocate(Propagated_Quantity_Adsorbed)
END SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY

! ==============================================================================

SUBROUTINE TEST_AND_ALLOCATE_WHAT_IS_NEEDED_FOR_MOMENT_PROPAGATION
  use system, only: lx, ly, lz
  implicit none
  ! Propagated_Quantity is the probability vector of arriving at r at time t
  if(allocated(Propagated_Quantity)) stop 'Propagated quantity should not be allocated in init_moment_propagation'
  allocate(Propagated_Quantity(x:z,lx,ly,lz,now:next), source=0.0_dp)
  if(allocated(Propagated_Quantity_Adsorbed)) stop 'Propagated_Quantity_Adsorbed should not be allocated yet'
  allocate(Propagated_Quantity_Adsorbed(x:z,lx,ly,lz,now:next), source=0.0_dp)
END SUBROUTINE TEST_AND_ALLOCATE_WHAT_IS_NEEDED_FOR_MOMENT_PROPAGATION

! ==============================================================================

LOGICAL PURE FUNCTION TEST(ka)
  implicit none
  real(dp), intent(in) :: ka
  test = .true.
  if( ka < 0.0_dp ) test = .false.
END FUNCTION TEST

! ==============================================================================

!LOGICAL PURE FUNCTION NOT_YET_CONVERGED(t)
!  implicit none
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
