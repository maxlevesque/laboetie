MODULE MOMENT_PROPAGATION

  USE precision_kinds
  USE constants, only: x, y, z
  USE system, only: tmax, tmom, pbc, supercell, node
  USE mod_lbmodel, only: lbm

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: init, propagate, deallocate_propagated_quantity

  REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: Propagated_Quantity, Propagated_Quantity_Adsorbed
  INTEGER(i2b), PARAMETER :: now=0, next=1, past=-1, tini=past
  REAL(dp), DIMENSION(x:z, past:next) :: vacf
  REAL(dp) :: lambda, lambda_s ! lambda bulk and surface
  TYPE type_tracer
    REAL(dp) :: ka, kd, K, z, Db, Ds !K=ka/kd, z=tracer charge
  END TYPE
  TYPE(type_tracer) :: tracer
  INTEGER(i2b), PRIVATE :: lx,ly,lz
  logical :: considerAdsorption

  CONTAINS

    ! ==============================================================================

    SUBROUTINE INIT

      use system, ONLY: phi, fluid, solid, n, node
      use input, ONLY: input_dp
      implicit none
      real(dp) :: boltz_weight, Pstat, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi, sum_of_boltz_weight, rho
      real(dp) :: pop(lbm%lmin:lbm%lmax)
      integer(i2b) :: i, j, k, l, l_inv, ip, jp, kp, i_sum

      tracer%ka = input_dp('tracer_ka') ! adsorption
      tracer%kd = input_dp('tracer_kd') ! desorption
      if (tracer%ka < epsilon(1._dp)) stop 'I detected tracer%ka to be <0 in module moment_propagation. STOP.'
      if (tracer%kd < epsilon(1._dp)) stop 'I detected tracer%kd to be <0 in module moment_propagation. STOP.'

      if ( abs(tracer%kd)<=epsilon(1._dp) ) then ! if kd=0
        tracer%K = 0.0_dp ! don't do the division by 0
      else
        tracer%K = tracer%ka / tracer%kd
      end if

      if ( abs(tracer%K)>epsilon(1.0_dp) ) then
        considerAdsorption = .true.
      else
        considerAdsorption = .false.
      end if

      tracer%z = input_dp('tracer_z') ! tracer's charge
      tracer%Db = input_dp('tracer_Db') ! bulk diffusion coefficient of tracer, i.e. the molecular diffusion coefficient
      tracer%Ds = input_dp('tracer_Ds') ! surface diffusion coefficient of tracer
      IF (tracer%Db <= 0.0_dp ) STOP 'tracer_Db as readen in input is invalid'
      IF (tracer%Ds /= 0.0_dp ) STOP "I've found a non-zero Ds (surface diffusion coefficient) in input file. Not implemented yet"

      lambda = calc_lambda(tracer%Db)      ! bulk diffusion
      lambda_s = calc_lambda(tracer%Ds)  ! surface diffusion. is 0 for Ds=0

      vacf = 0.0_dp ! vacf(x:z, past:next)

      lx = supercell%geometry%dimensions%indiceMax(x)
      ly = supercell%geometry%dimensions%indiceMax(y)
      lz = supercell%geometry%dimensions%indiceMax(z)
      call test_and_allocate_what_is_needed_for_moment_propagation

      ! the sum of all boltzman weights is the sum over all exp(-z*phi) where node%nature == fluid. Note that is_interfacial == fluid + at interface
      Pstat = sum(exp(-tracer%z*phi), mask=(node%nature==fluid))&
             +sum(exp(-tracer%z*phi)*tracer%K, mask=node%nature==fluid .and. node%isInterfacial)

      do concurrent (i=1:lx, j=1:ly, k=1:lz, node(i,j,k)%nature==fluid )
        boltz_weight = exp(-tracer%z*phi(i,j,k))/Pstat
        pop = n(i,j,k,:)
        rho = node(i,j,k)%solventDensity

        sum_of_boltz_weight = sum_of_boltz_weight + boltz_weight
        i_sum = i_sum +1

        do concurrent (l=lbm%lmin+1:lbm%lmax)
          ip = pbc (i+lbm%vel(l)%coo(x) ,x)
          jp = pbc (j+lbm%vel(l)%coo(y) ,y)
          kp = pbc (k+lbm%vel(l)%coo(z) ,z)

          if (node(ip,jp,kp)%nature==solid) cycle
          exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
          exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
          fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
          scattprop = calc_scattprop( pop(l), rho, lbm%vel(l)%a0, lambda, fermi)
          vacf(:,tini) = vacf(:,tini) + boltz_weight *scattprop *lbm%vel(l)%coo(:)**2

          l_inv = lbm%vel(l)%inv ! comes to r
          scattprop_p = calc_scattprop( &
            n(ip,jp,kp,l_inv), node(ip,jp,kp)%solventDensity, lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)
          Propagated_Quantity(:,i,j,k,tini+1) = Propagated_Quantity(:,i,j,k,tini+1) &
            + exp_min_dphi * scattprop_p * lbm%vel(l_inv)%coo(:) * boltz_weight
        end do
      end do

      PRINT*, 0, REAL(vacf(:,tini),sp)

      OPEN(99, file='output/vacf.dat')
      WRITE(99,*)'# time t, VACF_x(t), VACF_y(t), VACF_z(t)'
      WRITE(99,*) 0, vacf(x,tini), vacf(y,tini), vacf(z,tini)

      if (considerAdsorption) then
        OPEN(100, FILE='output/adsorbed_density.dat')
        WRITE(100,*)
        WRITE(100,*)"# time ",0
        DO i=1,lx; DO j=1,ly; DO k=1,lz;
          IF ( node(i,j,k)%isInterfacial .and. node(i,j,k)%nature==fluid ) THEN
            WRITE(100,*)i,j,k,SUM(Propagated_Quantity_Adsorbed(:,i,j,k,now))
          END IF
        END DO; END DO; END DO;
        CLOSE(100)
      end if

    END SUBROUTINE INIT

    ! ==============================================================================

    SUBROUTINE PROPAGATE(it, is_converged)

      use system, only: fluid, solid, n, node
      implicit none
      integer(i2b), intent(in) :: it
      real(dp) :: fermi, restpart, scattprop, scattprop_p, exp_dphi, rho, pop(lbm%lmin:lbm%lmax)
      integer(i2b), parameter :: now=0, next=1, past=-1
      real(dp), dimension(3) :: u_star
      integer(i2b) :: i, j, k, l, l_inv, ip, jp, kp, nature
      logical :: error, interfacial
      logical, intent(out) :: is_converged

      error=.false.

      do concurrent (i=1:lx, j=1:ly, k=1:lz, node(i,j,k)%nature==fluid )
        u_star = 0.0_dp ! the average velocity at r
        restpart = 1.0_dp ! fraction of particles staying at r; decreases in the loop over neighbours
        pop = n(i,j,k,:)
        rho = node(i,j,k)%solventDensity
        nature = node(i,j,k)%nature
        interfacial = node(i,j,k)%isInterfacial

        do l = lbm%lmin+1, lbm%lmax
          ip = pbc(i+lbm%vel(l)%coo(x) ,x)
          jp = pbc(j+lbm%vel(l)%coo(y) ,y)
          kp = pbc(k+lbm%vel(l)%coo(z) ,z)
          if ( node(ip,jp,kp)%nature /= fluid ) cycle
          fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
          scattprop = calc_scattprop( pop(l), rho, lbm%vel(l)%a0, lambda, fermi) ! scattering probability at r
          restpart = restpart - scattprop ! what is scattered away is not found anymore at r
          u_star = u_star + scattprop * lbm%vel(l)%coo(:)
          l_inv = lbm%vel(l)%inv
          scattprop_p = calc_scattprop( n(ip,jp,kp,lbm%vel(l)%inv), node(ip,jp,kp)%solventDensity,&
            lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)
          Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity(:,i,j,k,next) + Propagated_Quantity(:,ip,jp,kp,now)*scattprop_p
        end do

        vacf(:,now) = vacf(:,now) + Propagated_Quantity(:,i,j,k,now)*u_star(:)

        ! NOW, UPDATE THE PROPAGATED QUANTITIES

        if ( nature==fluid .and. .not.interfacial ) then
          Propagated_Quantity(:,i,j,k,next) = &
            Propagated_Quantity(:,i,j,k,next) + restpart*Propagated_Quantity(:,i,j,k,now)

        else if ( nature==fluid .and. interfacial .and. considerAdsorption ) then
          restpart = restpart - tracer%ka ! ICI JE METTRAI restpart*(1-ka) pour être plus physique, mais c'est peut être ça qui fait la discontinuité de JMV
          if (restpart<0.0_dp) error=.true.

          Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity (:,i,j,k,next) &
            + restpart * Propagated_Quantity (:,i,j,k,now) &
            + Propagated_Quantity_Adsorbed (:,i,j,k,now) * tracer%kd

          Propagated_Quantity_Adsorbed(:,i,j,k,next) = &
            Propagated_Quantity_Adsorbed(:,i,j,k,now) * (1.0_dp - tracer%kd) &
            + Propagated_Quantity(:,i,j,k,now)*tracer%ka

          ! do concurrent (l=lbm%lmin+1:lbm%lmax)
          !   ip = pbc(i+lbm%vel(l)%coo(x) ,x)
          !   jp = pbc(j+lbm%vel(l)%coo(y) ,y)
          !   kp = pbc(k+lbm%vel(l)%coo(z) ,z)
          !   if (.not. (node(ip,jp,kp)%isInterfacial .and. node(ip,jp,kp)%nature==fluid)) cycle ! is_interfacial is fluid AND interface
          !   fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp)) ! 1/2 when tracer has no charge
          !   scattprop = calc_scattprop( n(i,j,k,l), node(i,j,k)%solventDensity, lbm%vel(l)%a0, lambda_s, fermi)
          !   restpart = restpart - scattprop
          !   l_inv = lbm%vel(l)%inv
          !   scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv), node(ip,jp,kp)%solventDensity, lbm%vel(l_inv)%a0,&
          !     lambda_s, 1.0_dp-fermi)
          !   Propagated_Quantity_adsorbed (:,i,j,k,next) = &
          !     Propagated_Quantity_adsorbed (:,i,j,k,next) + Propagated_Quantity_adsorbed (:,ip,jp,kp,now)*scattprop_p
          !   u_star = u_star + scattprop* lbm%vel(l)%coo(:)
          ! end do

        end if
      end do

      if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

      if(modulo(it,10000)==0) print*,it,REAL(vacf(:,now),sp)

      ! back to the futur: the futur is now, and reinit futur
      Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
      Propagated_Quantity(:,:,:,:,next) = 0.0_dp
      if (considerAdsorption) then
        Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
        Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp
      end if

      WRITE(99,*) it, vacf(x,now), vacf(y,now), vacf(z,now)

      !~     OPEN(100, FILE='output/adsorbed_density.dat', ACCESS='append')
      !~         WRITE(100,*)
      !~         WRITE(100,*)"# time ",it
      !~         DO i=1,lx; DO j=1,ly; DO k=1,lz;
      !~             IF ( node(i,j,k)%isInterfacial .and. node(i,j,k)%nature==fluid ) THEN
      !~                 WRITE(100,*)i,j,k,SUM(Propagated_Quantity_Adsorbed(:,i,j,k,now))
      !~             END IF
      !~         END DO; END DO; END DO;
      !~     CLOSE(100)

      vacf(:,past) = vacf(:,now)
      vacf(:,now) = 0.0_dp

      IF( it>2 .and. all(abs(vacf)<1._dp/(2._dp*lx*ly*lz/tracer%Db)) .and. all(abs(vacf)<1.e-12) ) then ! TODO MAGIC NUMBER REMOVE THAT SOON
        is_converged = .true.
      ELSE
        is_converged = .false.
      END IF

    END SUBROUTINE PROPAGATE



!==============================================================================

  PURE FUNCTION CALC_EXP_DPHI( i, j, k, ip, jp, kp)
    implicit none
    REAL(DP) :: CALC_EXP_DPHI
    INTEGER(i2b), INTENT(IN) :: i, j, k, ip, jp, kp
    IF( abs(tracer%z) <= epsilon(1._dp) ) THEN
      calc_exp_dphi = 1.0_dp
    ELSE
      calc_exp_dphi = EXP( tracer%z * dphi(i,j,k,ip,jp,kp) )
    END IF
  END FUNCTION CALC_EXP_DPHI

  ! ==============================================================================

  PURE FUNCTION DPHI(i,j,k,ip,jp,kp)
    use system, only: phi, elec_slope
    implicit none
    REAL(DP) :: dphi
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

  PURE FUNCTION CALC_LAMBDA(d)
    use system, only: kBT
    implicit none
    real(dp) :: calc_lambda
    real(dp), intent(in) :: d ! diffusion coefficient to be used to compute lambda. See review by Ladd.
    calc_lambda = 4.0_dp*d/kBT
  END FUNCTION CALC_LAMBDA

! ==============================================================================

  pure function calc_scattprop (n,rho,w,lambda,fermi)
    implicit none
    real(dp) :: calc_scattprop
    real(dp), intent(in) :: n, rho, w, lambda, fermi
    calc_scattprop = n/rho - w + lambda*w*fermi
  end function calc_scattprop

! ==============================================================================

    SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY
      if (allocated(Propagated_Quantity)) deallocate(Propagated_Quantity)
      if (allocated(Propagated_Quantity_Adsorbed)) deallocate(Propagated_Quantity_Adsorbed)
    END SUBROUTINE DEALLOCATE_PROPAGATED_QUANTITY

! ==============================================================================

    SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation
      implicit none
      ! Propagated_Quantity is the probability vector of arriving at r at time t
      IF(ALLOCATED(Propagated_Quantity)) STOP 'Propagated quantity should not be allocated in init_moment_propagation'
      ALLOCATE(Propagated_Quantity(x:z,lx,ly,lz,now:next), source=0.0_dp)
      if (considerAdsorption) then
        IF(ALLOCATED(Propagated_Quantity_Adsorbed)) STOP 'Propagated_Quantity_Adsorbed should not be allocated yet'
        ALLOCATE(Propagated_Quantity_Adsorbed(x:z,lx,ly,lz,now:next), source=0.0_dp)
      end if
    END SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation

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
