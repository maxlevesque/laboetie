MODULE moment_propagation

    USE precision_kinds
    USE constants, only: x, y, z
    USE system, only: tmax, tmom, pbc, supercell, node
    USE mod_lbmodel, only: lbm

    IMPLICIT NONE

    PRIVATE
    PUBLIC :: init, propagate, deallocate_propagated_quantity

    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: Propagated_Quantity, Propagated_Quantity_Adsorbed
    INTEGER(1), PARAMETER :: now=0, next=1, past=-1, tini=past
    REAL(dp), DIMENSION(x:z, past:next) :: vacf
    REAL(dp), DIMENSION(x:z) :: Convergence, vacfOLD
    REAL(dp) :: lambda, lambda_s ! lambda bulk and surface
    TYPE type_tracer
      REAL(dp) :: ka, kd, K, z, Db, Ds !K=ka/kd, z=tracer charge
    END TYPE
    TYPE(type_tracer) :: tracer
    INTEGER, PRIVATE :: lx, ly, lz
    LOGICAL :: considerAdsorption, tracer_is_neutral
    REAL(dp), PARAMETER :: eps=EPSILON(1._dp)

    CONTAINS

    !
    !
    !
SUBROUTINE init( solventDensity )

    USE system, ONLY: phi, fluid, solid, n, node
    use module_input, ONLY: getinput

    IMPLICIT NONE

    real(dp), intent(in) :: solventDensity(:,:,:)
    real(dp) :: boltz_weight, Pstat, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi, sum_of_boltz_weight, rho
    real(dp) :: n_loc(lbm%lmin:lbm%lmax)
    integer :: i, j, k, l, l_inv, ip, jp, kp, i_sum

    tracer%ka = getinput%dp('tracer_ka', defaultValue=0._dp, assert=">=0") ! Adsorption coefficient of the tracer
    tracer%kd = getinput%dp('tracer_kd', defaultValue=0._dp, assert=">=0") ! Desorption coefficient of the tracer

    print*, '+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+'
    print*, 'init called 2'
    print*, '+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+'

    if ( abs(tracer%kd)<= eps ) then ! if kd=0
      tracer%K = 0.0_dp
    else
      tracer%K = tracer%ka / tracer%kd
    end if

    if ( abs(tracer%K) > eps ) then
      considerAdsorption = .true.
    else
      considerAdsorption = .false.
    end if

    tracer%Db = getinput%dp('tracer_Db', defaultValue=0._dp, assert=">=0") ! bulk diffusion coefficient of tracer, i.e. the molecular diffusion coefficient

    tracer%Ds = getinput%dp('tracer_Ds', defaultValue=0._dp, assert=">=0") ! surface diffusion coefficient of tracer

    IF ( ABS(tracer%Ds) > eps ) THEN
        PRINT*,"Tracers you defined have non-zero surface diffusion coefficient"
        PRINT*,"This is not implemented yet"
        ERROR STOP
    END IF

    lambda = calc_lambda(tracer%Db)      ! bulk diffusion
    lambda_s = calc_lambda(tracer%Ds)    ! surface diffusion. is 0 for Ds=0

    tracer%z = getinput%dp('tracer_z', defaultValue=0._dp, assert=">=0") ! tracer's charge
    !IF( ABS(tracer%z)>eps ) THEN
    !    IF( ALLOCATED(phi) ) THEN
    !        PRINT*,"Something is wrong (buuuug) line 79 of module_moment_propagation.f90"
    !        ERROR STOP
    !    END IF
    !    tracer_is_neutral = .FALSE.
    !    PRINT*,"charged tracers are not implemented"
    !    ERROR STOP "line 85 of module_moment_propagation.f90"
    !ELSE
    !    tracer_is_neutral = .TRUE.
    !END IF


    vacf = 0.0_dp ! vacf(x:z, past:next)

    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)
    call test_and_allocate_what_is_needed_for_moment_propagation



    ! the sum of all boltzman weights is the sum over all exp(-z*phi) where node%nature == fluid. 
    ! Note that is_interfacial == fluid + at interface

    IF( tracer_is_neutral ) THEN
        Pstat = COUNT( node%nature == fluid ) + tracer%K*COUNT(node%nature==fluid .AND. node%isinterfacial)
    ELSE
        Pstat = sum(exp(-tracer%z*phi), mask=(node%nature==fluid))&
           +sum(exp(-tracer%z*phi)*tracer%K, mask=node%nature==fluid .and. node%isInterfacial)
    END IF

    sum_of_boltz_weight = 0

    DO CONCURRENT (i=1:lx, j=1:ly, k=1:lz, node(i,j,k)%nature==fluid )

        IF( tracer_is_neutral ) THEN
            boltz_weight = 1._dp/Pstat
        ELSE
            boltz_weight = exp(-tracer%z*phi(i,j,k))/Pstat
        END IF

        n_loc(:) = n(:,i,j,k)
        rho = solventDensity(i,j,k)

        sum_of_boltz_weight = sum_of_boltz_weight + boltz_weight
        i_sum = i_sum +1

        DO CONCURRENT (l=lbm%lmin+1:lbm%lmax)
            ip = pbc (i+lbm%vel(l)%coo(x) ,x)
            jp = pbc (j+lbm%vel(l)%coo(y) ,y)
            kp = pbc (k+lbm%vel(l)%coo(z) ,z)

            if (node(ip,jp,kp)%nature==solid) cycle
            exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
            exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
            fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
            scattprop = calc_scattprop( n_loc(l), rho, lbm%vel(l)%a0, lambda, fermi)
            vacf(:,tini) = vacf(:,tini) + boltz_weight *scattprop *lbm%vel(l)%coo(:)**2

            l_inv = lbm%vel(l)%inv ! comes to r
            scattprop_p = calc_scattprop( &
              n(l_inv,ip,jp,kp), solventDensity(ip,jp,kp), lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)
            Propagated_Quantity(:,i,j,k,tini+1) = Propagated_Quantity(:,i,j,k,tini+1) &
              + exp_min_dphi * scattprop_p * lbm%vel(l_inv)%coo(:) * boltz_weight
        end do
    end do

    PRINT*, 0, REAL(vacf(:,tini),sp)


      if (getinput%log("print_vacf",.true.)) then
          OPEN(99, file='output/vacf.dat')
          WRITE(99,*)'# time t, VACF_x(t), VACF_y(t), VACF_z(t)'
          WRITE(99,*) 0, vacf(:,tini)
      end if

      ! if (considerAdsorption) then
      !   OPEN(100, FILE='output/adsorbed_density.dat')
      !   WRITE(100,*)
      !   WRITE(100,*)"# time ",0
      !   DO i=1,lx; DO j=1,ly; DO k=1,lz;
      !     IF ( node(i,j,k)%isInterfacial .and. node(i,j,k)%nature==fluid ) THEN
      !       WRITE(100,*)i,j,k,SUM(Propagated_Quantity_Adsorbed(:,i,j,k,now))
      !     END IF
      !   END DO; END DO; END DO;
      !   CLOSE(100)
      ! end if

    END SUBROUTINE INIT
    !
    !
    !
    SUBROUTINE PROPAGATE(it, is_converged, solventDensity)

      use system, only: fluid, solid, n, node
      use module_input, only: getinput
      implicit none
      real(dp), intent(in) :: solventDensity(:,:,:)
      integer(i2b), intent(in) :: it
      real(dp) :: fermi, fractionOfParticleRemaining, scattprop, scattprop_p, exp_dphi, rho, n_loc(lbm%lmin:lbm%lmax)
      integer(i2b), parameter :: now=0, next=1, past=-1
      real(dp) :: u_star(x:z), Propagated_Quantity_loc(x:z)
      integer(i2b) :: i, j, k, l, ip, jp, kp, ll, lu, l_inv_loc, ip_all(lbm%lmin:lbm%lmax),&
        jp_all(lbm%lmin:lbm%lmax), kp_all(lbm%lmin:lbm%lmax)
      integer(i2b), save, allocatable :: c(:,:), l_inv(:)
      integer(kind(fluid)), save, allocatable :: nature(:,:,:)
      real(dp), allocatable :: a0(:)
      logical :: error
      logical, save, allocatable :: interfacial(:,:,:)
      logical, intent(out) :: is_converged
      character(1) :: ompvar

      error=.false.

      vacfOLD = 0.0_dp ! Ade : 13/09/17 Initialisation of vacfOLD

      ll = lbm%lmin
      lu = lbm%lmax
      if (.not. allocated(c)) then
        allocate(c(x:z,ll:lu))
        c(x,:) = lbm%vel(:)%coo(x)
        c(y,:) = lbm%vel(:)%coo(y)
        c(z,:) = lbm%vel(:)%coo(z)
      end if
      if (.not. allocated(a0)) allocate(a0(ll:lu), source=lbm%vel(:)%a0)
      if (.not. allocated(l_inv)) allocate(l_inv(ll:lu), source=lbm%vel(:)%inv)
      if (.not. allocated(nature)) allocate (nature(lx,ly,lz), source=node(:,:,:)%nature)
      if (.not. allocated(interfacial)) allocate (interfacial(lx,ly,lz), source=node(:,:,:)%isInterfacial)

      ! ompvar = input_char("openmpover") ! this prepares the code to parallelize over x, y or z slices depending on lb.in. TODO

!$OMP PARALLEL DO PRIVATE(i,j,k,l,ip,jp,kp,fermi,scattprop,l_inv_loc,scattprop_p,n_loc,u_star) &
!$OMP PRIVATE(ip_all,jp_all,kp_all,fractionOfParticleRemaining,Propagated_Quantity_loc) &
!$OMP SHARED(nature,n,lambda,Propagated_Quantity,l_inv,a0,c,ll,lu,lx,ly,lz,considerAdsorption,tracer) &
!$OMP SHARED(Propagated_Quantity_Adsorbed,solventDensity,interfacial,error) &
!$OMP REDUCTION(+:vacf) &
!$OMP DEFAULT(NONE)
      do k=1,lz ! we parallelize over k. If system is 30x30x1 parallelization is useless!
        kp_all(:) = [( pbc( k+c(z,l) ,z) , l=ll,lu)]
        do j=1,ly
          jp_all(:) = [( pbc( j+c(y,l) ,y) , l=ll,lu)]
          do i=1,lx
            if (nature(i,j,k)/=fluid) cycle
            ip_all(:) = [( pbc( i+c(x,l) ,x) , l=ll,lu)]
            u_star(:) = 0.0_dp ! the average velocity at r
            fractionOfParticleRemaining = 1.0_dp ! fraction of particles staying at r; decreases in the loop over neighbours
            n_loc(:) = n(:,i,j,k) ! CTODO CHANGER ORDRE INDICES
            Propagated_Quantity_loc(:) = Propagated_Quantity(:,i,j,k,next)
            do l = ll+1, lu ! ll is velocity=0
              ip = ip_all(l)
              jp = jp_all(l)
              kp = kp_all(l)
              if ( nature(ip,jp,kp) /= fluid ) cycle
              fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp))
              scattprop = calc_scattprop( n_loc(l), solventDensity(i,j,k), a0(l), lambda, fermi) ! scattering probability at r
              fractionOfParticleRemaining = fractionOfParticleRemaining - scattprop ! what is scattered away is not found anymore at r
              u_star(:) = u_star(:) + scattprop*c(:,l)
              l_inv_loc = l_inv(l)
              scattprop_p = calc_scattprop( n(l_inv_loc,ip,jp,kp), solventDensity(ip,jp,kp), a0(l_inv_loc), lambda, 1.0_dp-fermi)
              Propagated_Quantity_loc(:) = Propagated_Quantity_loc(:) + Propagated_Quantity(:,ip,jp,kp,now)*scattprop_p
            end do
            Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity_loc(:)
            vacf(:,now) = vacf(:,now) + Propagated_Quantity(:,i,j,k,now)*u_star(:)
            Convergence = abs(vacf(:,now) - vacfOLD) ! Ade : 13/09/17 I introduced this variable to change the convergence criterion
            vacfOLD = vacf(:,now) ! Ade : Update of vacfOLD

            ! NOW, UPDATE THE PROPAGATED QUANTITIES
            if (   (.not.Interfacial(i,j,k) .and. considerAdsorption) &
              .or. (.not.considerAdsorption) )then
              Propagated_Quantity(:,i,j,k,next) = &
                Propagated_Quantity(:,i,j,k,next) + fractionOfParticleRemaining*Propagated_Quantity(:,i,j,k,now)
            else if ( Interfacial(i,j,k) .and. considerAdsorption ) then
              fractionOfParticleRemaining = fractionOfParticleRemaining - tracer%ka
              Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity (:,i,j,k,next) &
                + fractionOfParticleRemaining * Propagated_Quantity (:,i,j,k,now) &
                + Propagated_Quantity_Adsorbed (:,i,j,k,now) * tracer%kd
              Propagated_Quantity_Adsorbed(:,i,j,k,next) = &
                Propagated_Quantity_Adsorbed(:,i,j,k,now) * (1.0_dp - tracer%kd) &
                + Propagated_Quantity(:,i,j,k,now)*tracer%ka
            end if

            if( fractionOfParticleRemaining < epsilon(1._dp) ) error=.true.

          end do
        end do
      end do
!$OMP END PARALLEL DO


      if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi

      !if(modulo(it,10000)==0) print*,it,REAL(vacf(:,now),sp)
      !if(modulo(it,10000)==0) print*,it,REAL(vacf(:,now),sp)
      if(modulo(it,10000)==0) print*,it,vacf(1,now),' one'
      if(modulo(it,10000)==0) print*,it,REAL(vacf(2,now),sp), ' two'
      if(modulo(it,10000)==0) print*,it,REAL(vacf(3,now),sp), 'three'

      ! back to the futur: the futur is now, and reinit futur
      Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
      Propagated_Quantity(:,:,:,:,next) = 0.0_dp
      if (considerAdsorption) then
        Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
        Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp
      end if

      if (getinput%log("print_vacf",.true.)) write(99,*) it, vacf(:,now)

      !~     OPEN(100, FILE='output/adsorbed_density.dat', ACCESS='append')
      !~         WRITE(100,*)
      !~         WRITE(100,*)"# time ",it
      !~         DO i=1,lx; DO j=1,ly; DO k=1,lz;
      !~             IF ( node(i,j,k)%isInterfacial .and. node(i,j,k)%nature==fluid ) THEN
      !~                 WRITE(100,*)i,j,k,SUM(Propagated_Quantity_Adsorbed(:,i,j,k,now))
      !~             END IF
      !~         END DO; END DO; END DO;
      !~     CLOSE(100)

      !Convergence = Convergence + abs(vacf(:,now) - vacfOLD) ! Ade : 13/09/17 I introduced this variable to change the convergence criterion
      !print*,
      !print*, '************************************************************************'
      !print*, 'Convergence =', Convergence(:)
      !print*, '************************************************************************'
    
      !vacfOLD = vacf(:,now) ! Ade : Update of vacfOLD
      vacf(:,past) = vacf(:,now)
      vacf(:,now) = 0.0_dp

      !IF( it>2 .and. all(abs(vacf)<1._dp/(2._dp*lx*ly*lz/tracer%Db)) .and. all(abs(vacf)<1.e-12) ) then ! TODO MAGIC NUMBER REMOVE THAT SOON
      IF( it>2 .and. all(abs(vacf)<1._dp/(2._dp*lx*ly*lz/tracer%Db)) .and. all(abs(Convergence)<1.e-10) ) then ! TODO MAGIC NUMBER REMOVE THAT SOON
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
    IF( tracer_is_neutral ) THEN
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

END MODULE MOMENT_PROPAGATION
