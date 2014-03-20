MODULE MOMENT_PROPAGATION

    USE precision_kinds
    USE constants, only: x, y, z
    USE system, only: tmax, tmom, pbc, supercell
    USE mod_lbmodel, only: lbm
    
    IMPLICIT NONE
    
    PRIVATE
    PUBLIC :: init, propagate, deallocate_propagated_quantity
    
    REAL(dp), ALLOCATABLE, DIMENSION(:,:,:,:,:) :: Propagated_Quantity, Propagated_Quantity_Adsorbed
    INTEGER(i2b), PARAMETER :: now=0, next=1, past=-1
    REAL(dp), DIMENSION(x:z, past:next) :: vacf
    REAL(dp) :: lambda, lambda_s ! lambda bulk and surface
    TYPE type_tracer
        REAL(dp) :: ka, kd, K, z, Db, Ds !K=ka/kd, z=tracer charge
    END TYPE
    TYPE(type_tracer) :: tracer
    INTEGER(i2b), PRIVATE :: lx,ly,lz
    
    CONTAINS

! ==============================================================================

    SUBROUTINE INIT
    
        USE system, ONLY: phi, fluid, solid, n
        USE input, ONLY: input_dp
        REAL(dp) :: boltz_weight, Pstat, scattprop, scattprop_p, fermi, exp_dphi, exp_min_dphi
        INTEGER(i2b) :: i, j, k, l, l_inv, ip, jp, kp
        
        tracer%ka = input_dp('tracer_ka') ! adsorption 
        tracer%kd = input_dp('tracer_kd') ! desorption
        IF( .NOT. testPositivity(tracer%ka) ) STOP 'I detected tracer%ka to be <0 in module moment_propagation. STOP.'
        IF( .NOT. testPositivity(tracer%kd) ) STOP 'I detected tracer%kd to be <0 in module moment_propagation. STOP.'

        IF(tracer%kd==0.0_dp) THEN
            tracer%K = 0.0_dp
        ELSE
            tracer%K = (tracer%ka)/(tracer%kd)
        END IF
        
        tracer%z = input_dp('tracer_z') ! tracer's charge
        tracer%Db = input_dp('tracer_Db') ! bulk diffusion coefficient of tracer, i.e. the molecular diffusion coefficient
        tracer%Ds = input_dp('tracer_Ds') ! surface diffusion coefficient of tracer
        IF (tracer%Db <= 0.0_dp ) STOP 'tracer_Db as readen in input is invalid'
        IF (tracer%Ds /= 0.0_dp ) STOP "I've found a non-zero Ds (surface diffusion coefficient) in input file. Not implemented yet"
        
        lambda = calc_lambda()
        lambda_s = calc_lambda_s()
        
        vacf = 0.0_dp
        
        lx = supercell%geometry%dimensions%indiceMax(x)
        ly = supercell%geometry%dimensions%indiceMax(y)
        lz = supercell%geometry%dimensions%indiceMax(z)
        CALL test_and_allocate_what_is_needed_for_moment_propagation
        
        ! the sum of all boltzman weights is the sum over all exp(-z*phi) where supercell%node%nature == fluid. Note that is_interfacial == fluid + at interface
        Pstat = sum(exp(-tracer%z*phi), mask=(supercell%node%nature==fluid))&
                +sum(exp(-tracer%z*phi)*tracer%K, mask=supercell%node%nature==fluid .and. supercell%node%isInterfacial)
        
        DO concurrent (i=1:lx, j=1:ly, k=1:lz, supercell%node(i,j,k)%nature==fluid )
            boltz_weight = exp(-tracer%z*phi(i,j,k))/Pstat ! boltz_weight=1/Pstat if tracer%z=0
    
            DO concurrent (l=lbm%lmin+1:lbm%lmax)
                ip = pbc (i+lbm%vel(l)%coo(x) ,x)
                jp = pbc (j+lbm%vel(l)%coo(y) ,y)
                kp = pbc (k+lbm%vel(l)%coo(z) ,z)
                if (supercell%node(ip,jp,kp)%nature==solid) cycle
                exp_dphi = calc_exp_dphi( i, j, k, ip, jp, kp)
                exp_min_dphi = 1.0_dp/exp_dphi ! =1 if tracer%z=0
                fermi = 1.0_dp/(1.0_dp + exp_dphi) ! =0.5 if tracer%z=0
                scattprop = calc_scattprop( n(i,j,k,l), supercell%node(i,j,k)%solventDensity, lbm%vel(l)%a0, lambda, fermi)
                vacf(:,past) = vacf(:,past) + boltz_weight * scattprop * lbm%vel(l)%coo(:)**2
                l_inv = lbm%vel(l)%inv
                scattprop_p = calc_scattprop( &
                    n(ip,jp,kp,l_inv), supercell%node(ip,jp,kp)%solventDensity, lbm%vel(l_inv)%a0, lambda, 1.0_dp-fermi)
                Propagated_Quantity(:,i,j,k,now) = Propagated_Quantity(:,i,j,k,now) &
                            + exp_min_dphi * scattprop_p * lbm%vel(l_inv)%coo(:) * boltz_weight
            END DO
            
            IF(supercell%node(i,j,k)%isInterfacial .and. supercell%node(i,j,k)%nature==fluid) THEN
                Propagated_Quantity_Adsorbed(:,i,j,k,now) = 0.0_dp
            END IF
        end do
        
        PRINT*, 0, vacf(x,past), vacf(y,past), vacf(z,past)
        
        OPEN(99, file='output/vacf.dat')
            WRITE(99,*)'# time t, VACF_x(t), VACF_y(t), VACF_z(t)'
            WRITE(99,*) 0, vacf(x,past), vacf(y,past), vacf(z,past)
        CLOSE(99)
    
        OPEN(100, FILE='output/adsorbed_density.dat')
            WRITE(100,*)
            WRITE(100,*)"# time ",0
            DO i=1,lx; DO j=1,ly; DO k=1,lz;
                IF ( supercell%node(i,j,k)%isInterfacial .and. supercell%node(i,j,k)%nature==fluid ) THEN
                    WRITE(100,*)i,j,k,SUM(Propagated_Quantity_Adsorbed(:,i,j,k,now))
                END IF
            END DO; END DO; END DO;
        CLOSE(100)
    
    END SUBROUTINE INIT

! ==============================================================================

SUBROUTINE PROPAGATE(it, is_converged)
    
    use system, only: fluid, solid, n
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
            scattprop = calc_scattprop( n(i,j,k,l), supercell%node(i,j,k)%solventDensity, lbm%vel(l)%a0, lambda, fermi)
            restpart = restpart - scattprop
            scattprop_p = calc_scattprop( n(ip,jp,kp,lbm%vel(l)%inv), supercell%node(ip,jp,kp)%solventDensity,&
                    lbm%vel((lbm%vel(l)%inv))%a0, lambda, 1.0_dp-fermi)
            Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity(:,i,j,k,next) + Propagated_Quantity(:,ip,jp,kp,now)*scattprop_p
            u_star = u_star + scattprop * lbm%vel(l)%coo(:)
        end do
    
        vacf(:,now) = vacf(:,now) + Propagated_Quantity(:,i,j,k,now)*u_star(:)
    
        if (supercell%node(i,j,k)%isInterfacial .and. supercell%node(i,j,k)%nature==fluid) then
            restpart = restpart - tracer%ka
            if (restpart<0.0_dp) error=.true.
        
            Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity (:,i,j,k,next) &
                + restpart * Propagated_Quantity (:,i,j,k,now) &
                + Propagated_Quantity_Adsorbed (:,i,j,k,now) * tracer%kd
        
            Propagated_Quantity_Adsorbed(:,i,j,k,next) = &
                Propagated_Quantity_Adsorbed(:,i,j,k,now) * (1.0_dp - tracer%kd) &
                + Propagated_Quantity(:,i,j,k,now)*tracer%ka
        
            do concurrent (l=lbm%lmin+1:lbm%lmax)
                ip = pbc(i+lbm%vel(l)%coo(x) ,x)
                jp = pbc(j+lbm%vel(l)%coo(y) ,y)
                kp = pbc(k+lbm%vel(l)%coo(z) ,z)
                if (.not. (supercell%node(ip,jp,kp)%isInterfacial .and. supercell%node(ip,jp,kp)%nature==fluid)) cycle ! is_interfacial is fluid AND interface
                fermi = 1.0_dp/(1.0_dp + calc_exp_dphi(i,j,k,ip,jp,kp)) ! 1/2 when tracer has no charge
                scattprop = calc_scattprop( n(i,j,k,l), supercell%node(i,j,k)%solventDensity, lbm%vel(l)%a0, lambda_s, fermi)
                restpart = restpart - scattprop
                l_inv = lbm%vel(l)%inv
                scattprop_p = calc_scattprop( n(ip,jp,kp,l_inv), supercell%node(ip,jp,kp)%solventDensity, lbm%vel(l_inv)%a0,&
                                                                    lambda_s, 1.0_dp-fermi)
                Propagated_Quantity_adsorbed (:,i,j,k,next) = &
                    Propagated_Quantity_adsorbed (:,i,j,k,next) &
                + Propagated_Quantity_adsorbed (:,ip,jp,kp,now)*scattprop_p
                u_star = u_star + scattprop* lbm%vel(l)%coo(:)
            end do
        
        else if( .not. (supercell%node(i,j,k)%isInterfacial .and. supercell%node(i,j,k)%nature==fluid)) then
            Propagated_Quantity(:,i,j,k,next) = Propagated_Quantity(:,i,j,k,next) + restpart*Propagated_Quantity(:,i,j,k,now)
        end if
    
    end do fluidnodes
    
    if(error) stop 'somewhere restpart is negative' ! TODO one should find a better function for ads and des, as did Benjamin for pi
    
    if(modulo(it,1000)==0) print*,it,vacf(x,now),vacf(y,now),vacf(z,now) ! print to user every 1/100 steps. X should be read in input file. To be implemented.
    
    ! back to the futur: the futur is now, and reinit futur
    Propagated_Quantity(:,:,:,:,now) = Propagated_Quantity(:,:,:,:,next)
    Propagated_Quantity(:,:,:,:,next) = 0.0_dp
    Propagated_Quantity_Adsorbed(:,:,:,:,now) = Propagated_Quantity_Adsorbed(:,:,:,:,next)
    Propagated_Quantity_Adsorbed(:,:,:,:,next) = 0.0_dp

    OPEN(99, FILE='output/vacf.dat', ACCESS='append')
        WRITE(99,*) it, vacf(x,now), vacf(y,now), vacf(z,now)
    
!~     OPEN(100, FILE='output/adsorbed_density.dat', ACCESS='append')
!~         WRITE(100,*)
!~         WRITE(100,*)"# time ",it
!~         DO i=1,lx; DO j=1,ly; DO k=1,lz;
!~             IF ( supercell%node(i,j,k)%isInterfacial .and. supercell%node(i,j,k)%nature==fluid ) THEN
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
    
    
    
    ! ==============================================================================
    
    REAL(DP) PURE FUNCTION CALC_EXP_DPHI( i, j, k, ip, jp, kp)
        INTEGER(i2b), INTENT(IN) :: i, j, k, ip, jp, kp
        IF( tracer%z == 0.0_dp ) THEN
            calc_exp_dphi = 1.0_dp
        ELSE
            calc_exp_dphi = EXP( tracer%z * dphi(i,j,k,ip,jp,kp) )
        END IF
    END FUNCTION CALC_EXP_DPHI
    
    ! ==============================================================================
    
    REAL(DP) PURE FUNCTION DPHI(i,j,k,ip,jp,kp)
    use system, only: phi, elec_slope
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
    
    SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation
        ! Propagated_Quantity is the probability vector of arriving at r at time t
        IF(ALLOCATED(Propagated_Quantity)) STOP 'Propagated quantity should not be allocated in init_moment_propagation'
        ALLOCATE(Propagated_Quantity(x:z,lx,ly,lz,now:next), source=0.0_dp)
        IF(ALLOCATED(Propagated_Quantity_Adsorbed)) STOP 'Propagated_Quantity_Adsorbed should not be allocated yet'
        ALLOCATE(Propagated_Quantity_Adsorbed(x:z,lx,ly,lz,now:next), source=0.0_dp)
    END SUBROUTINE test_and_allocate_what_is_needed_for_moment_propagation
    
    ! ==============================================================================
    
    LOGICAL PURE FUNCTION testPositivity(ka)
        REAL(dp), INTENT(IN) :: ka
        LOGICAL, PARAMETER :: succeeded=.true.
        IF (ka<0.0_dp) THEN
            testPositivity = .NOT.succeeded
        ELSE
            testPositivity = succeeded
        END IF
    END FUNCTION testPositivity
    
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
    
