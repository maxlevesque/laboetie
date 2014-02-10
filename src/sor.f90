! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! It is now used for the computation of the phi electrostatic potential but it
! is a much more generic method, so it should be (TODO) an object called
! using a poisson solver.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

subroutine sor

    USE precision_kinds, ONLY: dp, i2b ! dp machine specific double precision, i2b simple precision int
    USE system, ONLY: bjl, sigma, pbc, anormf0, phi, kbt, c_plus, c_minus, supercell
    USE constants, ONLY: pi, x, y, z
    USE mod_lbmodel, ONLY: lbm
    USE myallocations
    
    IMPLICIT NONE

    integer(kind=i2b), parameter :: maxiterations = 500000
    real(kind=dp), parameter :: eps = 1.0e-5 ! convergence tolerance
    real(kind=dp), parameter :: omega = 1.4_dp ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
    real(kind=dp) :: factor
    real(kind=dp) :: anorm ! what we want to minimize, ie diff between phi and phi-old
    integer(kind=i2b) :: iter ! number of iterations to achieve the tolerance
    integer(kind=i2b) :: i, j, k, l, imin, jmin, kmin, imax, jmax, kmax
    real(kind=dp) :: phistar, phiold
    real(kind=dp), dimension(:,:,:), allocatable :: phitmp
    call allocateReal3D( phitmp )
    
    ! if the system wear no charge, the potential is zero.
    if( sigma==0.0_dp ) then
        if(.not.allocated(phi)) call allocateReal3D(phi)
        phi = 0.0_dp
        return ! phi has been computed, go on !
    end if
    
    ! Ben: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
    factor = 4.0_dp*pi*bjl*kbt/2.0_dp
    
    convergenceloop: do iter=1, maxiterations
    
        anorm = 0.0_dp ! cumulative diff between new and old phi
    
        imin = supercell%geometry%dimensions%indiceMin(x)
        imax = supercell%geometry%dimensions%indiceMax(x)
        jmin = supercell%geometry%dimensions%indiceMin(y)
        jmax = supercell%geometry%dimensions%indiceMax(y)
        kmin = supercell%geometry%dimensions%indiceMin(z)
        kmax = supercell%geometry%dimensions%indiceMax(z)
        do i = imin, imax
            do j = jmin, jmax
                do k = kmin, kmax
    
            phistar = 0.0_dp
            do l= lbm%lmin, lbm%lmax
                imin = pbc(i-lbm%vel(l)%coo(x),x)
                jmin = pbc(j-lbm%vel(l)%coo(y),y)
                kmin = pbc(k-lbm%vel(l)%coo(z),z)
                phistar = phistar + lbm%vel(l)%a0 * phi(imin,jmin,kmin)
            end do
            phistar = phistar +factor*(c_plus(i,j,k)-c_minus(i,j,k)) ! see PRE64, Horbach: Eq. 17
    
            phiold = phi(i,j,k)
            phitmp(i,j,k) = omega*phistar +(1.0_dp-omega)*phiold
            anorm = anorm + abs(phistar-phiold)
    
            end do
        end do
        end do
    
        ! replace phi by the newly calculated phitmp
        phi = phitmp
    
        ! inform user every 1000 steps
    !    if(modulo(iter,10000)==0) then
    !      print*,'SOR iter ',iter,' convergence at ',1-anorm/eps*anormf0,sum(phi)
    !    end if
    
        ! if convergence is found, exit SOR
        if(anorm <= eps*anormf0) exit convergenceloop
    
        ! change criteria after first step,
    
        ! tell user if maximum convergence steps is reached, ie if no convergence is found
        if( iter == maxiterations ) stop 'maximum iterations 500 000 reached without convergence in sor'
    
    end do convergenceloop
    
    print*,'SOR converged in',iter-1,'steps'
    print*,'with anormf0 = ',anormf0
    
    anormf0 = sum(abs(phi))
end subroutine sor
