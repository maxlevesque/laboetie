! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! It is now used for the computation of the phi electrostatic potential but it
! is a much more generic method, so it should be (TODO) an object called
! using a poisson solver.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

subroutine sor(timestep)

    USE precision_kinds, ONLY: dp, i2b ! dp machine specific double precision, i2b simple precision int
    USE system, ONLY: bjl, sigma, pbc, anormf0, phi, kbt, c_plus, c_minus, supercell, node
    USE constants, ONLY: pi, x, y, z, zerodp
    USE mod_lbmodel, ONLY: lbm
    USE myallocations
    use module_input, only: getinput

    IMPLICIT NONE

    integer(kind=i2b), parameter :: maxiterations = 3!500000
    integer, intent(in) :: timestep
    real(kind=dp), parameter :: eps = 1.0e-5 ! convergence tolerance
    real(kind=dp), parameter :: omega = 1.4_dp ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
    real(kind=dp) :: factor
    real(kind=dp) :: anorm ! what we want to minimize, ie diff between phi and phi-old
    real(kind=dp) :: anormf
    integer(kind=i2b) :: iter ! number of iterations to achieve the tolerance
    integer(kind=i2b) :: i, j, k, l, imin, jmin, kmin, imax, jmax, kmax, geometrie
    integer(kind=i2b) :: pmin, qmin, rmin
    real(kind=dp) :: phistar, phiold
    real(kind=dp), dimension(:,:,:), allocatable :: phitmp, phiR

    open(105, file='output/anorm.dat')
    open(106, file='output/phiAlongTIME.dat')
    call allocateReal3D( phitmp )
    geometrie = getinput%int('geometryLabel',-1)

    !print*, 'timestep = ', timestep
    ! if the system wear no charge, the potential is zero.
    if( sigma==0.0_dp ) then
        if(.not.allocated(phi)) call allocateReal3D(phi)
        phi = 0.0_dp
        return ! phi has been computed, go on !
    end if

    
    anormf = 0.0_dp
    ! Ade : 22/03/17 
    anormf = sum(abs(phi)) ! Ade : Why do we have this here, whilst it is absent in the C-code????
    !do i = imin, imax
    !    do j = jmin, jmax
    !        do k = kmin, kmax
    !            anormf = anormf + abs(phi(i,j,k))
    !        end do
    !    end do
    !end do
    ! Ade : 22/03/17
    !print*, 'anormf0 = ', anormf0
    !print*, 'timestep = ', timestep
    if(timestep==0) then
       anormf = anormf0
    print*, 'anorm = ', anorm
    end if
    print*, '-------------------------------------'
    print*, 'anorm = ', anorm
    ! Ade : 22/03/17 

    ! Ben: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
    factor = 4.0_dp*pi*bjl*kbt/2.0_dp

    convergenceloop: do iter=1, maxiterations

        anorm = 0.0_dp ! cumulative diff between new and old phi
        ! Ade : 22/03/17 
        !do i = imin, imax
        imin = supercell%geometry%dimensions%indiceMin(x)
        imax = supercell%geometry%dimensions%indiceMax(x)
        jmin = supercell%geometry%dimensions%indiceMin(y)
        jmax = supercell%geometry%dimensions%indiceMax(y)
        kmin = supercell%geometry%dimensions%indiceMin(z)
        kmax = supercell%geometry%dimensions%indiceMax(z)
        phitmp = zerodp
        do i = imin, imax
            do j = jmin, jmax
                do k = kmin, kmax
            phistar = 0.0_dp
            do l= lbm%lmin, lbm%lmax               ! Ade: 22/03/17 the following lines were using the imin, jmin, kmin variables
                                                   ! which were already being used (see lines above)  
                pmin = pbc(i-lbm%vel(l)%coo(x),x)  ! Ade : does fortran count from the first velocity or the second?
                qmin = pbc(j-lbm%vel(l)%coo(y),y)  ! Compare with C-code
                rmin = pbc(k-lbm%vel(l)%coo(z),z)
                phistar = phistar + lbm%vel(l)%a0 * phi(pmin,qmin,rmin)   ! Ade : 20/03/17 this phi is not initialised anywhere
            end do
            phistar = phistar + factor*(c_plus(i,j,k)-c_minus(i,j,k)) ! see PRE64, Horbach: Eq. 17

            phiold = phi(i,j,k)
            phitmp(i,j,k) = omega*phistar +(1.0_dp-omega)*phiold
            anorm = anorm + abs(phistar-phiold)
            end do
          end do
        end do
         
        do k = kmin, kmax
            write(106,*) k, 'phitmp =', phitmp(:,:,k)
        end do
        ! replace phi by the newly calculated phitmp
        phi = phitmp

        ! inform user every 1000 steps
    !    if(modulo(iter,10000)==0) then
    !      print*,'SOR iter ',iter,' convergence at ',1-anorm/eps*anormf0,sum(phi)
    !    end if

        ! if convergence is found, exit SOR
        if(anorm <= eps*anormf) exit convergenceloop

        write(105,*) '----------------------------------------------'
        write(105,*) 'anormf = ', anormf
        write(105,*) '----------------------------------------------'

        ! change criteria after first step,

        ! tell user if maximum convergence steps is reached, ie if no convergence is found
        if( iter == maxiterations ) stop 'maximum iterations 500 000 reached without convergence in sor'

    end do convergenceloop

    print*,'SOR converged in',iter-1,'steps'
    print*,'with anormf0 = ',anormf0

    !anormf0 = sum(abs(phi)) ! Ade : Why do we have this here, whilst it is absent in the C-code????
    close(105)
    close(106)
end subroutine sor
