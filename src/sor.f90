! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! It is now used for the computation of the phi electrostatic potential but it
! is a much more generic method, so it should be (TODO) an object called
! using a poisson solver.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

subroutine sor

    USE precision_kinds, ONLY: dp, i2b ! dp machine specific double precision, i2b simple precision int
    USE system, ONLY: bjl, tot_sol_charge, pbc, anormf0, phi, kbt, c_plus, c_minus, supercell, node, time
    USE constants, ONLY: pi, x, y, z, zerodp
    USE mod_lbmodel, ONLY: lbm
    USE myallocations
    use module_input, only: getinput

    IMPLICIT NONE

    integer(kind=i2b), parameter :: maxiterations = 500000
    real(kind=dp), parameter :: eps = 1.0e-4 ! convergence tolerance
    real(kind=dp), parameter :: omega = 1.4_dp ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
    real(kind=dp) :: factor
    real(kind=dp) :: anorm ! what we want to minimize, ie diff between phi and phi-old
    real(kind=dp) :: anormf, dphi
    integer(kind=i2b) :: iter ! number of iterations to achieve the tolerance
    integer(kind=i2b) :: i, j, k, l, imin, jmin, kmin, imax, jmax, kmax, geometrie
    integer(kind=i2b) :: pmin, qmin, rmin, timeTMP, n1
    integer(kind=i2b) :: ipass,isw,jsw,ksw 
    integer(kind=i2b) :: isw2,jsw2,ksw2
    real(kind=dp) :: phistar, phiold
    real(kind=dp), dimension(:,:,:), allocatable :: phitmp, phiR, phi_old

    open(105, file='output/anorm.dat')
    !open(106, file='output/phiAlongTIME.dat')
    !open(107, file='output/ELLE.dat')
    !open(108, file='output/PHITMP.dat')
    !open(109, file='output/PHIsor.dat')
    !open(110, file='output/PHIsor1.dat')
    
     ! Ade : The piece 7 lines below were at line 61. I moved them here 
    ! if the system wear no charge, the potential is zero.
    if( tot_sol_charge==0.0_dp ) then
        if(.not.allocated(phi)) call allocateReal3D(phi)
        phi = 0.0_dp
        return ! phi has been computed, go on !
    end if
    
    call allocateReal3D( phitmp )
    geometrie = getinput%int('geometryLabel',-1)
    if (.not. allocated(phi_old)) call allocateReal3D(phi_old)
    phi_old = phi

    ! Ade : the lines below are not working
    kmin = supercell%geometry%dimensions%indiceMin(z)
    kmax = supercell%geometry%dimensions%indiceMax(z)
        !do k = kmin, kmax
        !    write(110,*) k, 'phi = ', phi(:,:,k)
        !end do
    ! Ade : end

    !print*, 'timestep = ', timestep
    ! if the system wear no charge, the potential is zero.
    !if( tot_sol_charge==0.0_dp ) then
    !    if(.not.allocated(phi)) call allocateReal3D(phi)
    !    phi = 0.0_dp
    !    return ! phi has been computed, go on !
    !end if

    ! Ade : Block 1
    anormf = 0.0_dp
    anormf = sum(abs(phi)) 
    ! Ade : end Block 1

    ! Ade : Block 2
    if(time==0) then
       anormf = anormf0
    end if
    ! Ade : end Block 2

    ! Ben: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
    factor = 4.0_dp*pi*bjl*kbt/2.0_dp

    ! Ade : Block 3
    ! Ade : FIVE
    convergenceloop: do iter=1, maxiterations 

        anorm = 0.0_dp ! cumulative diff between new and old phi
        ksw = 0
        ksw2 = 0
        ! Ade : 22/03/17 
        imin = supercell%geometry%dimensions%indiceMin(x)
        imax = supercell%geometry%dimensions%indiceMax(x)
        jmin = supercell%geometry%dimensions%indiceMin(y)
        jmax = supercell%geometry%dimensions%indiceMax(y)
        kmin = supercell%geometry%dimensions%indiceMin(z)
        kmax = supercell%geometry%dimensions%indiceMax(z)
        !phitmp = zerodp
          do k = kmin, kmax
            do j = jmin, jmax
                do i = imin, imax 
                !write(108,*) 'i = ', i
            phitmp(i,j,k) = zerodp
            phistar = 0.0_dp
            !write(107,*) '# l       pmin     qmin     rmin     phistar'
            do l= lbm%lmin, lbm%lmax               ! Ade: 22/03/17 the following lines were using the imin, jmin, kmin variables
                                                   ! which were already being used (see lines above)  
                pmin = pbc(i-lbm%vel(l)%coo(x),x)  ! Ade : does fortran count from the first velocity or the second?
                qmin = pbc(j-lbm%vel(l)%coo(y),y)  ! Compare with C-code
                rmin = pbc(k-lbm%vel(l)%coo(z),z)
                phistar = phistar + lbm%vel(l)%a0 * phi(pmin,qmin,rmin)   
                !write(107,*) l, pmin, qmin, rmin, phistar, lbm%vel(l)%a0, phi(pmin,qmin,rmin)
            end do
            phistar = phistar + factor*(c_plus(i,j,k)-c_minus(i,j,k)) ! see PRE64, Horbach: Eq. 17

            phiold = phi(i,j,k)
            phitmp(i,j,k) = omega*phistar +(1.0_dp-omega)*phiold
            !write(108,*) 'phitmp = ', phitmp
            anorm = anorm + abs(phistar-phiold)
                end do
          end do
         end do
         do k=kmin, kmax
            do j= jmin, jmax
              do i= imin, imax 
                phi(i,j,k) = phitmp(i,j,k)
              end do
            end do
         end do

    
        ! ********** OUTPUT ******************
        !do k = kmin, kmax
        !    write(106,*) k, 'phitmp =', phitmp(:,:,k), 'phi = ', phi(:,:,k)
        !end do
        ! ********** end OUTPUT *************

        !if(anorm <= eps*anormf) exit convergenceloop
        ! **************************************************************************************
        ! Ade : new convergence criteria!
        dphi =0.0_dp
        ! if convergence is found, exit SOR
        ! count the number of times the array is not zero
        n1 = count(abs(phi_old)>1.0d-6)

        ! at each node, if phiold is different from 0, calculate the normalized relative differance 
        ! between new and old potential
        if(n1/=0) dphi = sum( abs(  (phi - phi_old)/phi_old ), mask= abs(phi_old)>eps) / real(n1,kind=dp) 
        phi_old = phi
        if(anorm <= eps*anormf .or. (iter>1 .and. dphi<1.0d-6) ) exit convergenceloop
        ! **************************************************************************************


        print*, iter,anorm,eps*anormf

        write(105,*) '----------------------------------------------'
        write(105,*) 'anormf = ', anormf
        write(105,*) '----------------------------------------------'
        
        !do k = kmin, kmax
        !    write(109,*) k, 'phitmp =', SUM(phitmp(:,:,k))
        !end do

        ! change criteria after first step,

        ! tell user if maximum convergence steps is reached, ie if no convergence is found
        if( iter == maxiterations ) stop 'maximum iterations 500 000 reached without convergence in sor'

    end do convergenceloop
    ! Ade : end FIVE
        !do k = kmin, kmax
            !write(109,*) k, (phitmp(1,1,k))
            !write(109,*) k, (phitmp(:,:,k))
        !end do

    print*,'SOR converged in',iter-1,'steps'
    print*,'with anormf0 = ',anormf0

    !anormf0 = sum(abs(phi)) ! Ade : Why do we have this here, whilst it is absent in the C-code????
    close(105)
    !close(106)
    !close(107)
    !close(108)
    !close(109)
    !close(110)

    ! ADE : Maybe add a charge test over here
end subroutine sor
