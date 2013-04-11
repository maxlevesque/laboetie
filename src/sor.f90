! This subroutine is a Successive Over Relaxation (SOR) method implementation.
! It is now used for the computation of the phi electrostatic potential but it
! is a much more generic method, so it should be (TODO) an object called
! using a poisson solver.
! see: Horbach and Frenkel, Phys. Rev. E 64, 061507 (2001)

subroutine sor

  use precision_kinds, only: dp, i2b ! dp machine specific double precision, i2b simple precision int
  use system, only: bjl, sigma, pbc, anormf0, phi, kbt, lx, ly, lz, c_plus, c_minus
  use constants, only: pi, x, y, z
  use mod_lbmodel, only: lbm
  implicit none

  integer(kind=i2b), parameter :: maxiterations = 500000
  real(kind=dp), parameter :: eps = 1.0e-5 ! convergence tolerance
  real(kind=dp), parameter :: omega = 1.4_dp ! the over-ralaxation factor. 1.45 proposed by Horbach & Frenkel, PRE64 Eq.17
  real(kind=dp) :: factor
  real(kind=dp) :: anorm ! what we want to minimize, ie diff between phi and phi-old
  integer(kind=i2b) :: iter ! number of iterations to achieve the tolerance
  integer(kind=i2b) :: i, j, k, l, imin, jmin, kmin
  real(kind=dp) :: phistar, phiold
  real(kind=dp), dimension(lx,ly,lz) :: phitmp
  

  ! if the system wear no charge, the potential is zero.
  if( sigma==0 ) then
    if(.not.allocated(phi)) allocate(phi(lx,ly,lz))
    phi = 0.0_dp
    return ! phi has been computed, go on !
  end if

  ! Ben: corresponds to (4*pi*bjl)*(cs^2 /2) in Eq (17) of PRE64, 061507) ok for cs^2 = 1/2 for D3Q18 only
  factor = 4.0_dp*pi*bjl*kbt/2.0_dp

  convergenceloop: do iter=1, maxiterations

    anorm = 0.0_dp ! cumulative diff between new and old phi

    do i=1, lx
      do j=1, ly
        do k=1, lz

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
