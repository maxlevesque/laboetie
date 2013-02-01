! Here is a smolu with diffusion + electrostatic part only
! The objective is to end up in a situation where equilibrium
! distribution of solutes (c_plus and c_minus) is found,
! ie in which the total flux of solutes is zero.

SUBROUTINE JUST_EQ_SMOLU
  use precision_kinds, only: i2b, dp
  use system, only: D_plus, D_minus, &
                     inside, lx, ly, lz, c, fluid, solid, phi, delta, c_plus, c_minus,&
                     sigma, NbVel, plusx, plusy, plusz, solute
  use constants, only: x, y, z
  use charges, only: my_lattice_nodes
  implicit none

  integer(i2b) :: iter, max_iter
  real(dp), dimension(lx,ly,lz) :: flux_site_plus, flux_site_minus
  integer(i2b) :: i, j, k, ip, jp, kp, l ! dummy
  real(dp) :: exp_dphi, exp_min_dphi ! Exp[phi site1 - phi site 2] and 1/Exp
  real(dp) :: flux_link_plus, flux_link_minus
  real(dp) :: tot_diff_plus, tot_diff_minus ! total flux of + and - solutes init at high values
  real(dp) :: eD_plus, eD_minus ! effective D_plus and D_minus
  real(dp), parameter :: convergence_criteria = 1.e-10
  integer(i2b), save :: move = 1


  ! init
  tot_diff_plus = huge(1.0_dp)
  tot_diff_minus = huge(1.0_dp)

  ! diffusion coefficients of solutes are read in input file.
  ! if we're in the first steps of equilibration, it speeds up convergence
  ! to use a smaller diffusion coefficient.
  ! i think there is a bug in C code here,
  ! especially considering the high amount of magic numbers and +- convention for time.
  call effective_diffusion(move)

  ! init convergence iterations
  iter = 0

  do while( (tot_diff_minus + tot_diff_plus) > convergence_criteria .and. iter < max_iter )

    ! flux_site_plus or minus is what makes the system evoluate.
    flux_site_plus = 0.0_dp
    flux_site_minus = 0.0_dp

    ! for all sites
    do concurrent( i=1:lx, j=1:ly, k=1:lz, inside(i,j,k)==fluid )
      ! and all neighbours of this site
      do l= 2, NbVel, 2 ! at once flux in both directions ! l=1 corresponds to no velocity ie flux toward itself ie delta(l)=0

        ! periodic boundary conditions to neighbours
        ip= plusx( i+ c(x,l))
        jp= plusy( j+ c(y,l))
        kp= plusz( k+ c(z,l))

        ! continue for fluid-fluid flux only
        if( inside(ip,jp,kp) == fluid ) then
          ! compute the difference in potential between sites i,j,k and ip,jp,kp
          exp_dphi = exp( phi(ip,jp,kp) - phi(i,j,k) ) ! be carefull to sign
          exp_min_dphi = 1.0_dp / exp_dphi ! dummy

          ! flux due to electrostatic and density gradients inside link i,j,k <-> ip,jp,kp
          flux_link_plus  = 0.5_dp * (1.0_dp+ exp_min_dphi)*( c_plus(ip,jp,kp)  * exp_dphi     - c_plus (i,j,k) )&
                            * eD_plus  / delta(l)
          flux_link_minus = 0.5_dp * (1.0_dp+ exp_dphi    )*( c_minus(ip,jp,kp) * exp_min_dphi - c_minus(i,j,k) )&
                            * eD_minus / delta(l)

          ! update flux on each site accordingly to flux inside link involving site
          flux_site_plus(i,j,k) = flux_site_plus(i,j,k) + flux_link_plus
          flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus
          flux_site_plus(ip,jp,kp) = flux_site_plus(ip,jp,kp) - flux_link_plus
          flux_site_minus(ip,jp,kp) = flux_site_minus(ip,jp,kp) - flux_link_minus
        end if
      end do
    end do

    ! no concentration should be lost, just transfered, ie sum of flux over all sites should be 0.
    if( modulo(move,100)==0 ) then
      if( abs(sum(flux_site_plus)) > 1.e-12 .or. abs(sum(flux_site_minus)) > 1.e-12) then
        print*,'abs(sum(flux_site_plus))=',abs(sum(flux_site_plus))
        print*,'abs(sum(flux_site_minus))=',abs(sum(flux_site_minus))
        stop 'the sum of all flux does not add up. problem in just_eq_smolu.f90'
      end if
    end if

    ! update concentrations (smolushowski part)
    where(inside==fluid)
      c_plus = c_plus + flux_site_plus
      c_minus = c_minus + flux_site_minus
    end where

    ! compute the total flux in this equilibration step one wants to minimize.
    if( sigma/=0 .and. eD_plus/=0.0_dp .and. eD_minus/=0.0_dp ) then
      ! the sum of all flux
      tot_diff_plus  = sqrt(sum(flux_site_plus**2,mask=(inside==fluid)))  / my_lattice_nodes%fluid &
             / (0.5_dp*solute%density*eD_plus ) / sigma ! 1st denominator is the number of fluid nodes)
      tot_diff_minus = sqrt(sum(flux_site_minus**2,mask=(inside==fluid))) / my_lattice_nodes%fluid &
             / (0.5_dp*solute%density*eD_minus) / sigma ! norm2 is the Fortran intrinsic for euclidean norm
      print*,'DIFF_plus =',tot_diff_plus,' DIFF_minus =',tot_diff_minus
    end if

    ! increment iteration
    iter = iter +1

  end do ! while loop about convergence on tot_diff_minus+tot_diff_plus

  move = move + 1 

CONTAINS

SUBROUTINE EFFECTIVE_DIFFUSION (move)
  implicit none
  integer(i2b), intent(in) :: move
  select case ( move )
  case ( 1:100 )
    eD_plus = 0.1_dp*D_plus
    eD_minus = 0.1_dp*D_minus
    max_iter = 10
  case ( 101:500 )
    eD_plus = D_plus
    eD_minus = D_minus
    max_iter = 1
  case default
    if( D_plus < 0.03_dp ) then
      eD_plus = 0.03_dp
    else
      eD_plus = D_plus
    end if
    if( D_minus < 0.03_dp ) then
      eD_minus = 0.03_dp
    else
      eD_minus = D_minus
    end if
    max_iter = 1
  end select
END SUBROUTINE EFFECTIVE_DIFFUSION

END SUBROUTINE JUST_EQ_SMOLU