! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE POISSON_NERNST_PLANCK

  use precision_kinds, only: i2b, dp
  use system, only: D_equil, sigma, time, c_minus, c_minus_old, c_plus, c_plus_old, phi, phi_old, lx, ly, lz
  use input, only: input_dp, input_int
  use charges, only: charge_conservation
  implicit none

  ! read the number of iterations one does for the first step of equilibration (D_iter)
  D_equil = input_int('D_equil')
  if( D_equil <= 0 ) stop 'D_equil should not be <= 0.'

  ! read electrostatic related stuff
  call read_lncb_slope
  call read_elec_slope

  call allocate_old ! allocate phi_old, c_plus_old, c_minus_old, which convergence with time are checked

  ! iterate until charges are equilibrated
  time = 0
  equiloop: do while(.not.solutes_distrib_converged(time))
    print*,'time =',time
    ! backup potential and solute concentrations from last step
    call backup_phi_c_plus_c_minus
    ! compute phi with the Successive Over-relaxation Routine (SOR)
    call sor
    ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
    call just_eq_smolu ! called just_equ_smolu in c code
    if( sigma == 0 ) exit equiloop
    ! check charge conservation every 1000 steps (arbitrary number)
    if(modulo(time,1000)==0 .and. .not. charge_conservation() ) stop 'solute charge not conserved l.36'
    ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
    ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
    time = time + 1
  end do equiloop

  if( .not. charge_conservation() ) stop 'solute charge not conserved l.35'   ! check charge conservation
  if( input_int('printphi') == 1 ) call print4darray(lx, ly, lz, phi, 'output/phi_of_x_y_z.dat') ! print potential to file
  if( .not. solutes_distrib_converged(time) ) stop 'at the end of poisson_nernst_planck.f90, charges are still not converged'
  if( allocated(c_minus_old) ) deallocate(c_minus_old)
  if( allocated(c_plus_old)  ) deallocate(c_plus_old)
  if( allocated(phi_old) ) deallocate(phi_old)

STOP 'MANU MILIT'
CONTAINS






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! this subroutine checks if the charge distribution has reached equilibrium.
! Equivalent to check_equil in para_lattice.c
LOGICAL FUNCTION SOLUTES_DISTRIB_CONVERGED (time)
  implicit none
  integer(i2b), intent(in) :: time
  real(dp) :: dphi, dcp, dcm ! diff_phi, diff_c_plus and diff_c_minus
  integer(i2b) :: n, p, q
  real(dp), parameter :: eps=1.0d-12
  integer(i2b), save :: move = 0 ! first time
  move = move+1
  ! in this file is printed the convergence of phi due to the Sucessive OverRelaxation method
  select case(move)
  case(1)
    open(unit=99, file='output/Diff_phi_c', status='unknown')
    write(99,*)'# time, diff_phi, diff_c_plus, diff_c_minus'
    solutes_distrib_converged = .false.
    return ! dont check converge as there has been no step for now.
  case(2)
    solutes_distrib_converged = .false.
    return
  case default
    open(unit=99, file='output/Diff_phi_c', position='append', status='unknown')
  end select

  ! init differences
  dphi= 0.0_dp
  dcp = 0.0_dp
  dcm = 0.0_dp

  ! count the number of times the array is not zero
  n = count(abs(phi_old)>eps)
  p = count(abs(c_plus_old)>eps)
  q = count(abs(c_minus_old)>eps)

  ! at each node, if phi_old is different from 0, calculate the normalized relative differance between new and old potential
  if(n/=0) dphi = sum( abs(  (    phi -     phi_old )/    phi_old ), mask= abs(phi_old)>eps) / real(n,dp) ! this is in fact the L_1 norm of (phi-phi_old)/phi_old
  if(p/=0) dcp  = sum( abs(  ( c_plus -  c_plus_old )/ c_plus_old ), mask= abs(c_plus_old)>eps) / real(p,dp) ! L_1 norm is also called Manhattan norm and Taxicab norm
  if(q/=0) dcm  = sum( abs(  (c_minus - c_minus_old )/c_minus_old ), mask= abs(c_minus_old)>eps) / real(q,dp) ! it may be worth implementing a function called norm1 somewhere when I have time

  ! check convergence between previous and current step
  if( (dphi<eps .or. n==0) .and. (dcp<eps .or. p==0) .and. (dcm<eps .or. q==0) ) then ! TODO here is a very very idiot test
    solutes_distrib_converged = .true.
  else
    solutes_distrib_converged = .false.
  end if

  ! inform user
  print*,'time, dphi, dcp, dcm =',time,dphi,dcp,dcm
  write(99,*) time, dphi, dcp, dcm
  close(99)
END FUNCTION SOLUTES_DISTRIB_CONVERGED

SUBROUTINE BACKUP_PHI_C_PLUS_C_MINUS
  use system, only: phi, phi_old, c_plus, c_plus_old, c_minus, c_minus_old, lx, ly, lz
  implicit none
  ! backup
  c_plus_old = c_plus
  c_minus_old = c_minus
  phi_old = phi
END SUBROUTINE BACKUP_PHI_C_PLUS_C_MINUS

SUBROUTINE ALLOCATE_OLD
  implicit none
  ! the original distribution of concentrations should already be initiated before backing it up
  ! if it is the first iteration in find_equilibrium_charge_distribution then allocate backup arrays
  if( .not. allocated( c_plus_old) ) allocate( c_plus_old (lx,ly,lz) )
  if( .not. allocated( c_minus_old) ) allocate( c_minus_old (lx,ly,lz) )
  if( .not. allocated( phi_old) ) allocate( phi_old (lx,ly,lz) )
END SUBROUTINE ALLOCATE_OLD


END SUBROUTINE POISSON_NERNST_PLANCK