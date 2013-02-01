! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

subroutine poisson_nernst_planck

  use precision_kinds, only: i2b, dp
  use system, only: D_equil, sigma, time, c_minus, c_minus_old, c_plus, c_plus_old, phi, phi_old, lx, ly, lz, inside, solid, fluid
  use input, only: input_dp, input_int
  use charges, only: charge_conservation
  implicit none

  ! read the number of iterations one does for the first step of equilibration (D_iter)
  D_equil = input_int('D_equil')
  if( D_equil <= 0 ) stop 'D_equil should not be <= 0.'

  ! read electrostatic related stuff
  call read_lncb_slope
  call read_elec_slope

  ! iterate until charges are equilibrated
  equiloop: do time = -D_equil, 0 ! time 0 being first step with flux

    print*,'time =',time

    ! backup potential and solute concentrations from last step
    call backup_phi_c_plus_c_minus

    ! compute phi with the Successive Over-relaxation Routine (SOR)
    call sor

    ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
    call just_eq_smolu ! called just_equ_smolu in c code

    ! check charge conservation every 1000 steps (arbitrary number)
    if(modulo(time,1000)==0 .and. .not. charge_conservation() ) stop 'solute charge not conserved l.36'

    if( sigma == 0 ) exit equiloop
    ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
    ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
    if((time/=-D_equil .and. modulo(time,1)==0) .and. solutes_distrib_converged(time) ) then
      print*,'solute distribution converged at step ',time,' after',abs(-D_equil-time),' steps'
      exit equiloop
    end if

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
  integer(kind=i2b), intent(in) :: time
  real(kind=dp) :: dphi, dcp, dcm ! diff_phi, diff_c_plus and diff_c_minus
  integer(kind=i2b) :: n, p, q
  real(kind=dp), parameter :: eps=1.0d-12
  logical, save :: is_firsttimehere = .true. ! first time 

  ! in this file is printed the convergence of phi due to the Sucessive OverRelaxation method
  if(is_firsttimehere) then
    open(unit=99, file='output/Diff_phi_c', status='unknown')
    write(99,*)'# time, diff_phi, diff_c_plus, diff_c_minus'
  else
    open(unit=99, file='output/Diff_phi_c', position='append', status='unknown')
  end if
  is_firsttimehere = .false.  ! now turn logical off to say the subroutine has already been called once

  ! init differences
  dphi =0.0_dp
  dcp =0.0_dp
  dcm =0.0_dp

  ! count the number of times the array is not zero
  n = count(abs(phi_old)>eps)
  p = count(abs(c_plus_old)>eps)
  q = count(abs(c_minus_old)>eps)

  ! at each node, if phi_old is different from 0, calculate the normalized relative differance between new and old potential
  if(n/=0) dphi = sum( abs(  (    phi -     phi_old )/    phi_old ), mask= abs(phi_old)>eps) / real(n,kind=dp) ! this is in fact the L_1 norm of (phi-phi_old)/phi_old
  if(p/=0) dcp  = sum( abs(  ( c_plus -  c_plus_old )/ c_plus_old ), mask= abs(c_plus_old)>eps) / real(p,kind=dp) ! L_1 norm is also called Manhattan norm and Taxicab norm
  if(q/=0) dcm  = sum( abs(  (c_minus - c_minus_old )/c_minus_old ), mask= abs(c_minus_old)>eps) / real(q,kind=dp) ! it may be worth implementing a function called norm1 somewhere when I have time

  ! check convergence between previous and current step
  if( (dphi<eps .or. n==0) .and. (dcp<eps .or. p==0) .and. (dcm<eps .or. q==0) ) then ! TODO here is a very very idiot test
    solutes_distrib_converged = .true.
  else
    solutes_distrib_converged = .false.
  end if

  ! inform user
  print*,'    diff_phi     =',dphi, n
  print*,'    diff_c_plus  =',dcp, p
  print*,'    diff_c_minus =',dcm, q

  ! write to diff_c
  write(99,*) time, dphi, dcp, dcm
  close(99)

END FUNCTION SOLUTES_DISTRIB_CONVERGED

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE POISSON_NERNST_PLANCK
