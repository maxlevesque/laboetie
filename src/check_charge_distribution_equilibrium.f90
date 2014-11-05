! this subroutine checks if the charge distribution has reached equilibrium.
! Equivalent to check_equil in para_lattice.c

subroutine check_charge_distribution_equilibrium( time, is_converged)

  use precision_kinds, only: dp, i2b
  use system, only: phi, phi_old, c_minus, c_minus_old, c_plus, c_plus_old, node
  implicit none
  integer(kind=i2b), intent(in) :: time
  logical, intent(out) :: is_converged
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
  if( dphi<eps .and. dcp<eps .and. dcm<eps ) then ! TODO here is a very very idiot test
    is_converged = .true.
  else
    is_converged = .false.
  end if

  ! inform user
  print*,'    diff_phi     =',dphi, n
  print*,'    diff_c_plus  =',dcp, p
  print*,'    diff_c_minus =',dcm, q

  ! write to diff_c
  write(99,*) time, dphi, dcp, dcm
  close(99)

  ! now turn logical off to say the subroutine has already been called
  is_firsttimehere = .false.

end subroutine check_charge_distribution_equilibrium
