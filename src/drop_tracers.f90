! Here tracers are droped in the fluid. They may be charged or not. They evoluate in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
! Imagine a very small droplet of radioactive particles, so few they do not change
! anything to the system, but numerous enough to be followed and make statistics.

subroutine drop_tracers

  use precision_kinds, only: dp, i2b
  use system, only: tmom, tmax, elec_slope
  use populations, only: calc_n_momprop
  use moment_propagation, only: init, propagate, deallocate_propagated_quantity!, print_vacf, integrate_vacf!, not_yet_converged
  use input, only: input_dp

  implicit none
  integer(kind=i2b) :: it
  logical :: is_converged

  print*,'       step           VACF(x)                   VACF(y)                   VACF(z)'
  print*,'       ----------------------------------------------------------------------------------'

  ! include elec_slope in population n
  call calc_n_momprop

  ! turn the electric field off for the rest of mom_prop (included in n)
  elec_slope = 0.0_dp

  ! add electrostatic potential computed by the SOR routine an external contribution
  ! elec_pot(singlx,ly,lz, ch, phi, phi2, t, t_equil);
  ! call elec_pot

  call init ! init moment propagation


  ! propagate in time
  momproploop: do it= 1, tmax-tmom
!  it=0
!  do while (not_yet_converged(it))
!   call elec_pot
    call propagate (it,is_converged) ! propagate the propagated quantity
!    if( modulo(it,50000)==0 ) print_vacf
!    it = it + 1
    if( is_converged ) exit momproploop
  end do momproploop

  print*,

!  call integrate_vacf ! compute the integral over time of vacf in each direction
!  call print_vacf ! print vacf to file vacf.dat
  call deallocate_propagated_quantity

end subroutine drop_tracers
