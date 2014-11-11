! now that the charges (if present) are equilibrated, one makes things move here.
! here we make things move but without any constraints, ie without external
! forces applied

subroutine equilibration_without_constraint

  use precision_kinds, only: i2b, dp
  use system, only: t_equil, f_ext, d_iter, time, elec_slope, lncb_slope
  use populations, only: update_populations
  use myallocations
  use constants, only: x, y, z

  implicit none
  integer(i2b) :: iteration

  print*,'       step       current(x)                current(y)                 current(z)'
  print*,'       ----------------------------------------------------------------------------------'

  f_ext = 0.0_dp
  elec_slope = 0.0_dp
  lncb_slope = 0.0_dp

  timeloop: do time = 1, t_equil

  call update_populations      ! get populations
  call propagation ! propagation of the distribution functions according to their direction to the next nodes
  call comp_rho    ! compute density
  call comp_j      ! compute momenta
  call advect      ! solute motion: advection step

  do iteration= 1, d_iter   ! solute motion: diffusion step
    call sor                ! compute phi with sucessive overrelaxation method
    call electrostatic_pot  ! compute phi + phi_external (due to electrostatic field elec_slope(x:z))
    call smolu              ! this time full smolu, not just_equ_smolu
    call charge_test        ! check charge conservation
  end do

  end do timeloop

end subroutine equilibration_without_constraint
