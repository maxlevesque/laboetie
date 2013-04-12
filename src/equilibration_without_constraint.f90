! Now that the charges (if present) are equilibrated, one makes things move here.
! Here we make things move but without any constraints, ie without external
! forces applied

subroutine equilibration_without_constraint

  use precision_kinds, only: i2b, dp
  use system, only: t_equil, f_ext, D_iter, time, elec_slope, lncb_slope
  use populations, only: calc_n
  use myallocations
  use constants, only: x, y, z
  implicit none
  integer(kind=i2b) :: iteration
  print*,'       step       current(x)                current(y)                 current(z)'
  print*,'       ----------------------------------------------------------------------------------'

  ! In these equilibration steps, we do not apply the external forces
  f_ext = 0.0_dp ! f_ext(x:y)
  elec_slope = 0.0_dp
  lncb_slope = 0.0_dp
  ! do not forget to read them when equilibrium is found ie in equilibration_with_constraints

  ! start timer
  timeloop: do time = 1, t_equil

    call calc_n ! populations

    call propagation ! propagatoin of the distribution functions according to their direction to the next nodes

    ! compute density
    call comp_rho

    ! compute momenta
    call comp_j

    ! solute motion: advection step
    call advect

    ! solute motion: diffusion step
    do iteration= 1, D_iter ! D_iter is read before in input file
      ! compute phi with sucessive overrelaxation method
      call sor
      ! compute phi + phi_external (due to electrostatic field elec_slope(x:z))
      call electrostatic_pot
      ! this time full smolu, not just_equ_smolu
      call smolu
      ! check charge conservation
      call charge_test
    end do

  end do timeloop
end subroutine equilibration_without_constraint
