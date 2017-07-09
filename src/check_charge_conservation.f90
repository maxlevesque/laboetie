! In this routine one wants to be sure that the total charges the user
! asks as input (tot_sol_charge) is kept constant during the simulation

pure subroutine check_charge_conservation

  use precision_kinds, only: dp, sp
  use system, only: c_plus, c_minus, tot_sol_charge, fluid, solid, supercell, node
  implicit none
  real(dp) :: tolerance ! how much change in the charge is allowed because of the floating point precision.
  real(dp) :: solid_charge, fluid_charge ! should be equal to tot_sol_charge (asked by user)

  if( dp /= sp ) then
    tolerance = 1.E-10
  else ! if dp is defined as single precision
    tolerance = 1.E-6
  end if

  ! the total charge is c_plus - c_minus
  ! the total charge in solid is
  solid_charge = sum(c_plus-c_minus, mask=(node%nature==solid)) ! I love Fortran. (The sum is applied only where the node is "solid")
  ! the total charge in fluid is
  fluid_charge = sum(c_plus-c_minus, mask=(node%nature==fluid))

  ! check conservation
  if( abs(solid_charge-tot_sol_charge) > tolerance .or. abs(fluid_charge+tot_sol_charge) > tolerance ) then
    ! print*,'abs(solid_charge-tot_sol_charge) = ',abs(solid_charge-tot_sol_charge)
    ! print*,'abs(fluid_charge+tot_sol_charge) = ',abs(fluid_charge+tot_sol_charge)
    error stop 'charge is not conserved. check check_charge_conservation.f90. stop'
  end if

end subroutine check_charge_conservation
