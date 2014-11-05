! In this routine one wants to be sure that the total charges the user
! asks as input (sigma) is kept constant during the simulation

subroutine charge_test

  use precision_kinds, only: dp
  use system, only: c_plus, c_minus, sigma, fluid, solid, supercell, node
  implicit none
  real(kind=dp), parameter :: tolerance = 1.0e-10 ! how much change due to numeric in charge is allowed TODO magic number
  real(kind=dp) :: solid_charge, fluid_charge ! should be equal to sigma (asked by user)


  ! the total charge is c_plus - c_minus
  ! the total charge in solid is
  solid_charge = sum(c_plus-c_minus,mask=(node%nature==solid)) ! I love Fortran
  ! the total charge in fluid is
  fluid_charge = sum(c_plus-c_minus,mask=(node%nature==fluid))

  ! check conservation
  if( abs(solid_charge-sigma) > tolerance .or. abs(fluid_charge+sigma) > tolerance ) then
    print*,'abs(solid_charge-sigma) = ',abs(solid_charge-sigma)
    print*,'abs(fluid_charge+sigma) = ',abs(fluid_charge+sigma)
    stop 'charge is not conserved. check charge_test.f90. stop'
  end if

end subroutine charge_test
