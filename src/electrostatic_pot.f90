! Here we compute the total electrostatic potential which is the sum
! of the internal potential (phi) computed in routine sor by the
! successive over relaxation method, and of the external electrostatic
! potential.

subroutine electrostatic_pot

  use precision_kinds, only: i2b, dp
  use system, only: elec_slope, phi, phi_tot, lx, ly, lz

  implicit none
  integer(kind=i2b) :: i, j, k ! dummy

  ! if no external electrostatics, don't loop: phi_tot = phi calculated by SOR
  if( all(elec_slope==0.0_dp) ) then
    phi_tot = phi
  else
    ! else compute phi_tot as the sum of phi calculated by SOR and external contribution
    do concurrent (i=1:lx, j=1:ly, k=1:lz)
      phi_tot(i,j,k) = phi(i,j,k) + sum( elec_slope * real([i,j,k],dp) )
    end do
  end if

end subroutine electrostatic_pot
