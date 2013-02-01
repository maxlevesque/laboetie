! Here we compute the total electrostatic potential which is the sum
! of the internal potential (phi) computed in routine sor by the
! successive over relaxation method, and of the external electrostatic
! potential.

subroutine electrostatic_pot

  use precision_kinds, only: i2b, dp
  use system, only: elec_slope_x, elec_slope_y, elec_slope_z, phi, phi_tot, lx, ly, lz

  implicit none
  integer(kind=i2b) :: i, j, k ! dummy

  ! if no external electrostatics, don't loop: phi_tot = phi calculated by SOR
  if( elec_slope_x==0.0_dp .and. elec_slope_y==0.0_dp .and. elec_slope_z==0.0_dp ) then
    phi_tot = phi
  else
    ! else compute phi_tot as the sum of phi calculated by SOR and external contribution
    do concurrent (i=1:lx, j=1:ly, k=1:lz)
      phi_tot(i,j,k) = phi(i,j,k) + elec_slope_x*real(i,kind=dp) + elec_slope_y*real(j,kind=dp) + elec_slope_z*real(k,kind=dp)
    end do
  end if

end subroutine electrostatic_pot
