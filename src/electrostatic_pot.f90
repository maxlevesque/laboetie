! Here we compute the total electrostatic potential.
! It is the sum of the internal potential, phi, that should have been computed before for instance by Successive Over Relaxation,
! and of an some external electric potential.

subroutine electrostatic_pot

    use precision_kinds, only: dp
    use system, only: elec_slope, phi, phi_tot, supercell

    implicit none

    integer :: i, j, k, lx, ly, lz
    real(dp) :: potx, poty, potz
    real(dp), parameter :: zerodp = 0._dp

    lx = supercell%geometry%dimensions%indiceMax(1)
    ly = supercell%geometry%dimensions%indiceMax(2)
    lz = supercell%geometry%dimensions%indiceMax(3)

    if( .not. allocated(phi_tot)) allocate( phi_tot(lx,ly,lz), source=zerodp)

    ! If no external electrostatics, don't loop: phi_tot = phi calculated by SOR
    ! else compute phi_tot as the sum of phi calculated by SOR and external contribution

    if( all( elec_slope == 0.0_dp ) ) then
        phi_tot = phi
    else
        do k = 1, lz
            potz = real(k-1,dp) * elec_slope(3)
            do j = 1, ly
                poty = real(j-1,dp) * elec_slope(2)
                do i = 1, lx
                    potx = real(i-1,dp) * elec_slope(1)
                    phi_tot(i,j,k) = phi(i,j,k) + potx + poty + potz
                end do
            end do
        end do
    end if

end subroutine electrostatic_pot
