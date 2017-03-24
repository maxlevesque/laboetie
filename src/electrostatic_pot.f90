! Here we compute the total electrostatic potential which is the sum
! of the internal potential (phi) computed in routine sor by the
! successive over relaxation method, and of the external electrostatic
! potential.

subroutine electrostatic_pot
    use precision_kinds, only: i2b, dp
    use system, only: elec_slope, phi, phi_tot, supercell
    use constants, only: x, y, z, zerodp
    USE myallocations
    implicit none
    integer(kind=i2b) :: i, j, k, lz, m, half, ly, lx ! dummy
    ! ----------- Ade : 21/02/17 ------------------------------
    if(.not.allocated(phi_tot)) call allocateReal3D(phi_tot)
    phi_tot(:,:,:) = zerodp 
    !open(375, file = "output/phiELECPOT.dat")
    lx = supercell%geometry%dimensions%indiceMax(x)
    ly = supercell%geometry%dimensions%indiceMax(y)
    lz = supercell%geometry%dimensions%indiceMax(z)
    half = lz/2
    ! ----------- Ade : 21/02/17 ------------------------------

    ! if no external electrostatics, don't loop: phi_tot = phi calculated by SOR

    if( all(elec_slope==0.0_dp) ) then
        phi_tot = phi
    else
    !DO m=1,lz
       !WRITE(375,*) 'Elec_slope = ', elec_slope(1), elec_slope(2), elec_slope(3)
       !WRITE(375,*) m, SUM(phi(:,:,m))
    !ENDDO
        ! else compute phi_tot as the sum of phi calculated by SOR and external contribution
        do concurrent( k= 1:supercell%geometry%dimensions%indiceMax(z),&
                       j= 1:supercell%geometry%dimensions%indiceMax(y),&
                       i= 1:supercell%geometry%dimensions%indiceMax(x) )
            phi_tot(i,j,k) = phi(i,j,k) + ( elec_slope(1) * real(i-1,dp) + elec_slope(2)*real(j-1,dp) + elec_slope(3)*real(k-1,dp) )
        end do
        !write(375,*) '# half = ', half
        !DO m=1,ly
        !    WRITE(375,*) m, (phi_tot(:,m,half))
        !ENDDO
    end if
    !close(375)
end subroutine electrostatic_pot
