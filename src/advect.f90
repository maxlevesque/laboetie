! Advect charge density using the volumetric interpretation of the lattice boltzmann
subroutine advect
    use precision_kinds
    use constants, only: x, y, z
    use system, only: c_plus, c_minus, flux_site_minus, flux_site_plus, el_curr_x, el_curr_y, el_curr_z, &
                    ion_curr_x, ion_curr_y, ion_curr_z, fluid, pbc, supercell, node
    use myallocations
    implicit none
    real(dp) :: c_minus_total_old, c_plus_total_old, c_minus_total_new, c_plus_total_new
    integer(i2b) :: i, j, k, ix, iy, iz, ip, jp, kp
    real(dp) :: vx, vy, vz, ax, ay, az
    real(dp) :: flux_link_minus, flux_link_plus
    call allocateReal3D( flux_site_minus )
    call allocateReal3D( flux_site_plus )

!   THIS IS VERY CRIPTIC...
!   What I am doing is what discussed with you.
!  Basically (when using this word it means that the explanation is a little bit mistic :-) )
!  you have to compute the portion of the volume which goes in one of the neighbouring 26 cells.
!  The part of the volume that goes in the next cell is ax*ay*az.
!  ax could have three values.
!  ax = (1 - vx)
!  ax = vx
!  ax = 0 (zero)
!  If I multiply by 0 I do not continue.
!  One complication of the algorithm is that if ix (the x coordinate of the next cell) is 0,
!  ax is 0 but I should transport on y and z direction.
!  Forget this explanation. It works!!!!

    ! convective contribution to the electric current
    el_curr_x  = -sum( node%solventFlux(x)*( c_plus - c_minus ) , mask=node%nature==fluid)
    el_curr_y  = -sum( node%solventFlux(y)*( c_plus - c_minus ) , mask=node%nature==fluid)
    el_curr_z  = -sum( node%solventFlux(z)*( c_plus - c_minus ) , mask=node%nature==fluid)
    ion_curr_x = -sum( node%solventFlux(x)*( c_plus + c_minus ) , mask=node%nature==fluid)
    ion_curr_y = -sum( node%solventFlux(y)*( c_plus + c_minus ) , mask=node%nature==fluid)
    ion_curr_z = -sum( node%solventFlux(z)*( c_plus + c_minus ) , mask=node%nature==fluid)
    ! compute concentration before advection step
    c_plus_total_old = sum( c_plus)
    c_minus_total_old = sum( c_minus)
    do i= supercell%geometry%dimensions%indiceMin(x), supercell%geometry%dimensions%indiceMax(x)
        do j= supercell%geometry%dimensions%indiceMin(y), supercell%geometry%dimensions%indiceMax(y)
            do k= supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
                ! velocities
                vx = node(i, j, k)%solventFlux(x)
                vy = node(i, j, k)%solventFlux(y)
                vz = node(i, j, k)%solventFlux(z)
                do ix= -1, 1 ! TODO improve here by getting things out of useless loops
                    do iy= -1, 1
                        do iz= -1, 1
                            ip = pbc( i+ ix,x)
                            jp = pbc( j+ iy,y)
                            kp = pbc( k+ iz,z)
                            ax = vx * real(ix,dp)
                            ay = vy * real(iy,dp)
                            az = vz * real(iz,dp)
                            ! check if link is accessible
                            if( ax >= 0.0_dp .and. ay >= 0.0_dp .and. az >= 0.0_dp .and. &
                                 node(i,j,k)%nature == node(ip,jp,kp)%nature ) then ! link is accessible
                                if( ix == 0 ) ax = 1.0_dp - abs(vx)
                                if( iy == 0 ) ay = 1.0_dp - abs(vy)
                                if( iz == 0 ) az = 1.0_dp - abs(vz)
                                flux_link_plus  = ax*ay*az* c_plus (i,j,k)
                                flux_link_minus = ax*ay*az* c_minus(i,j,k)
                            else ! link not accessible => no flux
                                flux_link_plus  = 0.0_dp
                                flux_link_minus = 0.0_dp
                            end if
                            ! remove what is leaving site i,j,k and add what is arriving at ip,jp,kp
                            flux_site_minus(i,j,k)    = flux_site_minus(i,j,k)    - flux_link_minus
                            flux_site_minus(ip,jp,kp) = flux_site_minus(ip,jp,kp) + flux_link_minus
                            flux_site_plus (i,j,k)    = flux_site_plus (i,j,k)    - flux_link_plus
                            flux_site_plus (ip,jp,kp) = flux_site_plus (ip,jp,kp) + flux_link_plus
                        end do
                    end do
                end do
            end do
        end do
    end do

  ! update concentrations accordingly to flux calculated previously
  where( node%nature == fluid )
    c_minus = c_minus + flux_site_minus
    c_plus  = c_plus  + flux_site_plus
  end where

  ! compute concentrations after advection
  c_minus_total_new = sum(c_minus)
  c_plus_total_new = sum(c_plus)

    ! check that total concentration has not changed
    if( abs( c_minus_total_new - c_minus_total_old ) > 1.e-8 .or. &
        abs( c_plus_total_new - c_plus_total_old ) > 1.e-8 ) then ! TODO magic number
        print*, 'Total concentration has changed in advect.f90. STOP but unsure STOP is needed'
        stop ! not sure it should stop every time it changes.
    end if
end subroutine advect
