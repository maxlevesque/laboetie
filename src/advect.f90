! Advect charge density using the volumetric interpretation of the lattice boltzmann

subroutine advect

  use precision_kinds
  use constants, only: x, y, z
  use system, only: c_plus, c_minus, flux_site_minus, flux_site_plus, el_curr_x, el_curr_y, el_curr_z, &
                    ion_curr_x, ion_curr_y, ion_curr_z, lx, ly, lz, inside, jx, jy, jz, fluid, pbc, supercell
  implicit none
  real(kind=dp) :: c_minus_total_old, c_plus_total_old, c_minus_total_new, c_plus_total_new
  integer(kind=i2b) :: i, j, k, ix, iy, iz, ip, jp, kp
  real(kind=dp) :: vx, vy, vz, ax, ay, az
  real(kind=dp) :: flux_link_minus, flux_link_plus

  supercell%node%nature = inside !TODO REMOVE

  ! set fluxes to zero. allocate them if necessary
  if (.not. allocated(flux_site_minus)) allocate( flux_site_minus(lx,ly,lz),source=0.0_dp )
  if (.not. allocated(flux_site_plus )) allocate( flux_site_plus (lx,ly,lz),source=0.0_dp )

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
  el_curr_x  = -sum( jx*( c_plus - c_minus ) , mask=supercell%node%nature==fluid)
  el_curr_y  = -sum( jy*( c_plus - c_minus ) , mask=supercell%node%nature==fluid)
  el_curr_z  = -sum( jz*( c_plus - c_minus ) , mask=supercell%node%nature==fluid)

  ion_curr_x = -sum( jx*( c_plus + c_minus ) , mask=supercell%node%nature==fluid)
  ion_curr_y = -sum( jy*( c_plus + c_minus ) , mask=supercell%node%nature==fluid)
  ion_curr_z = -sum( jz*( c_plus + c_minus ) , mask=supercell%node%nature==fluid)

  ! compute concentration before advection step
  c_plus_total_old = sum( c_plus)
  c_minus_total_old = sum( c_minus)

  do i= 1, lx
    do j= 1, ly
      do k= 1, lz

        ! velocities
        vx = jx( i, j, k)
        vy = jy( i, j, k)
        vz = jz( i, j, k)

        do ix= -1, 1 ! TODO improve here by getting things out of useless loops
          do iy= -1, 1
            do iz= -1, 1

              ip = pbc( i+ ix,x)
              jp = pbc( j+ iy,y)
              kp = pbc( k+ iz,z)

              ax = vx * real(ix,kind=dp)
              ay = vy * real(iy,kind=dp)
              az = vz * real(iz,kind=dp)

              ! check if link is accessible
              if( ax >= 0.0_dp .and. ay >= 0.0_dp .and. az >= 0.0_dp .and. &
                         supercell%node(i,j,k)%nature == supercell%node(ip,jp,kp)%nature ) then ! link is accessible

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
  where( supercell%node%nature == fluid )
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
