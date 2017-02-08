! Advect charge density using the volumetric interpretation of the lattice boltzmann
subroutine advect
  use precision_kinds, only: dp, i2b
  use constants, only: x, y, z
  use system, only: c_plus, c_minus, flux_site_minus, flux_site_plus, el_curr_x, el_curr_y, el_curr_z, &
  ion_curr_x, ion_curr_y, ion_curr_z, fluid, pbc, supercell, node
  use myallocations
  implicit none
  real(dp) :: c_minus_total_old, c_plus_total_old, c_minus_total_new, c_plus_total_new
  integer(i2b) :: i, j, k, ix, iy, iz, ip, jp, kp, ip_all(-1:1), jp_all(-1:1), kp_all(-1:1)
  real(dp) :: vx, vy, vz, ax, ay, az
  real(dp) :: flux_link_minus, flux_link_plus
  call allocateReal3D( flux_site_minus )
  call allocateReal3D( flux_site_plus )

  ! Ade : This is a temporary test. To be removed
  !-------------------- Ade --------------------------- 
  ! init ion (solute) concentrations
  if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
  if (.not. allocated(c_minus)) call allocateReal3D( c_minus)
  !-------------------- Ade --------------------------- 

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
  if ( all(abs(c_plus)<=epsilon(1._dp)) .and. all(abs(c_minus)<=epsilon(1._dp)) ) then ! advection of charge is useless since there is no charge!
    return
  end if

  if (all(abs(c_plus-c_minus)<=epsilon(1._dp))) then
    el_curr_x = 0._dp
    el_curr_y = 0._dp
    el_curr_z = 0._dp
  else
    el_curr_x  = -sum( node%solventFlux(x)*( c_plus - c_minus ) , mask=node%nature==fluid)
    el_curr_y  = -sum( node%solventFlux(y)*( c_plus - c_minus ) , mask=node%nature==fluid)
    el_curr_z  = -sum( node%solventFlux(z)*( c_plus - c_minus ) , mask=node%nature==fluid)
  end if

  if (all(abs(c_plus+c_minus)<=epsilon(1._dp))) then
    el_curr_x = 0._dp
    el_curr_y = 0._dp
    el_curr_z = 0._dp
  else
    ion_curr_x  = -sum( node%solventFlux(x)*( c_plus + c_minus ) , mask=node%nature==fluid)
    ion_curr_y  = -sum( node%solventFlux(y)*( c_plus + c_minus ) , mask=node%nature==fluid)
    ion_curr_z  = -sum( node%solventFlux(z)*( c_plus + c_minus ) , mask=node%nature==fluid)
  end if

  ! compute concentration before advection step
  c_plus_total_old = sum( c_plus)
  c_minus_total_old = sum( c_minus)



  !-------------------- Ade --------------------------- 
  ip_all = [0,0,0]
  jp_all = [0,0,0]
  kp_all = [0,0,0]

  ix = 0
  iy = 0
  iz = 0

  ip = 0
  jp = 0
  kp = 0
  !-------------------- Ade --------------------------- 
  
  do i= supercell%geometry%dimensions%indiceMin(x), supercell%geometry%dimensions%indiceMax(x)
    
    !do ix = -1,1
    !    ip_all(ix) = pbc( i+ix ,x)
    !    print*, 'ix dans la boucle vaut ', ix
    !end do     
    ip_all(:) = [(pbc(i+ix,x),ix=-1,1)]
    
    do j= supercell%geometry%dimensions%indiceMin(y), supercell%geometry%dimensions%indiceMax(y)
      
      !do iy = -1,1
      !    jp_all(iy) = pbc( j+iy ,y)
      !    print*, 'iy dans la boucle vaut ', iy
      !end do     
      jp_all(:) = [(pbc(j+iy,y),iy=-1,1)]
      
      do k= supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
        if( abs(c_plus(i,j,k))<=epsilon(1._dp) .and. abs(c_minus(i,j,k))<=epsilon(1._dp) ) cycle
        
        !do iz = -1,1
        !    kp_all(iz) = pbc( k+iz ,z)
        !    print*, 'iz dans la boucle vaut ', iz
        !end do     
        kp_all(:) = [(pbc(k+iz,z),iz=-1,1)]
        
        ! velocities
        vx = node(i,j,k)%solventFlux(x)
        vy = node(i,j,k)%solventFlux(y)
        vz = node(i,j,k)%solventFlux(z)
        do ix= -1, 1 ! TODO improve here by getting things out of useless loops
          ax = vx * real(ix,dp)
          do iy= -1, 1
            ay = vy * real(iy,dp)
            do iz= -1, 1
              ip = ip_all(ix) ! Ade: there is a problem with ip!! The number is huge
              jp = jp_all(iy)
              kp = kp_all(iz)
              az = vz * real(iz,dp)
              ! check if link is accessible
              if (all([ax,ay,az]>=0._dp) .and. node(i,j,k)%nature==node(ip,jp,kp)%nature) then ! link is accessible
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
  !if( abs( c_minus_total_new - c_minus_total_old ) > 1.e-8 .or. abs( c_plus_total_new - c_plus_total_old ) > 1.e-8 ) then ! TODO magic number
  !  stop 'Total concentration has changed in advect.f90. STOP but unsure STOP is needed'
  !end if
end subroutine advect
