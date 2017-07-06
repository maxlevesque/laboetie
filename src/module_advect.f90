! Advect charge density using the volumetric interpretation of the lattice boltzmann
module module_advect
    implicit none
    private
    public :: advect
contains

subroutine advect( solventDensity)
  
  use precision_kinds, only: dp
  use constants, only: x, y, z
  use system, only: c_plus, c_minus, flux_site_minus, flux_site_plus, el_curr_x, el_curr_y, el_curr_z, &
                    ion_curr_x, ion_curr_y, ion_curr_z, fluid, pbc, supercell, node
  use myallocations, only: allocateReal3D
  
  implicit none
  
  real(dp), intent(in) :: solventDensity(:,:,:)
  integer :: i, j, k, ix, iy, iz, ip, jp, kp, ip_all(-1:1), jp_all(-1:1), kp_all(-1:1)
  real(dp) :: vx, vy, vz, ax, ay, az
  real(dp) :: flux_link_minus, flux_link_plus
  integer :: natureNodeijk, natureNodeipjpkp ! dummy variables
  real(dp) :: c_plus_ijk, c_minus_ijk
  real(dp), parameter :: zerodp = 0._dp, epsdp = epsilon(1._dp)

  if (.not. allocated(flux_site_plus)) call allocateReal3D(flux_site_plus)
  if (.not. allocated(flux_site_minus)) call allocateReal3D(flux_site_minus)
  
  flux_site_plus = zerodp
  flux_site_minus = zerodp

  if (.not. allocated(c_plus))  call allocateReal3D( c_plus)
  if (.not. allocated(c_minus)) call allocateReal3D( c_minus)


  !  THIS IS VERY CRIPTIC...
  !  What I am doing is what discussed with you.
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
  if ( all(abs(c_plus)<=epsdp) .and. all(abs(c_minus)<=epsdp) ) then ! advection of charge is useless since there is no charge!
      return
  end if

  if (all(abs(c_plus-c_minus)<=epsdp)) then
      el_curr_x = 0._dp
      el_curr_y = 0._dp
      el_curr_z = 0._dp
  else
      el_curr_x  = -sum( node%solventFlux(1)/solventDensity*( c_plus - c_minus ) , mask=node%nature==fluid)
      el_curr_y  = -sum( node%solventFlux(2)/solventDensity*( c_plus - c_minus ) , mask=node%nature==fluid)
      el_curr_z  = -sum( node%solventFlux(3)/solventDensity*( c_plus - c_minus ) , mask=node%nature==fluid)
  end if

  if (all(abs(c_plus+c_minus)<=epsdp)) then
      ion_curr_x = 0._dp
      ion_curr_y = 0._dp
      ion_curr_z = 0._dp
  else
      ion_curr_x  = -sum( node%solventFlux(1)/solventDensity*( c_plus + c_minus ) , mask=node%nature==fluid)
      ion_curr_y  = -sum( node%solventFlux(2)/solventDensity*( c_plus + c_minus ) , mask=node%nature==fluid)
      ion_curr_z  = -sum( node%solventFlux(3)/solventDensity*( c_plus + c_minus ) , mask=node%nature==fluid)
  end if


  ! ! compute concentration before advection step
  ! c_plus_total_old  = sum( c_plus)
  ! c_minus_total_old = sum( c_minus)

  do i= supercell%geometry%dimensions%indiceMin(x), supercell%geometry%dimensions%indiceMax(x)
    ip_all(:) = [(pbc(i+ix,x),ix=-1,1)]
    
    do j= supercell%geometry%dimensions%indiceMin(y), supercell%geometry%dimensions%indiceMax(y)
      jp_all(:) = [(pbc(j+iy,y),iy=-1,1)]
      
      do k= supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
        kp_all(:) = [(pbc(k+iz,z),iz=-1,1)]

        if( abs(c_plus(i,j,k))<=epsilon(1._dp) .and. abs(c_minus(i,j,k))<=epsilon(1._dp) ) cycle
        
        natureNodeijk = node(i,j,k)%nature
        c_plus_ijk = c_plus(i,j,k)
        c_minus_ijk = c_minus(i,j,k)

        ! velocities
        if( natureNodeijk == fluid) then ! Avoid division by zero where solventDensity is 0, i.e. in solid nodes
            vx = node(i,j,k)%solventFlux(1) / solventDensity(i,j,k)
            vy = node(i,j,k)%solventFlux(2) / solventDensity(i,j,k)
            vz = node(i,j,k)%solventFlux(3) / solventDensity(i,j,k)
        else
            vx = node(i,j,k)%solventFlux(1)
            vy = node(i,j,k)%solventFlux(2)
            vz = node(i,j,k)%solventFlux(3)
        endif
 
        do ix= -1, 1 ! TODO improve here by getting things out of useless loops
            ax = vx * real(ix,dp)
            ip = ip_all(ix) 

            do iy= -1, 1
                ay = vy * real(iy,dp)
                jp = jp_all(iy)
            
                do iz= -1, 1
                    az = vz * real(iz,dp)
                    kp = kp_all(iz)

                    natureNodeipjpkp = node(ip,jp,kp)%nature
                    
                    ! check if link is accessible
                    if ( ax>=0._dp .and. ay>=0._dp .and. az>=0._dp .and. natureNodeijk==fluid .and. natureNodeipjpkp==fluid) then ! link is accessible
                        if( ix == 0 ) ax = 1.0_dp - abs(vx)
                        if( iy == 0 ) ay = 1.0_dp - abs(vy)
                        if( iz == 0 ) az = 1.0_dp - abs(vz)
                        flux_link_plus  = ax*ay*az* c_plus_ijk
                        flux_link_minus = ax*ay*az* c_minus_ijk
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
  c_minus = c_minus + flux_site_minus
  c_plus  = c_plus  + flux_site_plus

  ! compute concentrations after advection
  ! c_minus_total_new = sum(c_minus)
  ! c_plus_total_new = sum(c_plus)


  ! check that total concentration has not changed
  !if( abs( c_minus_total_new - c_minus_total_old ) > 1.e-8 .or. abs( c_plus_total_new - c_plus_total_old ) > 1.e-8 ) then ! TODO magic number
  !  stop 'Total concentration has changed in advect.f90. STOP but unsure STOP is needed'
  !end if

end subroutine advect
end module module_advect
