subroutine smolu
  use precision_kinds, only: dp, i2b
  use system, only: time, D_iter, D_plus, D_minus, solute_force,&
                     elec_slope, lncb_slope, phi_tot, &
                     kbt, fluid, solid, c_plus, c_minus, el_curr_x, el_curr_y, el_curr_z,&
                     ion_curr_x, ion_curr_y, ion_curr_z, pbc, supercell, node
  use constants, only: x, y, z, zerodp
  use mod_lbmodel, only: lbm
  use myallocations
  implicit none
  integer(kind=i2b) :: i, j, k, ip, jp, kp, l, m
  real(kind=dp) :: exp_dphi, exp_min_dphi, exp_dlncb, exp_min_dlncb
  real(dp), dimension(:,:,:), allocatable :: flux_site_minus, flux_site_plus, phiTEMP
  integer(kind=i2b) :: n_fluidsites ! total number of fluid side
  real(kind=dp) :: el_curr, ion_curr
  real(kind=dp) :: f_microions, f_plus, f_minus
  real(kind=dp) :: flux_link_plus, flux_link_minus
  integer(i2b) :: lx, ly, lz, half
  lx = supercell%geometry%dimensions%indiceMax(x)
  ly = supercell%geometry%dimensions%indiceMax(y)
  lz = supercell%geometry%dimensions%indiceMax(z)
  call allocateReal3D( flux_site_plus)
  call allocateReal3D( flux_site_minus)
  call allocateReal3D( phiTEMP)

  if(.not.allocated(solute_force)) allocate(solute_force(lx,ly,lz,x:z),source=0.0_dp)
  ! --------------------------- Ade 16/02/2017 ---------------------------------- 
  if (.not. allocated(c_plus)) call allocateReal3D(c_plus)
  if (.not. allocated(c_minus)) call allocateReal3D(c_minus)
  ! --------------------------- Ade 16/02/2017 ----------------------------------
  flux_site_plus = zerodp
  flux_site_minus = zerodp
  solute_force(:,:,:,1) = zerodp 
  solute_force(:,:,:,2) = zerodp
  solute_force(:,:,:,3) = zerodp

  !OPEN( 1616, FILE="output/TRUEFALSE.dat" )
  !OPEN( 1717, FILE="output/Debug.dat" )
  !OPEN( 1718, FILE="output/slp.dat" )
  
  !half = (lz/2)
  !DO m=1,ly
  !   WRITE(1718,*) m, (phi_tot(:,m,half))
  !ENDDO

  ! compute the flux between site i,j,k and site ip,jp,kp
  do i = 1, lx
    do j = 1, ly
      do k = 1, lz

        ! find arrival sites corresponding to each velocity
        do l= lbm%lmin+1, lbm%lmax ! not first l which is (ip,jp,kp)==(i,j,k). One may begin at l=1 but cycle if delta(l)==0. For instance in a where(delta/=0.0_dp)
          ! apply periodic boundary conditions
          ip= pbc( i+ lbm%vel(l)%coo(x) ,x)
          jp= pbc( j+ lbm%vel(l)%coo(y) ,y)
          kp= pbc( k+ lbm%vel(l)%coo(z) ,z)

           !write(1616,*) 'node(ip,jp,k)%nature == solid ', node(ip,jp,kp)%nature == solid
           !write(1616,*) 'node(i,j,k)%nature ==fluid ', node(i,j,k)%nature ==fluid
           !write(1616,*) 'Both ', node(i,j,k)%nature ==fluid .and. node(ip,jp,kp)%nature == solid

          ! if both nodes are in fluid
          if( node(i,j,k)%nature == fluid .and. node(ip,jp,kp)%nature == fluid ) then
            ! compute potential difference between sites
            exp_dphi = exp( phi_tot(ip,jp,kp) - phi_tot(i,j,k) ) ! arrival minus departure
            !if (l==lbm%lmin+3) phiTEMP(i,j,k) =  phi_tot(ip,jp,kp) - phi_tot(i,j,k)
            if (l==lbm%lmin+3) phiTEMP(i,j,k) = phi_tot(i,j,k)
            ! here is a very bizarre trick to correct for the jump in the external potential (elec_slope)
            if( i == lx .and. ip == 1  ) exp_dphi = exp_dphi* exp( elec_slope(x)* (lx+1) ) ! Ade : we added +1 here and below
            if( j == ly .and. jp == 1  ) exp_dphi = exp_dphi* exp( elec_slope(y)* (ly+1) )
            if( k == lz .and. kp == 1  ) exp_dphi = exp_dphi* exp( elec_slope(z)* (lz+1) )
            if( i == 1  .and. ip == lx ) exp_dphi = exp_dphi* exp(-elec_slope(x)* (lx+1) )
            if( j == 1  .and. jp == ly ) exp_dphi = exp_dphi* exp(-elec_slope(y)* (ly+1) )
            if( k == 1  .and. kp == lz ) exp_dphi = exp_dphi* exp(-elec_slope(z)* (lz+1) )
            ! phiTEMP = log(exp_dphi)
            ! inverse
            exp_min_dphi = 1.0_dp / exp_dphi

            ! contribution of the applied salt gradient (contribution of Magali)
            exp_dlncb = 1.0_dp !exp (sum (lncb_slope*lbm%vel(l)%coo(:)))! Ade : we modified this part and set it equal to one instead of the rest
                                                                        ! however, it seems to correpond to the C-code version given by Amael at 
                                                                        ! line 317 of smolu.c
            exp_min_dlncb = 1.0_dp / exp_dlncb  ! inverse

            ! flux due to electrostatic potential, density gradients and applied salt gradient
            flux_link_plus  = 0.5_dp*( 1.0_dp + exp_min_dphi * exp_min_dlncb )&
                                    *( c_plus(ip,jp,kp) *exp_dphi *exp_dlncb -c_plus(i,j,k) )
            flux_link_minus = 0.5_dp*( 1.0_dp + exp_dphi     * exp_min_dlncb )&
                                    *( c_minus(ip,jp,kp)*exp_min_dphi*exp_dlncb -c_minus(i,j,k) )
            el_curr  = flux_link_plus * D_plus - flux_link_minus * D_minus
            ion_curr = flux_link_plus * D_plus + flux_link_minus * D_minus
            ! forces exerted by solute on fluid
            f_plus  = 0.5_dp * ( 1.0_dp + exp_min_dphi )*( c_plus(ip,jp,kp) * exp_dphi - c_plus(i,j,k) )
            f_minus = 0.5_dp * ( 1.0_dp + exp_dphi )*( c_minus(ip,jp,kp) * exp_min_dphi - c_minus(i,j,k) )
            f_microions = kBT * (f_plus + f_minus)
            ! Subtract a term equal to the gradient of the charged densities
            f_microions = f_microions -kBT*(c_plus(ip,jp,kp)-c_plus(i,j,k)+c_minus(ip,jp,kp)-c_minus(i,j,k))

            flux_link_plus = flux_link_plus * (D_plus / lbm%vel(l)%delta )
            flux_link_minus = flux_link_minus * (D_minus / lbm%vel(l)%delta )

            flux_site_plus(i,j,k) = flux_site_plus(i,j,k) + flux_link_plus
            flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus
            ! -------------------- Ade : 8-02-17 -------------------------------------
            flux_site_plus(ip,jp,kp) = flux_site_plus(ip,jp,kp) - flux_link_plus
            flux_site_minus(ip,jp,kp) = flux_site_minus(ip,jp,kp) - flux_link_minus
            ! -------------------- Ade : 8-02-17 -------------------------------------
            ! flux. real flux is arriving at node. Real flux and el_curr are opposite.
            el_curr_x  = el_curr_x + lbm%vel(l)%a1 *lbm%vel(l)%coo(x) *el_curr / D_iter
            el_curr_y  = el_curr_y + lbm%vel(l)%a1 *lbm%vel(l)%coo(y) *el_curr / D_iter
            el_curr_z  = el_curr_z + lbm%vel(l)%a1 *lbm%vel(l)%coo(z) *el_curr / D_iter
            ion_curr_x = ion_curr_x + lbm%vel(l)%a1 *lbm%vel(l)%coo(x) *ion_curr / D_iter
            ion_curr_y = ion_curr_y + lbm%vel(l)%a1 *lbm%vel(l)%coo(y) *ion_curr / D_iter
            ion_curr_z = ion_curr_z + lbm%vel(l)%a1 *lbm%vel(l)%coo(z) *ion_curr / D_iter
            !rite(1717,*) 'exp_dphi1 =', flux_link_plus, flux_site_plus(:,:,k), f_plus, exp_min_dphi, c_plus(:,:,k), exp_dphi

            ! force exerted on fluid
            !print*, '!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!'
            !DO m=1,lz
            !      print*, 'FIRST solute_force = ', solute_force(:,:,k,3) ! Ade : The fluid is moving in the y-direction whenever a slit 
            !END DO
            !print*, '!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!'
            solute_force(i,j,k,:) = solute_force(i,j,k,:) + lbm%vel(l)%a1 *lbm%vel(l)%coo(:) *f_microions/D_iter

          else if( node(i,j,k)%nature ==fluid .and. node(ip,jp,kp)%nature == solid) then
             exp_dphi = exp(phi_tot(ip,jp,kp)-phi_tot(i,j,k))
            ! here is a very bizarre trick to correct for the jump in the external potential (elec_slope)
            if( i == lx .and. ip == 1  ) exp_dphi = exp_dphi* exp( elec_slope(x)* (lx+1) ) ! Ade : we added +1 here and below
            if( j == ly .and. jp == 1  ) exp_dphi = exp_dphi* exp( elec_slope(y)* (ly+1) )
            if( k == lz .and. kp == 1  ) exp_dphi = exp_dphi* exp( elec_slope(z)* (lz+1) )
            if( i == 1  .and. ip == lx ) exp_dphi = exp_dphi* exp(-elec_slope(x)* (lx+1) )
            if( j == 1  .and. jp == ly ) exp_dphi = exp_dphi* exp(-elec_slope(y)* (ly+1) )
            if( k == 1  .and. kp == lz ) exp_dphi = exp_dphi* exp(-elec_slope(z)* (lz+1) )

!            dphi = phi(ip,jp,kp)-phi(i,j,k)

!            ! If needed correct for the jump in the external potential (slope)
!            if(i==lx .and. ip==1) dphi = dphi+elec_slope(x)*lx
!            if(i==1 .and. ip==lx) dphi = dphi-elec_slope(x)*lx
!            if(j==ly .and. jp==1) dphi = dphi+elec_slope(y)*ly
!            if(j==1 .and. jp==ly) dphi = dphi-elec_slope(y)*ly
!            if(k==lz .and. kp==1) dphi = dphi+elec_slope(z)*lz
!            if(k==1 .and. kp==lz) dphi = dphi-elec_slope(z)*lz
            !write(1717,*) 'exp_dphi2 =', exp_dphi
            !write(1717,*) 'f_plus(1,1,3) =', f_plus(1,1,3)
            ! force acting on the solid. NO contribution of lncb_slope
            exp_min_dphi = 1.0_dp/exp_dphi
            f_plus = c_plus(i,j,k)*(1.0_dp-exp_min_dphi)
            f_minus = c_minus(i,j,k)*(1.0_dp-exp_dphi)
            f_microions = kBT*(f_plus + f_minus)
            solute_force(i,j,k,:) = solute_force(i,j,k,:) + lbm%vel(l)%a1*lbm%vel(l)%coo(:)*f_microions/D_iter
            !print*, '!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!'
            !DO m=1,lz
            !      print*, 'SECOND solute_force = ', solute_force(:,:,k,3) ! Ade : The fluid is moving in the y-direction whenever a slit 
            !END DO
            !print*, '!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!'

          end if
 
        end do ! loop over neighbours
       
      end do
    end do
  end do
  !half = (lz/2)
  !DO m=1,ly
  !   WRITE(1718,*) m, (phiTEMP(:,m,half))
  !ENDDO
  !close(1616)
  !lose(1717)
  !close(1718)

  ! update concentrations. Smoluchowski part.
  where(node%nature==fluid)
    c_plus = c_plus + flux_site_plus
    c_minus = c_minus + flux_site_minus
  end where

  ! normalize by volume of fluid phase
  if( time > 0 ) then
    n_fluidsites = sum(node%nature, mask=node%nature==fluid)
    el_curr_x = el_curr_x / n_fluidsites
    el_curr_y = el_curr_y / n_fluidsites
    el_curr_z = el_curr_z / n_fluidsites
    ion_curr_x = ion_curr_x / n_fluidsites
    ion_curr_y = ion_curr_y / n_fluidsites
    ion_curr_z = ion_curr_z / n_fluidsites
  end if

  ! check that the sum of all fluxes is zero
  if( abs(sum( flux_site_plus,mask=node%nature==fluid)) > 1.0e-12 .or. &
      abs(sum( flux_site_minus,mask=node%nature==fluid)) > 1.0e-12 ) then
      stop 'in smolu.f90 the sum of all fluxes does not add up to zero! stop.'
  end if
end subroutine smolu
