subroutine smolu

  use precision_kinds, only: dp, i2b
  use system, only: time, D_iter, inside, a1, c, D_plus, D_minus, delta, solute_force,&
                     elec_slope_x, elec_slope_y, elec_slope_z, lncb_slope_x, lncb_slope_y, lncb_slope_z,&
                     NbVel, phi_tot, lx, ly, lz,&
                     kbt, fluid, solid, c_plus, c_minus, el_curr_x, el_curr_y, el_curr_z,&
                     ion_curr_x, ion_curr_y, ion_curr_z, plusx, plusy, plusz
  use constants, only: x, y, z

  implicit none
  integer(kind=i2b) :: i, j, k, ip, jp, kp, l
  real(kind=dp) :: exp_dphi, exp_min_dphi, exp_dlncb, exp_min_dlncb
  real(kind=dp), dimension(lx, ly, lz) :: flux_site_minus, flux_site_plus
  integer(kind=i2b) :: n_fluidsites ! total number of fluid side
  real(kind=dp) :: el_curr, ion_curr
  real(kind=dp) :: f_microions, f_plus, f_minus
  real(kind=dp) :: flux_link_plus, flux_link_minus


  flux_site_plus = 0.0_dp
  flux_site_minus = 0.0_dp

  if(.not.allocated(solute_force)) allocate(solute_force(lx,ly,lz,x:z),source=0.0_dp)


  ! compute the flux between site i,j,k and site ip,jp,kp
  do i = 1, lx
    do j = 1, ly
      do k = 1, lz

        ! find arrival sites corresponding to each velocity
        do l = 2, NbVel ! not l==1 which is (ip,jp,kp)==(i,j,k). One may begin at l=1 but cycle if delta(l)==0. For instance in a where(delta/=0.0_dp)

          ! apply periodic boundary conditions
          ip= plusx( i+ c(x,l) )
          jp= plusy( j+ c(y,l) )
          kp= plusz( k+ c(z,l) )

          ! if both nodes are in fluid
          if( inside(i,j,k) == fluid .and. inside(ip,jp,kp) == fluid ) then

            ! compute potential difference between sites
            exp_dphi = exp( phi_tot(ip,jp,kp) - phi_tot(i,j,k) ) ! arrival minus departure

            ! here is a very bizarre trick to correct for the jump in the external potential (elec_slope)
            if( i == lx .and. ip == 1  ) exp_dphi = exp_dphi* exp( elec_slope_x* lx )
            if( j == ly .and. jp == 1  ) exp_dphi = exp_dphi* exp( elec_slope_y* ly )
            if( k == lz .and. kp == 1  ) exp_dphi = exp_dphi* exp( elec_slope_z* lz )
            if( i == 1  .and. ip == lx ) exp_dphi = exp_dphi* exp(-elec_slope_x* lx )
            if( j == 1  .and. jp == ly ) exp_dphi = exp_dphi* exp(-elec_slope_y* ly )
            if( k == 1  .and. kp == lz ) exp_dphi = exp_dphi* exp(-elec_slope_z* lz )

            ! inverse
            exp_min_dphi = 1.0_dp / exp_dphi

            ! contribution of the applied salt gradient (contribution of Magali)
            exp_dlncb = exp( lncb_slope_x*c(x,l) +lncb_slope_y*c(y,l) +lncb_slope_z*c(z,l) )
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

            flux_link_plus = flux_link_plus * (D_plus / delta(l))
            flux_link_minus = flux_link_minus * (D_minus / delta(l))

            flux_site_plus(i,j,k) = flux_site_plus(i,j,k) + flux_link_plus
            flux_site_minus(i,j,k) = flux_site_minus(i,j,k) + flux_link_minus

            ! flux. real flux is arriving at node. Real flux and el_curr are opposite. 
            el_curr_x  = el_curr_x +a1(l) *c(x,l) *el_curr / D_iter
            el_curr_y  = el_curr_y +a1(l) *c(y,l) *el_curr / D_iter
            el_curr_z  = el_curr_z +a1(l) *c(z,l) *el_curr / D_iter
            ion_curr_x = ion_curr_x +a1(l) *c(x,l) *ion_curr / D_iter
            ion_curr_y = ion_curr_y +a1(l) *c(y,l) *ion_curr / D_iter
            ion_curr_z = ion_curr_z +a1(l) *c(z,l) *ion_curr / D_iter

            ! force exerted on fluid
            solute_force(i,j,k,:) = solute_force(i,j,k,:) + a1(l)*c(:,l) *f_microions/D_iter

          else if( inside(i,j,k)==fluid .and. inside(ip,jp,kp)==solid) then
            
!            dphi = phi(ip,jp,kp)-phi(i,j,k)

!            ! If needed correct for the jump in the external potential (slope) 
!            if(i==lx .and. ip==1) dphi = dphi+elec_slope_x*lx
!            if(i==1 .and. ip==lx) dphi = dphi-elec_slope_x*lx
!            if(j==ly .and. jp==1) dphi = dphi+elec_slope_y*ly
!            if(j==1 .and. jp==ly) dphi = dphi-elec_slope_y*ly
!            if(k==lz .and. kp==1) dphi = dphi+elec_slope_z*lz
!            if(k==1 .and. kp==lz) dphi = dphi-elec_slope_z*lz

            ! force acting on the solid. NO contribution of lncb_slope
            exp_dphi = exp(phi_tot(ip,jp,kp)-phi_tot(i,j,k))
            exp_min_dphi = 1.0_dp/exp_dphi
            f_plus = c_plus(i,j,k)*(1.0_dp-exp_min_dphi)
            f_minus = c_minus(i,j,k)*(1.0_dp-exp_dphi)
            f_microions = kBT*(f_plus + f_minus)
            solute_force(i,j,k,:) = solute_force(i,j,k,:) +a1(l)*c(:,l)*f_microions/D_iter

          end if

        end do ! loop over neighbours

      end do
    end do
  end do

  ! update concentrations. Smoluchowski part.
  where(inside==fluid)
    c_plus = c_plus + flux_site_plus
    c_minus = c_minus + flux_site_minus
  end where

  ! normalize by volume of fluid phase
  if( time > 0 ) then
    n_fluidsites = sum(inside,mask=inside==fluid)
    el_curr_x = el_curr_x / n_fluidsites
    el_curr_y = el_curr_y / n_fluidsites
    el_curr_z = el_curr_z / n_fluidsites
    ion_curr_x = ion_curr_x / n_fluidsites
    ion_curr_y = ion_curr_y / n_fluidsites
    ion_curr_z = ion_curr_z / n_fluidsites
  end if

  ! check that the sum of all fluxes is zero
  if( abs(sum( flux_site_plus,mask=inside==fluid)) > 1.0e-12 .or. abs(sum( flux_site_minus,mask=inside==fluid)) > 1.0e-12 ) then
    stop 'in smolu.f90 the sum of all fluxes does not add up to zero! stop.'
  end if

end subroutine smolu
