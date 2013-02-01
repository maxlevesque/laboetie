subroutine equilibration_with_constraints

  use precision_kinds, only: i2b, dp
  use system, only: tmom, D_iter, t_equil, time, jx, jy, jz, rho, inside, fluid, sigma, lx, ly, lz
  use populations, only: calc_n
  use charges, only: charge_conservation

  implicit none
  integer(i2b) :: i
  integer(i2b) :: fluid_nodes

  fluid_nodes = count(inside==fluid)

  print*,'       step       current(x)                current(y)                 current(z)          density(debug purp.)'
  print*,'       --------------------------------------------------------------------------------------------------------'

  ! read external forces
  call read_f_ext

  ! continue time
  timeloop: do time = t_equil, tmom

    if( modulo(time, 10000) == 0) &
      print*,time,&
             sum(jx/rho, mask=inside==fluid)/fluid_nodes, &
             sum(jy/rho, mask=inside==fluid)/fluid_nodes, &
             sum(jz/rho, mask=inside==fluid)/fluid_nodes, &
             sum(rho)/(lx*ly*lz)

    ! populations
    call calc_n

    ! print velocity profiles
    if( modulo(time, 10000) == 0) call velocity_profiles( time)

    ! fluid motion
    call propagation

    ! fluid density
    call comp_rho

    ! momenta
    call comp_j



    ! solute motion: advection step
    call advect

    ! solute motion: diffusion step
    if( sigma /= 0.0_dp ) then
      do i= 1, D_iter
        call sor ! TODO ! compute phi with sucessive overrelaxation method     
        call electrostatic_pot ! sum the electrostatic potential due to internal charges to the external field imposed by elec_slope_x:z
        call smolu ! Smoluchowski
        if( .not. charge_conservation() ) stop 'charge not conserved'
      end do
    end if

  end do timeloop

end subroutine equilibration_with_constraints
