SUBROUTINE equilibration_with_constraints
  use precision_kinds, only: i2b, dp
  use system, only: tmom, d_iter, t_equil, time, fluid, sigma, supercell, f_ext, node
  use populations, only: update_populations
  use input, only: input_dp3
  use constants, only: x, y, z
  implicit none
  integer(i2b) :: i
  integer(i2b) :: fluid_nodes
  fluid_nodes = count( node%nature==fluid )
  print*,'       step       current(x)                current(y)                 current(z)          density(debug purp.)'
  print*,'       --------------------------------------------------------------------------------------------------------'
  f_ext = input_dp3("f_ext") ! read external forces
  do time = t_equil, tmom
    if( modulo(time, 10000) == 0) then
      print*,time,&
        sum(node%solventflux(x)/node%solventdensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventflux(y)/node%solventdensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventflux(z)/node%solventdensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventdensity) / real(product(supercell%geometry%dimensions%indicemax(:)))
    end if
    call update_populations ! populations
    if( modulo(time, 10000) == 0) call velocity_profiles( time) ! print velocity profiles
    call propagation    ! fluid motion
    call comp_rho    ! fluid density
    call comp_j    ! momenta
    call advect    ! solute motion: advection step
    ! solute motion: diffusion step
    if( abs(sigma) > epsilon(1._dp) ) then
        do i= 1, d_iter
            call sor                ! compute phi with sucessive overrelaxation method
            call electrostatic_pot  ! sum the electrostatic potential due to internal charges to the external field imposed by elec_slope(x:z))
            call smolu              ! smoluchowski
            call charge_test        ! make sure charge is kept constant during simulation
        end do
    end if
  end do
end subroutine equilibration_with_constraints
