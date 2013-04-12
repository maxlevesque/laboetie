subroutine equilibration_with_constraints

  use precision_kinds, only: i2b, dp
  use system, only: tmom, D_iter, t_equil, time, fluid, sigma, f_ext, supercell
  use populations, only: calc_n
  use input
  use constants, only: x, y, z
  implicit none
  integer(i2b) :: i
  integer(i2b) :: fluid_nodes

  fluid_nodes = count(supercell%node%nature==fluid)



  print*,'       step       current(x)                current(y)                 current(z)          density(debug purp.)'
  print*,'       --------------------------------------------------------------------------------------------------------'

  ! read external forces
  f_ext = input_dp3("f_ext")

  ! continue time
  timeloop: do time = t_equil, tmom

    if( modulo(time, 10000) == 0) &
      print*,time,&
             sum(supercell%node%solventFlux(x)/supercell%node%solventDensity, mask=supercell%node%nature==fluid)/fluid_nodes, &
             sum(supercell%node%solventFlux(y)/supercell%node%solventDensity, mask=supercell%node%nature==fluid)/fluid_nodes, &
             sum(supercell%node%solventFlux(z)/supercell%node%solventDensity, mask=supercell%node%nature==fluid)/fluid_nodes, &
             sum(supercell%node%solventDensity) / real(product(supercell%geometry%dimensions%indiceMax(:)))

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
        call electrostatic_pot ! sum the electrostatic potential due to internal charges to the external field imposed by elec_slope(x:z))
        call smolu ! Smoluchowski
        call charge_test ! make sure charge is kept constant during simulation
      end do
    end if

  end do timeloop

end subroutine equilibration_with_constraints
