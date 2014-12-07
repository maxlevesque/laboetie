SUBROUTINE equilibration_with_constraints
  USE precision_kinds, ONLY: i2b, dp
  USE system, ONLY: tmom, D_iter, t_equil, time, fluid, sigma, supercell, f_ext, node
  USE populations, ONLY: update_populations
  USE input, only: input_dp3
  USE constants, ONLY: x, y, z
  IMPLICIT NONE
  INTEGER(i2b) :: i
  INTEGER(i2b), parameter :: fluid_nodes
  fluid_nodes = COUNT( node%nature==fluid )
  PRINT*,'       step       current(x)                current(y)                 current(z)          density(debug purp.)'
  PRINT*,'       --------------------------------------------------------------------------------------------------------'
  f_ext = input_dp3("f_ext") ! read external forces
  DO time = t_equil, tmom
    IF( MODULO(time, 10000) == 0) THEN
      PRINT*,time,&
        sum(node%solventFlux(x)/node%solventDensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventFlux(y)/node%solventDensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventFlux(z)/node%solventDensity, mask=node%nature==fluid)/fluid_nodes, &
        sum(node%solventDensity) / real(product(supercell%geometry%dimensions%indiceMax(:)))
    END IF
    CALL update_populations ! populations
    IF( MODULO(time, 10000) == 0) CALL velocity_profiles( time) ! print velocity profiles
    CALL propagation    ! fluid motion
    CALL comp_rho    ! fluid density
    CALL comp_j    ! momenta
    CALL advect    ! solute motion: advection step
    ! solute motion: diffusion step
    IF( abs(sigma) > epsilon(1._dp) ) THEN
        DO i= 1, D_iter
            CALL sor                ! compute phi with sucessive overrelaxation method
            CALL electrostatic_pot  ! sum the electrostatic potential due to internal charges to the external field imposed by elec_slope(x:z))
            CALL smolu              ! Smoluchowski
            CALL charge_test        ! make sure charge is kept constant during simulation
        END DO
    END IF
  END DO
END SUBROUTINE equilibration_with_constraints
