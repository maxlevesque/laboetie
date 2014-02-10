! Now that the charges (if present) are equilibrated, one makes things move here.
! Here we make things move but without any constraints, ie without external
! forces applied

SUBROUTINE equilibration_without_constraint

    USE precision_kinds, ONLY: i2b, dp
    USE system, ONLY: t_equil, f_ext, D_iter, time, elec_slope, lncb_slope
    USE populations, ONLY: calc_n
    USE myallocations
    USE constants, ONLY: x, y, z
    
    IMPLICIT NONE
    INTEGER(i2b) :: iteration
    
    PRINT*,'       step       current(x)                current(y)                 current(z)'
    PRINT*,'       ----------------------------------------------------------------------------------'

    f_ext = 0.0_dp
    elec_slope = 0.0_dp
    lncb_slope = 0.0_dp

    timeloop: DO time = 1, t_equil
    
        CALL calc_n      ! get populations
        CALL propagation ! propagation of the distribution functions according to their direction to the next nodes
        CALL comp_rho    ! compute density
        CALL comp_j      ! compute momenta
        CALL advect      ! solute motion: advection step

        DO iteration= 1, D_iter     ! solute motion: diffusion step
            CALL sor                ! compute phi with sucessive overrelaxation method
            CALL electrostatic_pot  ! compute phi + phi_external (due to electrostatic field elec_slope(x:z))
            CALL smolu              ! this time full smolu, not just_equ_smolu
            CALL charge_test        ! check charge conservation
        END DO
    
    END DO timeloop

END SUBROUTINE equilibration_without_constraint
