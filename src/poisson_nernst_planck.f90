! This is the first part of the LaBo code. Here, no forces are applied on solutes,
! and no flux is allowed. The objective here is to find the initial equilibrium
! distribution of the charges.

SUBROUTINE poisson_nernst_planck

  USE precision_kinds, ONLY: i2b, dp
  USE system, ONLY: D_equil, sigma, time, elec_slope, lncb_slope
  USE input
  USE io, ONLY: print_everything_related_to_charge_equil

  IMPLICIT NONE
  INTEGER(i2b) :: timestep, timestepmax
  LOGICAL :: is_converged = .FALSE.

  sigma = input_dp('sigma')
  IF (sigma==0.0_dp) PRINT*,'salt free system.' ! sigma has already been read in charge initialisation

  ! I THINK it remains important to go through all next steps in order to insure the reading of everything, etc.

  ! read the number of iterations one does for the first step of equilibration (D_iter)

  IF (D_equil<0) STOP "D_equil should be >= 0"

  CALL get_timestepmax (timestepmax)

  ! read electrostatic related stuff
  lncb_slope = input_dp3("lncb_slope")
  elec_slope = input_dp3("elec_slope")

  ! iterate until charges are equilibrated
    timestep = 0
    DO WHILE ((.NOT. is_converged) .AND. (timestep<=timestepmax))
        CALL backup_phi_c_plus_c_minus ! backup potential and solute concentrations from last step
        CALL sor ! TODO    ! compute phi with the Successive Overrelation Routine (SOR)
        CALL just_eq_smolu ! solve smoluchowski (diffusion + electrostatic part) ie not a full smolu
        ! monitor evolution of phi, c_plus, c_minus w.r.t. the backup at the beginning of the iteration
        ! this is done only every 10 loops in order not to waste too much time. arbitrary number.
        IF ((time/=-timestepmax .and. modulo(time,10)==0) .or. sigma==0.0_dp ) THEN
            CALL check_charge_distribution_equilibrium (time, is_converged)
        END IF
    END DO

    call charge_test ! check charge conservation ! TODO rename to call check_charge_conservation ! check charge conservation every 1000 steps (arbitrary number)
    if( sigma == 0.0_dp ) return
    call print_everything_related_to_charge_equil

    IF (is_converged) THEN
        PRINT*,'Convergence found at step ',time-1,' after',timestep,' steps'
    ELSE
        STOP 'Equilibrium distribution of salts not found'
    END IF



    CONTAINS

        SUBROUTINE get_timestepmax (a)
            INTEGER(i2b) :: a
            a = input_int("timestepmax_for_PoissonNernstPlanck")
            IF (a<=0) STOP "timestepmax_for_PoissonNernstPlanck must be >=1"
        END SUBROUTINE



end subroutine poisson_nernst_planck
