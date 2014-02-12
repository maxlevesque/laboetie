! Here tracers are droped in the fluid. They may be charged or not. They evoluate in the
! equilibrated fluid and its solutes. They do not change the potential even if they
! have a charge. The idea is to follow them in order to get the velocity auto correlation
! function, while making them not to change the equilibrium properties of the system.
! Imagine a very small droplet of radioactive particles, so few they do not change
! anything to the system, but numerous enough to be followed and make statistics.

SUBROUTINE drop_tracers

    USE precision_kinds, only: dp, i2b
    USE system, only: tmom, tmax, elec_slope
    USE populations, only: calc_n_momprop
    USE moment_propagation, only: init, propagate, deallocate_propagated_quantity!, print_vacf, integrate_vacf!, not_yet_converged
    USE input, only: input_dp
    
    IMPLICIT NONE
    INTEGER(i2b) :: it
    LOGICAL :: is_converged

    PRINT*,'       step           VACF(x)                   VACF(y)                   VACF(z)'
    PRINT*,'       ----------------------------------------------------------------------------------'

    CALL calc_n_momprop ! include elec_slope in population n
    elec_slope = 0.0_dp ! turn the electric field off for the rest of mom_prop (included in n)

    ! add electrostatic potential computed by the SOR routine an external contribution
    ! elec_pot(singlx,ly,lz, ch, phi, phi2, t, t_equil);
    ! call elec_pot
    CALL init ! init moment propagation

    ! propagate in time
    momproploop: DO it= 1, tmax-tmom
    !  it=0
    !  do while (not_yet_converged(it))
    !   call elec_pot
        CALL propagate (it,is_converged) ! propagate the propagated quantity
    !    if( modulo(it,50000)==0 ) print_vacf
    !    it = it + 1
        IF( is_converged ) EXIT momproploop
    END DO momproploop

    PRINT*,

    !  CALL integrate_vacf ! compute the integral over time of vacf in each direction
    !  CALL print_vacf ! print vacf to file vacf.dat
    CALL deallocate_propagated_quantity

END SUBROUTINE drop_tracers
