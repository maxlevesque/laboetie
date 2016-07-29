program main

    use precision_kinds, only: dp
    use mod_time, only: tick, tock
    use io, ONLY: print_tail
    use module_init_simu, only: init_simu
    use module_equilibration, only: equilibration
    implicit none
    INTEGER :: t

    !
    ! The heart of the lattice boltzmann method lies in the populations, n(x,y,z,l).
    ! where l is the index of the discret velocities.
    ! populations evolve with time.
    ! This variable is the most important want, the one you should follow with great attention.
    !
    real(dp), allocatable :: n(:,:,:,:)

    !
    ! Start the chronometre
    !
    call tick(t)


    !
    ! init the system: read input, geometry, composition, external forces, periodicity, ...
    !
    call init_simu(n)

    !
    ! system equilibration: -D_equil <= time <= 0
    ! one solves coupled Poisson and Nernst-Planck equations without solvant flux nor external forces
    ! one then arrives at Poisson-Boltzmann distribution. For now, Poisson is solved using SOR and
    ! Nernst-Planck by Link-Flux without advection.
    ! THIS DONT WORK RIGHT NOW
    !
    ! CALL poisson_nernst_planck

    !
    ! Make things move by applying forces
    !
    call equilibration(n)

  CALL drop_tracers

  CALL print_tail
END PROGRAM main
