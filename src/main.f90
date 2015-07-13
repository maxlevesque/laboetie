PROGRAM main

  use mod_time, only: tick, tock

  IMPLICIT NONE
  
  INTEGER :: t

  CALL tick(t)

  !
  ! init the system: read input, geometry, composition, external forces, periodicity, ...
  !
  CALL init_simu

  !
  ! system equilibration: -D_equil <= time <= 0
  ! one solves coupled Poisson and Nernst-Planck equations without solvant flux nor external forces
  ! one then arrives at Poisson-Boltzmann distribution. For now, Poisson is solved using SOR and
  ! Nernst-Planck by Link-Flux without advection.
  ! THIS DONT WORK RIGHT NOW
  !
  CALL poisson_nernst_planck


  CALL equilibration
stop "entre equilibration_new et drop_tracers"
  PRINT*
  PRINT*,'Moment propagation'
  PRINT*,'=================='
  CALL drop_tracers
  PRINT*

  print*,"Execution time:", real(tock(t),4),"sec"
END PROGRAM main
