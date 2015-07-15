PROGRAM main

  use mod_time, only: tick, tock
  use io, ONLY: print_tail

  IMPLICIT NONE
  
  INTEGER :: t
  character(8)  :: date
  character(10) :: time

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

  CALL drop_tracers

  CALL print_tail
  print*,"Execution time:", real(tock(t),4),"sec"
END PROGRAM main
