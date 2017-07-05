PROGRAM main

  use mod_time, only: tick, tock
  use io, ONLY: print_tail
  use module_input, only: getinput

  implicit none
  
  integer :: t
  character(8)  :: date
  character(10) :: time
  logical :: RestartPNP

  CALL tick(t)

  !
  ! init the system: read input, geometry, composition, external forces, periodicity, ...
  !
  CALL init_simu

  !
  ! one solves coupled Poisson and Nernst-Planck equations without solvant flux nor external forces
  ! one then arrives at Poisson-Boltzmann distribution. For now, Poisson is solved using SOR and
  ! Nernst-Planck by Link-Flux without advection.
  ! THIS DONT WORK RIGHT NOW
  !
  RestartPNP = getinput%log("RestartPNP", .TRUE.)
  print*, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  print *, ' RestartPNP = ', RestartPNP
  print*, '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'
  if (RestartPNP) CALL poisson_nernst_planck

  CALL equilibration

  CALL drop_tracers

  CALL print_tail
  print*,"Execution time:", real(tock(t),4),"sec"
END PROGRAM main
