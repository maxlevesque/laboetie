PROGRAM main

  use precision_kinds, only: dp
  use mod_time, only: tick, tock
  use io, ONLY: print_tail
  use module_input, only: getinput
  use system, only: supercell
  use module_equilibration, only: equilibration
  use module_tracers, only: drop_tracers

  implicit none
  
  integer :: t, lx, ly, lz
  character(8)  :: date
  character(10) :: time
  logical :: RestartPNP
  real(dp), allocatable, dimension(:,:,:) :: jx, jy, jz

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

  lx = supercell%geometry%dimensions%indiceMax(1)
  ly = supercell%geometry%dimensions%indiceMax(2)
  lz = supercell%geometry%dimensions%indiceMax(3)
  allocate( jx(lx, ly, lz), source=0._dp)
  allocate( jy(lx, ly, lz), source=0._dp)
  allocate( jz(lx, ly, lz), source=0._dp)
  
  CALL equilibration( jx, jy, jz)

  CALL drop_tracers( jx, jy, jz)

  CALL print_tail
  print*,"Execution time:", real(tock(t),4),"sec"
END PROGRAM main
