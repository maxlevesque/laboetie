PROGRAM main

  IMPLICIT NONE
  ! init the system: read input, geometry, composition, external forces, periodicity, ...
  CALL init_simu

  ! system equilibration: -D_equil <= time <= 0
  ! one solves coupled Poisson and Nernst-Planck equations without solvant flux nor external forces
  ! one then arrives at Poisson-Boltzmann distribution. For now, Poisson is solved using SOR and
  ! Nernst-Planck by Link-Flux without advection.
  PRINT*
  PRINT*,'Poisson + Nernst-Plack'
  PRINT*,'======================'
  CALL poisson_nernst_planck
  PRINT*

  PRINT*
  PRINT*,'Equilibration new'
  PRINT*,'==========================='
  CALL equilibration_new
  PRINT*

   !0 <= time <= Tequil
  !PRINT*
  !PRINT*,'Unconstrained equilibration'
  !PRINT*,'==========================='
  !CALL equilibration_without_constraint
  !PRINT*

   !Tequil < time < Tmom
   !One applies external forces and solves coupled Navier-Stokes and Nernst-Planck equations.
   !At each timestep, one calculates populations via LB using forces calculated on each node
   !at previous step and deduces local velocities. Ionic flux are then calculated using link-flux.
   !Link flux gives access to force created on fluid for following step. At each new local ionic
   !concentration, one calculates the electrostatic potential solving Poisson equation.
  !PRINT*
  !PRINT*,'Constrained equilibration'
  !PRINT*,'========================='
  !CALL equilibration_with_constraints
  !PRINT*

  ! Tmom <= time <= Tmax
  ! once stationary point is reached under constraints, populations and electrostatic potential
  ! are used in moment propagation for charged tracers. At this step, system is no more
  ! evolving via LB and link-flux: only P(r,t) of moment propagation method evolves and
  ! is used to compute the velocity autocorrelation function.
  PRINT*
  PRINT*,'Moment propagation'
  PRINT*,'=================='
  CALL drop_tracers
  PRINT*

  ! close everything nicely
  CALL close_simu

END PROGRAM main
