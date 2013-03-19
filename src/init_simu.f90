subroutine init_simu

  implicit none

  ! print header
  call print_header

  ! check that input, output folder and file ./lb.in exist
block
  logical :: file_exists
  inquire(file="./output/.", exist=file_exists)
  if( .not. file_exists) then
    print*,"./output folder does not exist. I create one."
    ! the system should be linux TODO TEST LINUX OR MAC DERIVATIVE HERE
    call system('mkdir -p output') ! -p do not print error if exist and create parent directory if needed
  end if
  inquire(file="./lb.in", exist=file_exists)
  if( .not. file_exists) stop "./lb.in does not exist"
end block

  ! go read the input file(s) and put everything in a character array
  call put_input_in_character_array

  ! Print input parameters found to output folder
  call print_input_in_output_folder

  ! prepare supercell geometry
  call supercell_definition ! porousmedia, inside_model, lx, ly, lz

  ! schedule simulation
  call scheduler ! t_equil, tmom, tmax

  ! init moments for Lattice Boltzmann
  call init_moments_for_LB

  ! init charge distribution
  call charges_init

end subroutine init_simu







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_moments_for_LB
  use precision_kinds, only: dp
  use system, only: rho_0, jx, jy, jz, lx, ly, lz, rho, n, NbVel
  use constants, only: x, z
  use input, only: input_dp
  implicit none
  rho_0 = input_dp('rho_0') ! read the initial, homogeneous, solvent density in input file
  allocate( rho(lx,ly,lz), source=rho_0 ) ! the density will evolve with time. It is initiated with rho_0 from input file
  if (.not. allocated(n)) allocate(n(lx,ly,lz,NbVel), source=0.0_dp)  ! zeroth order moment == population(r,v) == mass density
end subroutine init_moments_for_LB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine scheduler
  use precision_kinds, only: i2b
  use system, only: t_equil, tmom, tmax, D_iter, time
  use input, only: input_int
  implicit none

  ! 4 times are important :
  ! - 0 at which simulation starts
  ! - t_equil which is the time of equilibration ending
  ! - tmom at which we're looking at tracer moment propagation
  ! - tmax at which simulation stops

  ! init to non-physical value catched later in order to be sure they are modified
  time = 0

  D_iter = input_int('D_iter')
  tmax = input_int('tmax')
  tmom = input_int('tmom')
  t_equil = input_int('t_equil')

  ! check coherence
  if( tmax <= 0 .or. tmom <= 0 .or. t_equil <= 0 ) then
    stop 'in scheduler. no time should be negative or zero'
  end if

  ! check tmax is last 
  if( tmom > tmax .or. t_equil > tmax ) then
    stop 'equilibration and moment propagation cannot start after simulation end. check input.'
  end if

  ! tracer moment preparation should come after equilibration
  if( tmom < t_equil ) then
    stop 'tracer moment propagation should come after equilibration.tmom should be >= t_equil. check input file.'
  end if

end subroutine scheduler
