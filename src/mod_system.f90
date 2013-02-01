! Here we define everything related to the LaBo code. Should be done in separate modules.
! By the way, that's a lot of work and I'm not sure how pertinent it is.
! For now we're stuck to D3Q19

module system

  use precision_kinds, only: i2b, dp
  use constants, only : x, y, z
  implicit none


  ! time
  integer(i2b) :: time ! time in simu, from 0 to tmax
  integer(i2b) :: t_equil, tmom, tmax ! 0->t_equil; t_equil->tmom; t_mom->tmax

  ! equilibration
  integer(i2b) :: D_equil ! number of equilibration steps whithout constraints nor flux
  integer(i2b) :: D_iter


  ! supercell geometry related
  integer(i2b) :: PorousMedia
  integer(i2b) :: wall
  integer(i2b) :: inside_model
  integer(i2b) :: lx, ly, lz
  integer(i2b), parameter :: fluid=0, solid=1
  integer(i2b), parameter :: NbVel=19 ! D3Q19 lattice

  ! fluid related
  integer(i2b), allocatable, dimension(:,:,:) :: inside ! fluid or solid at each node
!  real(dp), allocatable, dimension(:,:,:,:) :: normal_c ! TODO REMOVE CERTAINLY. DONT KNOW WHAT IS FOR.
  real(dp), allocatable, dimension(:,:,:,:) :: normal ! nfft1, nfft2, nfft3, 3
  real(dp) :: rho_0 ! solvent density in the bulk
  real(dp), allocatable, dimension(:,:,:,:) :: n ! population :(lx,ly,lz,velocities)
  real(dp), allocatable, dimension(:,:,:) :: rho ! rho(i,j,k) = sum_l n(i,j,k,l)
  real(dp), allocatable, dimension(:,:,:) :: jx, jy, jz ! jx(i,j,k) = sum_l c_x(l) * n(i,j,k,l)    where c_x(l) = c(1,l)


  ! moments propagation
  real(dp), allocatable, dimension(:,:,:) :: flux_site_plus, flux_site_minus

  ! integration quadrature related
  real(dp), private, parameter :: a_00 = 1.0_dp/3.0_dp ! should be calculated in the code. Easy TODO. NO MAGIC NUMBERS ALLOWED
  real(dp), private, parameter :: a_01 = 1.0_dp/18.0_dp
  real(dp), private, parameter :: a_02 = 1.0_dp/36.0_dp
  real(dp), private, parameter :: a_10 = 0.0_dp
  real(dp), private, parameter :: a_11 = 1.0_dp/6.0_dp
  real(dp), private, parameter :: a_12 = 1.0_dp/12.0_dp

  integer(i2b), dimension(NbVel), parameter :: vel_inv &
                    = (/ 1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14, 17, 16, 19, 18 /) ! should be calculated in the code. Easy TODO. NO MAGIC NUMBERS ALLOWED

  integer(i2b), dimension( x:z, NbVel), parameter :: c = & ! c is one of the worst names i've ever seen ... change it as soon as possible
      transpose(reshape( source=(/ 0, 1,-1, 0, 0, 0, 0, 1, -1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0 ,&
                                   0, 0, 0, 1,-1, 0, 0, 1, -1,-1, 1, 0, 0, 0, 0, 1,-1, 1,-1 ,&
                                   0, 0, 0, 0, 0, 1,-1, 0,  0, 0, 0, 1,-1,-1, 1, 1,-1,-1, 1 /),&
                shape= (/size(c,2),size(c,1)/)  ) ) ! I use transpose because of the difference in row/column order between C and Fortran.

  real(dp) :: anormf0

  real(dp), dimension(NbVel), parameter :: a0 = &
                              (/ a_00, a_01, a_01, a_01, a_01, a_01, a_01, &
                                       a_02, a_02, a_02, a_02, a_02, a_02, &
                                       a_02, a_02, a_02, a_02, a_02, a_02  /)
  real(dp), dimension(NbVel), parameter :: a1 = &
                              (/ a_10, a_11, a_11, a_11, a_11, a_11, a_11, &
                                       a_12, a_12, a_12, a_12, a_12, a_12, &
                                       a_12, a_12, a_12, a_12, a_12, a_12  /)

  real(dp), private, parameter :: itself=0.0_dp, nn=1.0_dp, nnn=sqrt(2.0_dp) ! nearest and next nearest neighbour distance
  real(dp), dimension( 1:NbVel), parameter :: delta = &
                              (/ itself, nn, nn, nn, nn, nn, nn, nnn, nnn, nnn, &
                                 nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn /)

  ! electrostatic related
  real(dp) :: bjl ! Bjerum length
  real(dp) :: lambda_D ! debye length. <0 means salt free fluid
  real(dp) :: sigma ! sigma is the number of charges the user want to put in solid and surface. read from input file.
  real(dp) :: D_plus, D_minus ! diffusion coefficient of ions
  real(dp), allocatable, dimension(:,:,:) :: c_plus, c_minus ! concentrations
  real(dp), allocatable, dimension(:,:,:) :: phi ! internal potential (calculated by the sucessive overrelaxation methode)
  real(dp), allocatable, dimension(:,:,:) :: phi2 ! DONT KNOW WHAT IT IS FOR. IT IS USED IN ELEC_POT, CALLED DURING TRACER MOMENT PROPAGATION, BUT COMPLETELY USELESS
  real(dp), allocatable, dimension(:,:,:) :: phi_tot ! total electrostatic potential = phi + phi external
  real(dp), allocatable, dimension(:,:,:) :: phi_old, c_plus_old, c_minus_old ! backups of phi, c_plus and c_minus at previous steps for checking convergence
  character(len=3) :: charge_distrib
  real(dp) :: lncb_slope_x, lncb_slope_y, lncb_slope_z
  real(dp) :: elec_slope_x, elec_slope_y, elec_slope_z
  real(dp) :: el_curr_x, el_curr_y, el_curr_z
  real(dp) :: ion_curr_x, ion_curr_y, ion_curr_z
  type type_solute
    real(dp) :: density
  end type
  type(type_solute) :: solute


  ! external force given by user
  real(dp), dimension(x:z) :: f_ext ! external force (constraints) applied to flux
  real(dp), allocatable, dimension(:,:,:,:) :: solute_force ! (lx, ly, lz, 3)

  ! tracer related
  real(dp) :: D_tracer ! diffusion coefficient of tracers for moment propagation
  real(dp) :: z_tracer ! charge of tracer


  ! thermo
  real(dp), parameter :: kBT = 1.0_dp/3.0_dp ! boltzmann constant * temperature

  ! to be put somewhere else, plusx finds the right arrival site
  integer(i2b), allocatable, dimension(:) :: plusx, plusy, plusz ! THIS HAS TO BE CODED IN ANOTHER WAY AS SOON AS POSSIBLE



end module system
