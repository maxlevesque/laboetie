! Here we define everything related to the LaBo code. Should be done in separate modules.
! By the way, that's a lot of work and I'm not sure how pertinent it is.
! For now we're stuck to D3Q19

module system
    use precision_kinds
    use constants, only: x, y, z
    use mod_lbmodel, only: lbm
    implicit none  
    ! time
    integer(i2b) :: time ! time in simu, from 0 to tmax
    integer(i2b) :: t_equil, tmom, tmax ! 0->t_equil; t_equil->tmom; t_mom->tmax
    ! equilibration
    integer(i2b) :: D_equil ! number of equilibration steps whithout constraints nor flux
    integer(i2b) :: D_iter
    ! supercell geometry related
    integer(i2b) :: wall
    integer(i2b), public :: lx, ly, lz
    integer(i0b), parameter, public :: fluid=0, liquid=0, solid=1

    ! fluid related
    type type_latticenode
        integer(kind(fluid)) :: nature ! solid liquid; replace latticenode(i,j,k)%nature = inside(i,j,k)
        real(dp), dimension(x:z) :: normal ! vector normal to interface if interfacial site
    end type
    type type_supercell
        integer(i2b), dimension(x:z) :: length, lengthmin, lengthmax
        type (type_latticenode), public, allocatable, dimension(:,:,:) :: node
    end type
    type (type_supercell), public :: supercell ! supercell%node(i,j,k)%nature  supercell%node(i,j,k)%normal(x)   supercell%length   supercell%lengthmin:supercell%lengthmax
    integer(i0b), allocatable, dimension(:,:,:) :: inside ! fluid or solid at each node
    !  real(dp), allocatable, dimension(:,:,:,:) :: normal_c
    real(dp), allocatable, dimension(:,:,:,:) :: normal ! nfft1, nfft2, nfft3, 3
    real(dp) :: rho_0 ! solvent density in the bulk
    real(dp), allocatable, dimension(:,:,:,:) :: n ! population :(lx,ly,lz,velocities)
    real(dp), allocatable, dimension(:,:,:) :: rho ! rho(i,j,k) = sum_l n(i,j,k,l)
    real(dp), allocatable, dimension(:,:,:) :: jx, jy, jz ! jx(i,j,k) = sum_l c_x(l) * n(i,j,k,l)    where c_x(l) = c(1,l)

    ! moments propagation
    real(dp), allocatable, dimension(:,:,:) :: flux_site_plus, flux_site_minus
    real(dp) :: anormf0

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
    real(dp), dimension(x:z) :: elec_slope, lncb_slope
    real(dp) :: el_curr_x, el_curr_y, el_curr_z
    real(dp) :: ion_curr_x, ion_curr_y, ion_curr_z
    real(dp) :: rho_ch ! charge density?
    ! external force given by user
    real(dp), dimension(x:z) :: f_ext ! external force (constraints) applied to flux
    real(dp), allocatable, dimension(:,:,:,:) :: solute_force ! (lx, ly, lz, 3)
    ! thermo
    real(dp), parameter :: kBT = 1.0_dp/3.0_dp ! boltzmann constant * temperature
    real(dp), parameter :: beta = 1.0_dp/kBT

    contains
    ! periodic boundary conditions
    pure function pbc (i,direction)
        integer(i2b) :: pbc
        integer(i2b), intent(in) :: i, direction
        select case (direction)
            case (x)
                if (i==0) then
                    pbc = lx
                else if (i==lx+1) then
                    pbc = 1
                else
                    pbc = i
                end if
            case (y)
                if (i==0) then
                    pbc = ly
                else if (i==ly+1) then
                    pbc = 1
                else
                    pbc = i
                end if
            case (z)
                if (i==0) then
                    pbc = lz
                else if (i==lz+1) then
                    pbc = 1
                else
                    pbc = i
                end if
        end select
    end function
end module system
