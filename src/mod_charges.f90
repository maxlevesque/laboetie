MODULE CHARGES

  use precision_kinds, only: dp, i2b
  use system, only: lambda_D, lx, ly, lz, c_plus, c_minus, phi, charge_distrib, sigma, &
                     inside, normal, fluid, solid, anormf0, bjl, kBT, solute, D_plus, D_minus
  use constants, only: pi
  use input, only: input_dp

  public :: charges_init, charge_conservation
  type lattice_nodes
    integer(i2b) :: solid, fluid, solid_int, fluid_int
    real(dp) :: sigma_solid, sigma_fluid
  end type
  type(lattice_nodes) my_lattice_nodes


! Here we initiate the charged densities (solutes + and -)

CONTAINS

subroutine charges_init
  use supercell, only: count_interfacial
  implicit none
  real(kind=dp) :: in_c_plus_solid, in_c_plus_fluid, in_c_minus_solid, in_c_minus_fluid

  ! init ion (solute) concentrations
  if (.not. allocated(c_plus)) allocate( c_plus(lx,ly,lz), source=0.0_dp )
  if (.not. allocated(c_minus)) allocate( c_minus(lx,ly,lz), source=0.0_dp )

  ! init potential
  allocate( phi(lx,ly,lz), source=0.0_dp )

  ! read
  bjl = input_dp('bjl')
  sigma = input_dp('sigma')
  lambda_d = input_dp('lambda_D')
  
  if( lambda_D <= 0.0_dp ) then
    solute%density = 0.0_dp ! salt free fluid
  else
    solute%density = 1.0_dp / (4.0_dp*pi*bjl*lambda_D**2)
  end if

  ! read where are distributed the charges
  call read_charge_distrib ! int or sol

  my_lattice_nodes%solid = count(inside==solid) ! total number of solid nodes
  my_lattice_nodes%fluid = count(inside==fluid) ! total number of liquid nodes
  my_lattice_nodes%solid_int = count_interfacial(solid)
  my_lattice_nodes%fluid_int = count_interfacial(fluid)

  ! distribute charge, on all solid nodes (sol) or interfacial solid nodes (int)
  if( charge_distrib(1:3) == 'sol') then
    my_lattice_nodes%sigma_solid = sigma / my_lattice_nodes%solid ! charges distributed in all solid nodes
  else if( charge_distrib(1:3) == 'int') then
    my_lattice_nodes%sigma_solid = sigma / my_lattice_nodes%solid_int ! charges distributed in interfacial solid nodes only
  end if

  my_lattice_nodes%sigma_fluid = -1.0_dp * sigma / my_lattice_nodes%fluid ! charges distributed in all fluid nodes, no other option

  in_c_plus_solid  = +0.5_dp*my_lattice_nodes%sigma_solid;
  in_c_minus_solid = -0.5_dp*my_lattice_nodes%sigma_solid;

  if( my_lattice_nodes%sigma_fluid > 0.0_dp ) then
    in_c_plus_fluid  = 0.5_dp*solute%density + my_lattice_nodes%sigma_fluid
    in_c_minus_fluid = 0.5_dp*solute%density
  else
    in_c_plus_fluid  = 0.5_dp*solute%density
    in_c_minus_fluid = 0.5_dp*solute%density - my_lattice_nodes%sigma_fluid
  end if

  print*,'Salt concentration ',0.5*solute%density
  print*,'Init density values :'
  print*,'p_solid =',in_c_plus_solid
  print*,'p_fluid =',in_c_plus_fluid
  print*,'m_solid =',in_c_minus_solid
  print*,'m_fluid =',in_c_minus_fluid
  print*,'*********************************************************************'

! compute solute concentration in fluid
  where( inside == fluid )
    c_plus = in_c_plus_fluid
    c_minus = in_c_minus_fluid
  end where

! compute solute concentration in solid
  call distrib_charges(charge_distrib(1:3), in_c_plus_solid, in_c_minus_solid)


  anormf0 = 4.0_dp*pi*bjl*kBT/2.0_dp *sum(c_plus+c_minus)

  if( .not. charge_conservation() ) stop 'charge not conserved l.09'

  ! read diffusion coefficients of solute + and solute -
  D_plus = input_dp('D_plus')
  if( D_plus < 0.0_dp ) stop 'D_plus detected as negative in charges_init.f90. critical.'
  D_minus = input_dp('D_minus')
  if( D_minus < 0.0_dp ) stop 'D_minus detected as negative in charges_init.f90. critical.'

END SUBROUTINE CHARGES_INIT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! In this routine one wants to be sure that the total charges the user
! asks as input (sigma) is kept constant during the simulation
LOGICAL  FUNCTION CHARGE_CONSERVATION() ! TODO REMOVE PRINT AND TRANSFORM TO PURE FUNCTION
  implicit none
  real(dp), parameter :: tolerance = 1.0e-10 ! how much change due to numeric in charge is allowed TODO magic number
  real(dp) :: solid_charge, fluid_charge ! should be equal to sigma (asked by user)
  ! the total charge is c_plus - c_minus
  solid_charge = sum(c_plus-c_minus ,mask=(inside==solid))
  fluid_charge = sum(c_plus-c_minus ,mask=(inside==fluid))
  if( abs(solid_charge-sigma) > tolerance .or. abs(fluid_charge+sigma) > tolerance ) then
    charge_conservation = .false.
    PRINT*,'a99123 fluid charge=',fluid_charge,' solid charge=',solid_charge
    STOP
  else
    charge_conservation = .true.
    PRINT*,'a99123 fluid charge=',fluid_charge,' solid charge=',solid_charge
  end if
END FUNCTION CHARGE_CONSERVATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE DISTRIB_CHARGES(dis,cp,cm)
  use supercell, only: is_an_interfacial_node
  implicit none
  character(len=3), intent(in) :: dis
  real(dp), intent(in) :: cp, cm
  integer(i2b) :: i, j, k
  if( dis == 'int' ) then
    do concurrent( i=1:lx, j=1:ly, k=1:lz, inside(i,j,k) == solid )
      if( is_an_interfacial_node(i,j,k) ) then
        c_plus(i,j,k) = cp
        c_minus(i,j,k) = cm
      else
        c_plus(i,j,k) = 0.0_dp
        c_minus(i,j,k) = 0.0_dp
      end if
    end do
  else if( dis == 'sol' ) then
    where( inside == solid )
      c_plus = cp
      c_minus = cm
    end where
  end if
END SUBROUTINE DISTRIB_CHARGES


END MODULE CHARGES
