! Here we define the supercell: where are solid nodes and where are liquid nodes.
subroutine supercell_definition
  use precision_kinds!, only: i2b, dp
  use constants, only: x, y, z
  use system, only: wall, fluid, solid, lx, ly, lz, &
                    pbc, supercell, inside, normal
  use supercell, only: where_is_it_fluid_and_interfacial,&
                       check_that_at_least_one_node_is_fluid,&
                       check_that_all_nodes_are_wether_fluid_or_solid
  use input, only: input_int
  use geometry, only: construct_slit, construct_cylinder, construct_cc, construct_disc_benichou,&
                      construct_sinusoidal_walls_2d, CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D,&
                      CONSTRUCT_TUBE_WITH_VARYING_DIAMETER, CONSTRUCT_SPHERICAL_CAVITY,&
                      CONSTRUCT_SPHERE_BENICHOU,&
                      CONSTRUCT_XUDONG_VINCENT_MARIE_CYL_BETWEEN_WALLS
  use io, only: print_supercell_xsf
  use mod_lbmodel, only: lbm
  implicit none
  integer(kind=i2b) :: i, j, k, ip, jp, kp, l !dummy
  character(len=150) :: filename
  real(dp), dimension(lx,ly,lz) :: testarray

  ! geometry
  wall = input_int("wall")

  ! define grid
  lx = input_int("lx")
  ly = input_int("ly")
  lz = input_int("lz")

  ! define which nodes are fluid and solid
  ! begins with fluid everywhere. Remember one defined fluid=0 and solid=1 as parameters.
  allocate( inside( lx, ly, lz), source=fluid )
  allocate( supercell%node(lx,ly,lz) )

  ! construct medium geometry
  select case (wall)
  case (1) ! wall = 1 is two solid walls normal to Z axis.
    call construct_slit
  case (2) ! wall = 2 => cylinder around Z axis.
    call construct_cylinder
  case (3) ! wall = 3 => solid spheres in cubic face centred unit cell with at contact
    call construct_cc
  case (4) ! wall = 4 => benichou disc with exists
    call construct_disc_benichou
  case (5)
    call construct_sinusoidal_walls_2d
  case (6)
    call CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D
  case (7)
    call CONSTRUCT_TUBE_WITH_VARYING_DIAMETER
  case (8)
    call CONSTRUCT_SPHERICAL_CAVITY
  case (9)
    call CONSTRUCT_SPHERE_BENICHOU
  case (10)
    call CONSTRUCT_XUDONG_VINCENT_MARIE_CYL_BETWEEN_WALLS
  case default
    stop 'wall tag in input file is invalid'
  end select

  call print_supercell_xsf

  call define_normal_to_surfaces

  ! check that at least one node (!!) is of fluid type
  call check_that_at_least_one_node_is_fluid

  ! check that the sum of all nodes is the sum of fluid nodes and of solid nodes, ie that no node has been forgotten somewhere
  call check_that_all_nodes_are_wether_fluid_or_solid

  ! give a table that tells if you're interfacial or not
  call where_is_it_fluid_and_interfacial
  
  contains
    ! here we define normal, which is the normal of the surface. It is not obvious to me (Maximilien Levesque) what it is usefull for. I even suspect it to be tricky if you're a single node vacancy, for instance, where normal will be zero but you're still in an interfacial node.
      subroutine define_normal_to_surfaces
          integer(i2b) :: i, j, k, l, ip, jp, kp
          allocate( normal( lx, ly, lz, x:z), source=0.0_dp ) ! vector normal to interface
          do concurrent ( i=1:lx, j=1:ly, k=1:lz, l=lbm%lmin:lbm%lmax )
             ip = pbc (i+ lbm%vel(l)%coo(x) ,x)
             jp = pbc (j+ lbm%vel(l)%coo(y) ,y)
             kp = pbc (k+ lbm%vel(l)%coo(z) ,z)
             normal(i,j,k,:) = normal(i,j,k,:) - lbm%vel(l)%a1 * lbm%vel(l)%coo(:)*(inside(ip,jp,kp) - inside(i,j,k))
             if (any(normal(i,j,k,:)/=0.0_dp)) normal(i,j,k,:) = normal(i,j,k,:)/norm2(normal(i,j,k,:))
          end do

      end subroutine define_normal_to_surfaces
end subroutine supercell_definition
