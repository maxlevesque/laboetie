! Here we define the supercell: where are solid nodes and where are liquid nodes.
subroutine supercell_definition

  use precision_kinds, only: i2b, dp
  use constants, only: x, y, z
  use system, only: wall, fluid, solid, lx, ly, lz, inside, &
                    normal, a1, c, NbVel,&
                    plusx, plusy, plusz!, normal_c
  use supercell, only: where_is_it_fluid_and_interfacial,&
                       check_that_at_least_one_node_is_fluid,&
                       check_that_all_nodes_are_wether_fluid_or_solid,&
                       define_periodic_boundary_conditions,&
                       printgeometryxyz
  use input, only: input_int
  use geometry, only: construct_wall, construct_cylinder, construct_bcc, construct_fcc

  implicit none
  integer(kind=i2b) :: i, j, k, ip, jp, kp, l !dummy

  ! geometry
  wall = input_int("wall")

  ! define grid
  lx = input_int("lx")
  ly = input_int("ly")
  lz = input_int("lz")

  ! define which nodes are fluid and solid
  ! begins with fluid everywhere. Remember one defined fluid=0 and solid=1 as parameters.
  allocate( inside( 1:lx, 1:ly, 1:lz), source=fluid )

  ! construct medium geometry
  select case (wall)
  case (1) ! wall = 1 is two solid walls normal to Z axis.
    call construct_wall(inside)

  case (2) ! wall = 2 => cylinder around Z axis.
    call construct_cylinder(inside)

  case (3) ! wall = 3 => solid spheres in cubic face centred unit cell with at contact
    call construct_bcc(inside)

  case (4) ! wall = 4 => solid spheres in cubic face centred unit cell with at contact
    call construct_fcc(inside)

  end select

  ! counts number of solid and fluid nodes
!  print*, 'number of solid nodes / fluid nodes = ', count(inside==solid),' / ', count(inside==fluid)

  ! localise interface and define its normal vector
  allocate( normal( lx, ly, lz, x:z), source=0.0_dp ) ! vector normal to interface
!  allocate( normal_c( lx, ly, lz, NbVel ), source=0.0_dp )

  call define_periodic_boundary_conditions

block
  real(dp) :: norm_of_normal

  ! for each node
  do concurrent( i= 1: lx, j= 1: ly, k= 1: lz)

    ! each arrival site
    do concurrent( l= 1: NbVel )
      ip= plusx( i+ c( x, l))
      jp= plusy( j+ c( y, l))
      kp= plusz( k+ c( z, l))
      normal( i, j, k, :) = normal( i, j, k, :) - a1( l)* c( :, l)* ( inside( ip, jp, kp)- inside( i, j, k) )
    end do

    ! normalization
    norm_of_normal = norm2(normal(i,j,k,:))
    if(norm_of_normal/=0.0_dp) normal(i,j,k,:) = normal(i,j,k,:)/norm_of_normal

!    ! don't know what normal_c is used for
!    do concurrent( l= 1: NbVel )
!      normal_c( i, j, k, l) = sum( normal( i, j, k, :)* c( :, l) ) 
!    end do
 
  end do
end block

  ! check that at least one node (!!) is of fluid type
  call check_that_at_least_one_node_is_fluid

  ! check that the sum of all nodes is the sum of fluid nodes and of solid nodes, ie that no node has been forgotten somewhere
  call check_that_all_nodes_are_wether_fluid_or_solid

  ! give a table that tells if you're interfacial or not
  call where_is_it_fluid_and_interfacial

  ! print geometry in .xyz file
  call printgeometryxyz

end subroutine supercell_definition
