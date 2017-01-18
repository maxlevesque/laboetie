! Here we define the supercell: where are solid nodes and where are liquid nodes.
subroutine supercell_definition
  use precision_kinds!, only: i2b, dp
  use constants, only: x, y, z
  use system, only: fluid, solid, pbc, supercell, node
  use look_at_supercell, only: check_that_at_least_one_node_is_fluid,&
                       check_that_all_nodes_are_wether_fluid_or_solid
  use input, only: input_int
  use geometry, only: construct_slit, construct_cylinder, construct_cc, construct_disc_benichou,&
                      construct_sinusoidal_walls_2d, CONSTRUCT_PLANES_WITH_VARIOUS_RADIUS_2D,&
                      CONSTRUCT_TUBE_WITH_VARYING_DIAMETER, CONSTRUCT_SPHERICAL_CAVITY,&
                      CONSTRUCT_SPHERE_BENICHOU,&
                      CONSTRUCT_XUDONG_VINCENT_MARIE_CYL_BETWEEN_WALLS,&
                      construct_custom
  use io, only: print_supercell_xsf
  use mod_lbmodel, only: lbm
  implicit none

  integer ::  nb_node_inter, nb_node_solid,nb_node_fluid
  integer(i2b) :: lx, ly, lz

  ! geometry
  supercell%geometry%label = input_int("geometryLabel")

  ! define grid
  lx = input_int("lx")
  ly = input_int("ly")
  lz = input_int("lz")
  supercell%geometry%dimensions%indiceMax(x) = lx
  supercell%geometry%dimensions%indiceMax(y) = ly
  supercell%geometry%dimensions%indiceMax(z) = lz
  supercell%geometry%dimensions%indiceMin(x:z) = [1, 1, 1] ! Fortran style. C style arrays would start at 0.

  ! define which nodes are fluid and solid
  ! begins with fluid everywhere. Remember one defined fluid=0 and solid=1 as parameters.
  allocate( node(lx,ly,lz) )
  node%nature = fluid

  ! construct medium geometry
  select case (supercell%geometry%label)
  case (0) ! custom geometry
    call construct_custom
  case (1) ! supercell%geometry%label = 1 is two solid walls normal to Z axis.
    call construct_slit
  case (2) ! supercell%geometry%label = 2 => cylinder around Z axis.
    call construct_cylinder
  case (3) ! supercell%geometry%label = 3 => solid spheres in cubic face centred unit cell with at contact
    call construct_cc
  case (4) ! supercell%geometry%label = 4 => benichou disc with exists
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
    stop 'supercell%geometry%label tag in input file is invalid'
  end select

  call detectInterfacialNodes

  !Write the number of fluid, solid and interfacial node in file 'Nb_node_type.dat'
  call count_node_type

  call print_supercell_xsf

!  call defineNormalToSurface

  ! check that at least one node (!!) is of fluid type
  call check_that_at_least_one_node_is_fluid

  ! check that the sum of all nodes is the sum of fluid nodes and of solid nodes, ie that no node has been forgotten somewhere
  call check_that_all_nodes_are_wether_fluid_or_solid

  ! give a table that tells if you're interfacial or not
!  call where_is_it_fluid_and_interfacial

    contains
    ! here we define the normal to the surface. It is not obvious to me (Maximilien Levesque) what it is useful for. I even suspect it to be tricky if you're a single node vacancy, for instance, where normal will be zero but you're still in an interfacial node.
!        subroutine defineNormalToSurface
!            integer(i2b) :: i, j, k, l, iNext, jNext, kNext
!            node%normal(x) = 0._dp ! initialize to zero
!            node%normal(y) = 0._dp
!            node%normal(z) = 0._dp
!            do concurrent( i=1:lx, j=1:ly, k=1:lz, l=lbm%lmin:lbm%lmax )
!                iNext = pbc( i + lbm%vel(l)%coo(x), x)
!                jNext = pbc( j + lbm%vel(l)%coo(y), y)
!                kNext = pbc( k + lbm%vel(l)%coo(z), z)
!                node(i,j,k)%normal(:) = node(i,j,k)%normal(:)  &
!                   - lbm%vel(l)%a1 * lbm%vel(l)%coo(:) * (node(iNext,jNext,kNext)%nature - node(i,j,k)%nature)
!                if (any(node(i,j,k)%normal(:)/=0.0_dp)) then
!                    node(i,j,k)%normal(:) = node(i,j,k)%normal(:)/norm2(node(i,j,k)%normal(:))
!                end if
!            end do
!        end subroutine
        subroutine detectInterfacialNodes
        ! interfaces are here defined when one of two neighboring nodes is fluid when the other is solid
            integer(i2b) :: i, j, k, l, iNext, jNext, kNext
            node%isInterfacial = .false. ! init
            do concurrent( i=1:lx, j=1:ly, k=1:lz )
                velocityloop: do l = lbm%lmin+1, lbm%lmax ! lbm%lmin is zero velocity. arrival and departure cannot be different in nature
                    iNext = pbc( i + lbm%vel(l)%coo(x), x)
                    jNext = pbc( j + lbm%vel(l)%coo(y), y)
                    kNext = pbc( k + lbm%vel(l)%coo(z), z)
                    if( node(i,j,k)%nature /= node(iNext,jNext,kNext)%nature ) then
                        node(i,j,k)%isInterfacial = .true.
                        exit velocityloop
                    end if
                end do velocityloop
            end do
        end subroutine

    subroutine count_node_type     
            implicit none 
            integer :: i,j,k
            double precision :: nb_node_inter_dble, nb_node_fluid_dble
            nb_node_inter = 0
            nb_node_solid = 0
            nb_node_fluid = 0
        
            do concurrent(i=1:lx, j=1:ly, k=1:lz)
                if (node(i,j,k)%nature == solid) then
                    nb_node_solid = nb_node_solid + 1
                elseif (node(i,j,k)%nature == fluid) then
                    nb_node_fluid = nb_node_fluid + 1
                    if (node(i,j,k)%isInterfacial .eqv. .true.) then
                        nb_node_inter = nb_node_inter + 1
                    end if
                end if
            end do 
            nb_node_inter_dble = DBLE(nb_node_inter)
            nb_node_fluid_dble = DBLE(nb_node_fluid)
        
            open(unit=1, file='output/Nb_node_type.dat')
                write(1,*) '&INPUT_NML'
                write(1,*) 'Nn = ', lx*ly*lz, '!Total number of nodes'
                write(1,*) 'Ns = ', nb_node_solid, '!Number of solid nodes'
                write(1,*) 'Nf = ', nb_node_fluid, '!Number of fluid nodes'
                write(1,*) 'Nads = ', nb_node_inter, '!Number of interfacial nodes'
                write(1,*) 'ratio = ', nb_node_inter_dble/nb_node_fluid_dble, '!Number of interfacial nodes/Number of fluid nodes'
                write(1,*) '/'
            close(1)

        end subroutine
end subroutine supercell_definition
