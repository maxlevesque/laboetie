! where one defines everything related to the supercell. Not the fluid or anything else. For now, does not contain
! its most evident constants lx ly and lz which are the supercell sizes (in LB units). TODO
MODULE supercell
  use precision_kinds
  use system, only: lx, ly, lz, inside, fluid, solid
  implicit none
  private ! everything private by default
  public :: is_interfacial,&
            where_is_it_fluid_and_interfacial,&
            check_that_at_least_one_node_is_fluid,&
            check_that_all_nodes_are_wether_fluid_or_solid,&
            define_periodic_boundary_conditions,&
            printgeometryxyz,&
            is_an_interfacial_node,&
            count_interfacial
  logical, allocatable, dimension(:,:,:) :: is_interfacial ! true if the node is fluid .AND. at an interface
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

LOGICAL PURE FUNCTION is_an_interfacial_node (i,j,k)
  use system, only: nbvel, c, plusx, plusy, plusz
  implicit none
  integer(i2b), intent(in) :: i, j, k
  integer(i2b) :: ip, jp, kp, l
  integer(i2b) :: central_node
  is_an_interfacial_node = .false. ! default value false
  central_node = inside(i,j,k)
  do l = 1, NbVel
    ip = plusx(i+c(1,l)) ; jp = plusy(j+c(2,l)) ; kp = plusz(k+c(3,l))
    if( inside(ip,jp,kp) /= central_node ) then
      is_an_interfacial_node = .true.
      exit
    end if
  end do
END FUNCTION is_an_interfacial_node 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE where_is_it_fluid_and_interfacial
  use constants, only: x, y, z
  implicit none
  integer(i2b) :: i, j, k
  allocate( is_interfacial(lx,ly,lz) )
  do k = 1, lz ; do j = 1, ly ; do i = 1, lx
    if( inside(i,j,k)==fluid .and. is_an_interfacial_node(i,j,k) ) is_interfacial(i,j,k) = .true.
  end do ; end do ; end do
END SUBROUTINE where_is_it_fluid_and_interfacial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PURE INTEGER(I2B) FUNCTION COUNT_INTERFACIAL(state)
  implicit none
  integer(i2b), intent(in) :: state
  integer(i2b) :: i, j, k
  count_interfacial = 0 ! total number of solid nodes at interfaces
  do k=1,lz ; do j=1,ly ; do i=1,lx
    if( inside(i,j,k)==state .and. is_an_interfacial_node(i,j,k) ) count_interfacial = count_interfacial + 1
  end do ; end do ; end do
END FUNCTION COUNT_INTERFACIAL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_that_at_least_one_node_is_fluid
  implicit none
  if(count(inside==fluid)==0) stop 'subroutine check_that_at_least_one_node_is_fluid identified no fluid node. fatal error.'
END SUBROUTINE check_that_at_least_one_node_is_fluid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_that_all_nodes_are_wether_fluid_or_solid
  use system, only: inside, fluid, solid, lx, ly, lz
  implicit none
  if((count(inside==fluid)+count(inside==solid)) /= lx*ly*lz) stop &
       'subroutine check_that_all_nodes_are_wether_fluid_or_solidsome identified some nodes as neither solid nor fluid.'
END SUBROUTINE check_that_all_nodes_are_wether_fluid_or_solid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE define_periodic_boundary_conditions
  use system, only: plusx, plusy, plusz
  implicit none
  integer(i2b) :: i
  allocate( plusx(0:lx+1), plusy(0:ly+1), plusz(0:lz+1) )
  do concurrent (i = 1: lx)
    plusx(i) = i
  end do
  do concurrent (i = 1: ly)
    plusy(i) = i
  end do
  do concurrent (i = 1: lz)
    plusz(i) = i
  end do
  plusx(0) = lx
  plusy(0) = ly
  plusz(0) = lz
  plusx(lx+1) = 1
  plusy(ly+1) = 1
  plusz(lz+1) = 1
END SUBROUTINE define_periodic_boundary_conditions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PRINTGEOMETRYXYZ
  ! print solid nodes as atoms in a .xyz file so that any atomic visualisation tool can make solid nodes visible
  implicit none
  integer(i2b) :: i, j, k
  open( unit=99, file='output/geometry.xyz', iostat=i)
  if( i /= 0 ) stop 'problem in init_simu.f90 to opening output/solidliquid.xyz'
  write(99,*) lx*ly*lz
  write(99,*) 'supercell from laboetie'
  do i= 1, lx ; do j= 1, ly ; do k = 1, lz
    if( inside( i, j, k) == solid ) then
      write(99,*)'C ', real(i-1,dp), real(j-1,dp), real(k-1,dp)
    else if( inside( i, j, k) == fluid .and. .not.is_an_interfacial_node(i,j,k) ) then
      write(99,*)'H ', real(i-1,dp), real(j-1,dp), real(k-1,dp)
    else if( inside( i, j, k) == fluid .and. is_an_interfacial_node(i,j,k) ) then
      write(99,*)'U ', real(i-1,dp), real(j-1,dp), real(k-1,dp)
    else
      stop 'unattended stop at subroutine supercell_definition in charges_init.f90'
    end if
  end do ; end do ; end do
  close( 99)
END SUBROUTINE PRINTGEOMETRYXYZ

END MODULE supercell
