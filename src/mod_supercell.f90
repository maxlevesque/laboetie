! where one defines everything related to the supercell. Not the fluid or anything else. For now, does not contain
! its most evident constants lx ly and lz which are the supercell sizes (in LB units). TODO

MODULE look_at_supercell
  use precision_kinds
  use system, only: supercell
  implicit none
  private ! everything private by default
  public :: is_interfacial,&
            where_is_it_fluid_and_interfacial,&
            check_that_at_least_one_node_is_fluid,&
            check_that_all_nodes_are_wether_fluid_or_solid
  logical, allocatable, dimension(:,:,:) :: is_interfacial ! true if the node is fluid .AND. at an interface
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE where_is_it_fluid_and_interfacial
  use system, only: lx, ly, lz, fluid, normal
  allocate( is_interfacial(lx,ly,lz) )
  where(norm2(normal,4)/=0.0_dp .and. supercell%node%nature==fluid) ! where the norm of normal(i,j,k,:)
    is_interfacial = .true.
  else where
    is_interfacial = .false.
  end where
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_that_at_least_one_node_is_fluid
  use system, only: fluid
  if(count(supercell%node%nature==fluid)==0) &
      stop 'STOP subroutine check_that_at_least_one_node_is_fluid identified no fluid node. fatal error.'
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_that_all_nodes_are_wether_fluid_or_solid
  use system, only: fluid, solid, lx, ly, lz
  if((count(supercell%node%nature==fluid)+count(supercell%node%nature==solid)) /= lx*ly*lz) stop &
       'subroutine check_that_all_nodes_are_wether_fluid_or_solidsome identified some nodes as neither solid nor fluid.'
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE
