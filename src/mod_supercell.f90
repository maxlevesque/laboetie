! where one defines everything related to the supercell. Not the fluid or anything else. For now, does not contain
! its most evident constants lx ly and lz which are the supercell sizes (in LB units). TODO
module look_at_supercell
    use precision_kinds
    use system, only: supercell, fluid, solid
    use constants, only: x, y, z
    implicit none
    private
    public :: check_that_at_least_one_node_is_fluid,&
              check_that_all_nodes_are_wether_fluid_or_solid
    contains
        subroutine check_that_at_least_one_node_is_fluid
            if(count(supercell%node%nature==fluid)==0) &
                stop 'STOP subroutine check_that_at_least_one_node_is_fluid identified no fluid node. fatal error.'
        end subroutine
        subroutine check_that_all_nodes_are_wether_fluid_or_solid
            integer(i2b) :: volume
            volume = product( supercell%geometry%dimensions%indiceMax(:) )
            if((count(supercell%node%nature==fluid)+count(supercell%node%nature==solid)) /= volume ) then
                print*,count(supercell%node%nature==fluid)+count(supercell%node%nature==solid), volume
                stop 'STOP some nodes are neither solid nor fluid.'
            end if
        end subroutine
end module
