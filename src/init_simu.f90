SUBROUTINE init_simu

    USE precision_kinds
    USE mod_lbmodel, ONLY: init_everything_related_to_lb_model => initialize, lbm
    USE system, ONLY: node, n, solid
    USE io, ONLY: print_header, print_input_in_output_folder, inquireNecessaryFilesExistence
    USE myallocations
    use module_input, ONLY: getinput

    IMPLICIT NONE

    REAL(dp) :: initialSolventDensity
    integer :: l

    CALL print_header
    CALL inquireNecessaryFilesExistence  ! check that input, output folder and file ./lb.in exist
    CALL init_everything_related_to_lb_model ! init everything related to D3Q15 or D3Q19 etc ie LB models
    CALL supercell_definition    ! prepare supercell geometry
    CALL scheduler ! t_equil, tmom, tmax! schedule simulation

    !
    ! Init solvent populations
    !
    IF( .NOT. ALLOCATED(n)) CALL allocatereal4D(n)
    n = 0

    !
    ! Init Solvent density
    ! but in the solid, where it is 0
    !
    initialSolventDensity = getinput%dp("initialSolventDensity", 1._dp)
    do l = lbm%lmin, lbm%lmax
        WHERE( node%nature /= solid )
            n(:,:,:,l) = initialSolventDensity * lbm%vel(l)%a0
        ELSEWHERE
            n(:,:,:,l) = 0
        END WHERE
    end do

    CALL charges_init    ! init charge distribution

CONTAINS
    !
    !
    !
SUBROUTINE scheduler
        use system, only: t_equil, tmom, tmax, D_iter, time
        use module_input, only: getinput
        ! 4 times are important :
        ! - 0 at which simulation starts
        ! - t_equil which is the time of equilibration ending
        ! - tmom at which we're looking at tracer moment propagation
        ! - tmax at which simulation stops
        ! init to non-physical value catched later in order to be sure they are modified
        time = 0
        D_iter = getinput%int('D_iter',1) ! ADE : the default value should be 1 (check)
        tmax = getinput%int('tmax',-1)
        tmom = getinput%int('tmom',-1)
        t_equil = getinput%int('t_equil',-1)

        !! check coherence
        !if( tmax <= 0 .or. tmom <= 0 .or. t_equil <= 0 ) then
        !    stop 'in scheduler. no time should be negative or zero'
        !end if
        !! check tmax is last
        !if( tmom > tmax .or. t_equil > tmax ) then
        !    stop 'equilibration and moment propagation cannot start after simulation end. check input.'
        !end if
        !! tracer moment preparation should come after equilibration
        !if( tmom < t_equil ) then
        !    stop 'tracer moment propagation should come after equilibration.tmom should be >= t_equil. check input file.'
        !end if
end subroutine
    !
    !
    !
end subroutine
