SUBROUTINE init_simu

    USE precision_kinds
    USE mod_lbmodel, ONLY: init_everything_related_to_lb_model => initialize
    USE io, ONLY: manageInputs => manage, print_header, print_input_in_output_folder, inquireNecessaryFilesExistence
    USE myallocations

    IMPLICIT NONE

    CALL print_header  ! print header
    CALL inquireNecessaryFilesExistence  ! check that input, output folder and file ./lb.in exist
    CALL manageInputs  ! go read the input file(s) and put everything in a character array
    CALL init_everything_related_to_lb_model ! init everything related to D3Q15 or D3Q19 etc ie LB models
    CALL supercell_definition    ! prepare supercell geometry
    CALL scheduler ! t_equil, tmom, tmax! schedule simulation
    CALL init_moments_for_LB    ! init moments for Lattice Boltzmann
    CALL charges_init    ! init charge distribution

    CONTAINS
    !
    !
    !
    SUBROUTINE init_moments_for_LB
        use precision_kinds, only: dp
        use system, only: n, node
        use input, only: input_dp
        use mod_lbmodel, only: lbm
        node%solventDensity = input_dp('initialSolventDensity',1.0_dp) ! read the initial, homogeneous, solvent density in input file
        if (.not. allocated(n)) call allocateReal4D(n)  ! zeroth order moment == population(r,v) == mass density
    END SUBROUTINE init_moments_for_LB
    !
    !
    !
    SUBROUTINE scheduler
        use system, only: t_equil, tmom, tmax, D_iter, time
        use input, only: input_int
        ! 4 times are important :
        ! - 0 at which simulation starts
        ! - t_equil which is the time of equilibration ending
        ! - tmom at which we're looking at tracer moment propagation
        ! - tmax at which simulation stops
        ! init to non-physical value catched later in order to be sure they are modified
        time = 0
        D_iter = input_int('D_iter',-1)
        tmax = input_int('tmax',-1)
        tmom = input_int('tmom',-1)
        t_equil = input_int('t_equil',-1)

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
