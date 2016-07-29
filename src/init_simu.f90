module module_init_simu
    implicit none
    private
    public :: init_simu
contains

SUBROUTINE init_simu(n)

    USE precision_kinds, only: dp
    USE mod_lbmodel, ONLY: init_everything_related_to_lb_model => initialize
    USE system, ONLY: node, solid
    USE io, ONLY: print_header, print_input_in_output_folder, inquireNecessaryFilesExistence
    USE myallocations
    use module_input, ONLY: getinput
    IMPLICIT NONE
    real(dp), allocatable, intent(inout) :: n(:,:,:,:) ! xyz;l
    integer :: l

    !
    ! laboetie doesnt work for charged solutes
    !
    IF( ABS(getinput%dp('sigma', defaultvalue=0._dp)) > 0._dp ) THEN
        print*,"ERROR: laboetie can only consider uncharged systems."
        print*,"===== Dont tell Benjamin you'd like to see such feature in Laboetie :)"
        print*,"Hi Benjamin. I'm sure it is you testing this! grrrr :))"
        ERROR STOP
      END IF

    CALL print_header
    CALL inquireNecessaryFilesExistence  ! check that input, output folder and file ./lb.in exist
    CALL init_everything_related_to_lb_model ! init everything related to D3Q15 or D3Q19 etc ie LB models
    CALL supercell_definition    ! prepare supercell geometry
    CALL scheduler ! t_equil, tmom, tmax! schedule simulation

    !
    ! Init solvent populations
    !
    if(allocated(n)) then
        error stop "n is already alloc in init_simu"
    else
        if(lbm%lmin/=1) error stop "lmin is not 1 in init_simu"
        associate( nx => getinput%int("lx", assert='>0') ,&
                   ny => getinput%int("ly", assert='>0') ,&
                   nz => getinput%int("lz", assert='>0') )
        allocate( n(nx,ny,nz,lbm%lmax), source=0._dp )
        end associate
    end if

    !
    ! Init Solvent density
    ! It is 0. at the solid nodes where solvent cannot penetrate.
    !
    block
        real(dp) :: initialSolventDensity
        initialSolventDensity = getinput%dp("initialSolventDensity", 1._dp, assert=">=0")
        WHERE( node%nature /= solid )
            node%solventdensity = initialSolventDensity
        ELSEWHERE
            node%solventdensity = 0
        END WHERE
    end block
    do l= lbm%lmin, lbm%lmax
        n(:,:,:,l) = node%solventdensity * lbm%vel(l)%a0
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
        D_iter = getinput%int('D_iter',-1)
        tmax = getinput%int('tmax',-1)
        tmom = getinput%int('tmom',-1)
        t_equil = getinput%int('t_equil',-1)

end subroutine
    !
    !
    !
end subroutine
end module module_init_simu
