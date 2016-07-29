module myAllocations
    use precision_kinds
    use system, only: supercell
    use mod_lbmodel, only: lbm
    use constants, only: x, y, z
    implicit none
    contains
        subroutine allocateReal3D (array)
            real(dp), dimension(:,:,:), allocatable, intent(inout) :: array
            integer(i2b) :: imin, imax, jmin, jmax, kmin, kmax
            if( allocated(array) ) return
            imin = supercell%geometry%dimensions%indiceMin(x)
            jmin = supercell%geometry%dimensions%indiceMin(y)
            kmin = supercell%geometry%dimensions%indiceMin(z)
            imax = supercell%geometry%dimensions%indiceMax(x)
            jmax = supercell%geometry%dimensions%indiceMax(y)
            kmax = supercell%geometry%dimensions%indiceMax(z)
            allocate( array( imin:imax, jmin:jmax, kmin:kmax), source=0._dp )
        end subroutine
        subroutine allocateReal4D (array)
            real(dp), dimension(:,:,:,:), allocatable, intent(inout) :: array
            integer(i2b) :: imin, imax, jmin, jmax, kmin, kmax, lmin, lmax
            if( allocated(array) ) return
            imin = supercell%geometry%dimensions%indiceMin(x)
            jmin = supercell%geometry%dimensions%indiceMin(y)
            kmin = supercell%geometry%dimensions%indiceMin(z)
            imax = supercell%geometry%dimensions%indiceMax(x)
            jmax = supercell%geometry%dimensions%indiceMax(y)
            kmax = supercell%geometry%dimensions%indiceMax(z)
            lmin = lbm%lmin
            lmax = lbm%lmax
            allocate( array( imin:imax, jmin:jmax, kmin:kmax, lmin:lmax), source=0._dp )
        end subroutine
end module
