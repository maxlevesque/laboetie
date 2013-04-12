subroutine velocity_profiles (time)
    use precision_kinds, only : i2b, dp
    use system, only: supercell
    use constants, only : x, y, z
    implicit none
    integer(i2b), intent(in) :: time
    open(11, file= 'output/velocity_profile')
    write(11,*)
    write(11,*)'# timestep = ',time
    call twoWallsNormalToZFluxAlongX
    close(11)
    contains
        ! velocity profile v_x(z) for
        !   z
        !   | 
        !   |
        !   |======================= infinite wall along x
        !   |
        !   |       flux along x only
        !   |
        !   |---------------------------> x
        !   |
        !   |       velocity profile v_y=v_z=0, v_x function of z
        !   |
        !   |======================= infinite wall along x
        !   |
        subroutine twoWallsNormalToZFluxAlongX
            integer(i2b) :: imin, jmin, k
            imin = supercell%geometry%dimensions%indiceMin(x)
            jmin = supercell%geometry%dimensions%indiceMin(y)
            do k = supercell%geometry%dimensions%indiceMin(z), supercell%geometry%dimensions%indiceMax(z)
                write(11,*) k, supercell%node(imin,jmin,k)%solventFlux(x)/supercell%node(imin,jmin,k)%solventDensity
            end do
        end subroutine
end subroutine
