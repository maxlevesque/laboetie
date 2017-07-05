module module_bounceback

    implicit none
    private
    public bounceback    

contains

    subroutine bounceback( n, nodeNature)
        use precision_kinds, only: dp
        use system, only: lbm, pbc
        implicit none
        integer :: i, j, k, l, ip, jp, kp, lx, ly, lz, ios, lmin, lmax
        integer(1), intent(in) :: nodeNature(:,:,:)
        real(dp), intent(inout) :: n(:,:,:,:) ! population(x,y,z,l)
        real(dp) :: n_loc
        integer, allocatable :: il(:,:), jl(:,:), kl(:,:)
        lmin = lbm%lmin
        lmax = lbm%lmax
        lx = ubound(n,1)
        ly = ubound(n,2)
        lz = ubound(n,3)
        !
        ! Tabulate the index of the node one finishes if one starts from a node and a velocity index l
        ! per direction
        !
        allocate( il(lbm%lmin:lbm%lmax, 1:lx), stat=ios)
        if (ios /= 0) stop "il: Allocation request denied"
        allocate( jl(lbm%lmin:lbm%lmax, 1:ly), stat=ios)
        if (ios /= 0) stop "jl: Allocation request denied"
        allocate( kl(lbm%lmin:lbm%lmax, 1:lz), stat=ios)
        if (ios /= 0) stop "kl: Allocation request denied"
        do l= lbm%lmin, lbm%lmax
            il(l,:) = [( pbc(i+lbm%vel(l)%coo(1), 1) ,i=1,lx )]
            jl(l,:) = [( pbc(j+lbm%vel(l)%coo(2), 2) ,j=1,ly )]
            kl(l,:) = [( pbc(k+lbm%vel(l)%coo(3), 3) ,k=1,lz )]
        end do

        do l = lbm%lmin, lbm%lmax, 2
            do k = 1, lz
                kp = kl(l,k)
                do j = 1, ly
                    jp = jl(l,j)
                    do i = 1, lx
                        ip = il(l,i)
                        if( nodeNature(i,j,k) /= nodeNature(ip,jp,kp) ) then
                            n_loc = n(i,j,k,l)
                            n(i,j,k,l) = n(ip,jp,kp,lbm%vel(l)%inv)
                            n(ip,jp,kp,lbm%vel(l)%inv) = n_loc
                        end if
                    end do
                end do
            end do
        end do
        deallocate( il, jl, kl)
    end subroutine bounceback

end module module_bounceback