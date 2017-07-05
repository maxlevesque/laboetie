module module_propagation
    implicit none
    private
    public propagation
contains
    subroutine propagation(n, lmin, lmax, lx, ly, lz, il, jl, kl)
        use precision_kinds, only: dp
        implicit none
        real(dp), intent(inout) :: n(:,:,:,:)
        integer, intent(in) :: lmin, lmax, lx, ly, lz
        integer, intent(in) :: il(:,:), jl(:,:), kl(:,:)
        real(dp) :: n_old(lx,ly,lz)
        integer :: l, k, j, i, ip, jp, kp

        !$OMP PARALLEL DO DEFAULT(NONE) &
        !$OMP SHARED(n,lz,ly,lx,lmin,lmax,il,jl,kl) &
        !$OMP PRIVATE(l,k,j,i,ip,jp,kp,n_old)
        do l=lmin,lmax
            n_old = n(:,:,:,l)
            do k=1,lz
                kp = kl(l,k)
                do j=1,ly
                    jp = jl(l,j)
                    do i=1,lx
                        ip = il(l,i)
                        n(ip,jp,kp,l) = n_old(i,j,k)
                    end do
                end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine propagation
end module module_propagation
