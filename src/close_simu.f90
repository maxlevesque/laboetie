subroutine close_simu
    use precision_kinds
    implicit none
    ! total execution time
    call timeExecution
    contains
        subroutine timeExecution
            real(sp), dimension(2) :: tarray
            real(sp) :: time
            call ETIME(tarray,time)
            print*,'User time in seconds =',tarray(1)
            print*,'System time in seconds =',tarray(2)
            print*,'Run time since start in seconds =',time
        end subroutine timeExecution
end subroutine close_simu
