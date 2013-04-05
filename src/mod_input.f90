! This module contains the array that contains the input file.
! It is easy to read this input_array and modify it, in opposition to the
! input data file that one does not want to modify.
! TODO transform to a single polymorphic function, with optional arguments for the number of values to be readen
! TODO input_dp(that)+input_dp3(that) could be fusionned in a input_dp(that,hm) with hm declared as optional intent(in)
module input
    use precision_kinds
    implicit none
    character (len=100), allocatable, dimension(:), public :: input_line ! array containing all input lines
    contains
        pure function input_dp(that) ! reads hm (how many) reals. hm is optional, 
            real(dp) :: input_dp
            character(*), intent(in) :: that
            integer(i2b) :: i, j, imax
            j=len(that)
            imax = size(input_line)
            do concurrent (i=1:imax)
                if( input_line(i)(1:j) == that) read ( input_line (i) (j+4:j+50) , * ) input_dp
            end do
        end function
        pure function input_dp3(that) ! reads hm (how many) reals. hm is optional, 
            real(dp), dimension(3) :: input_dp3
            character(*), intent(in) :: that
            integer(i2b) :: i, j, imax
            j=len(that)
            imax = size(input_line)
            do concurrent (i=1:imax)
                if( input_line(i)(1:j) == that) read ( input_line (i) (j+4:j+50) , * ) input_dp3
            end do
        end function
        pure function input_int(that) ! reads an integer
            integer(i2b) :: input_int
            character(*), intent(in) :: that
            integer(i2b) :: i, j, imax
            imax = size(input_line)
            j=len(that)
            do concurrent (i=1:imax)
                if( input_line( i)( 1:j) == that ) read ( input_line (i) (j+4:j+50) , * ) input_int
            end do
        end function
        pure function input_int3(that) ! reads an integer
            integer(i2b), dimension(3) :: input_int3
            character(*), intent(in) :: that
            integer(i2b) :: i, j, imax
            imax = size(input_line)
            j=len(that)
            do concurrent (i=1:imax)
                if( input_line( i)( 1:j) == that ) read ( input_line (i) (j+4:j+50) , * ) input_int3
            end do
        end function
        pure function input_ch(that) ! reads a character
            character(50) :: input_ch
            character(*), intent(in) :: that
            integer(i2b) :: i, j, imax
            imax = size(input_line)
            j=len(that)
            do concurrent (i=1:imax)
                if( input_line( i)( 1:j) == that ) read(input_line (i) (j+4:j+50) , '(A)' ) input_ch
            end do
        end function
end module
