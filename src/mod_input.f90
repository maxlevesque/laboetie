! This module contains the array that contains the input file.
! It is easy to read this input_array and modify it, in opposition to the
! input data file that one does not want to modify.

MODULE INPUT
  use precision_kinds
  implicit none
  character (len=100), allocatable, dimension(:) :: input_line ! array containing all input lines
  private
  public :: input_line, input_dp, input_int
  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

REAL(DP) PURE FUNCTION INPUT_DP( That)
  implicit none
  character(*), intent(in) :: That
  integer(i2b) :: i, j
  j=len(That)
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) input_dp
  end do
END FUNCTION INPUT_DP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

INTEGER(I2B) PURE FUNCTION INPUT_INT( That)
  implicit none
  character(*), intent(in) :: That
  integer(i2b) :: i, j
  j=len(That)
  do concurrent( i = 1: size( input_line) )
    if( input_line( i)( 1:j) == That ) read ( input_line (i) (j+4:j+50) , * ) input_int
  end do
END FUNCTION INPUT_INT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE INPUT
