! precision_kinds is the module defining the precision of each variables.
! This allows architecture independant programming.
! defining the precision kind this way we say to the computer : use the default single precision and double precision.
! remember that this default value may change depending on the architecture of your CPU.
! but one might one day wish to change that. It would be such a mess do change every number in the program !
! One would just have to change here the definition of single and double precision.
! For a very nice lesson about this : http://www.owlnet.rice.edu/~ceng303/manuals/fortran/FOR2_5.html

module precision_kinds
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, parameter, public :: i0b = INT8  !  8 bits integer, from -128 to 127, ie 2^8 -1 values
    integer, parameter, public :: i1b = INT16 ! 16 bits integer, from -32,768 to 32,767
    integer, parameter, public :: i2b = INT32 ! 32 bits integer, from -2,147,483,648 to 2,147,483,647
    integer, parameter, public :: i4b = INT64 ! 64 bits integer, from -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
    integer, parameter, public :: sp = REAL32 ! 32 bits real number (often called single precision)
    integer, parameter, public :: dp = REAL64 ! 64 bits real number (often called double precision)
    integer, parameter, public :: qp = REAL128 ! 128 bits real number (often called quad precision)
end module precision_kinds
