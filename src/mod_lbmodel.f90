!... Where we define lattices using DnQm classification

module mod_lbmodel

    use precision_kinds
    use constants, only: x, z

    implicit none

    type type_velocity
        integer(i2b), dimension(x:z) :: coo
        real(dp) :: a0, a1, delta
        integer(i2b) :: inv
    end type
    type type_lbmodel
        character(len=50) :: name
        integer(i2b) :: dimension ! e.g. 3 in D3Q19
        integer(i2b) :: nvel ! number of velocities, e.g. 19 in D3Q19
        integer(i2b) :: lmin, lmax ! velmin can be 0 or 1, so is velmax 18 or 19 for instance for D3Q19
        type (type_velocity), dimension(:), allocatable :: vel
    end type
    type(type_lbmodel), public :: lbm
    contains

!===================================================================================================================================

        subroutine initialize
            use input
            implicit none
            lbm%name = input_char("lbmodel")
            if (lbm%name(1:1)/="D") then
                stop "critical error. LB model name should begin by D"
            else if (lbm%name(2:2)/="3") then
                stop "critical error. LB model name should have dimension 3"
            else if (lbm%name(3:3)/="Q") then
                stop "critical error. LB model should be of the form D3Q.."
            end if
            read(lbm%name(2:2),'(I4)') lbm%dimension
            read(lbm%name(4:10),'(I4)') lbm%nvel
            lbm%lmin = 1 ! fortran style
            lbm%lmax = lbm%lmin + lbm%nvel -1
            call init_velocities
            call init_weight_factors
            call determine_velocity_inverse
        end subroutine

!===================================================================================================================================

        subroutine init_velocities
            implicit none
            allocate (lbm%vel(lbm%lmin:lbm%lmax))
            select case (lbm%nvel)
            case (15) ! D3Q15
                lbm%vel(lbm%lmin+0)%coo = [0,0,0]
                lbm%vel(lbm%lmin+1)%coo = [1,0,0]
                lbm%vel(lbm%lmin+2)%coo = [-1,0,0]
                lbm%vel(lbm%lmin+3)%coo = [0,1,0]
                lbm%vel(lbm%lmin+4)%coo = [0,-1,0]
                lbm%vel(lbm%lmin+5)%coo = [0,0,1]
                lbm%vel(lbm%lmin+6)%coo = [0,0,-1]
                lbm%vel(lbm%lmin+7)%coo = [1,1,1]
                lbm%vel(lbm%lmin+8)%coo = [-1,1,1]
                lbm%vel(lbm%lmin+9)%coo = [1,-1,1]
                lbm%vel(lbm%lmin+10)%coo = [1,1,-1]
                lbm%vel(lbm%lmin+11)%coo = [-1,-1,1]
                lbm%vel(lbm%lmin+12)%coo = [-1,1,-1]
                lbm%vel(lbm%lmin+13)%coo = [1,-1,-1]
                lbm%vel(lbm%lmin+14)%coo = [-1,-1,-1]
            case (19) ! D3Q19
                lbm%vel(lbm%lmin+0)%coo = [0,0,0]
                lbm%vel(lbm%lmin+1)%coo = [1,0,0]
                lbm%vel(lbm%lmin+2)%coo = [-1,0,0]
                lbm%vel(lbm%lmin+3)%coo = [0,1,0]
                lbm%vel(lbm%lmin+4)%coo = [0,-1,0]
                lbm%vel(lbm%lmin+5)%coo = [0,0,1]
                lbm%vel(lbm%lmin+6)%coo = [0,0,-1]
                lbm%vel(lbm%lmin+7)%coo = [1,1,0]
                lbm%vel(lbm%lmin+8)%coo = [-1,1,0]
                lbm%vel(lbm%lmin+9)%coo = [1,-1,0]
                lbm%vel(lbm%lmin+10)%coo = [-1,-1,0]
                lbm%vel(lbm%lmin+11)%coo = [1,0,1]
                lbm%vel(lbm%lmin+12)%coo = [-1,0,1]
                lbm%vel(lbm%lmin+13)%coo = [1,0,-1]
                lbm%vel(lbm%lmin+14)%coo = [-1,0,-1]
                lbm%vel(lbm%lmin+15)%coo = [0,1,1]
                lbm%vel(lbm%lmin+16)%coo = [0,-1,1]
                lbm%vel(lbm%lmin+17)%coo = [0,1,-1]
                lbm%vel(lbm%lmin+18)%coo = [0,-1,-1]
            case (27) ! D3Q27
                lbm%vel(lbm%lmin+0)%coo = [0,0,0]
                lbm%vel(lbm%lmin+1)%coo = [1,0,0]
                lbm%vel(lbm%lmin+2)%coo = [-1,0,0]
                lbm%vel(lbm%lmin+3)%coo = [0,1,0]
                lbm%vel(lbm%lmin+4)%coo = [0,-1,0]
                lbm%vel(lbm%lmin+5)%coo = [0,0,1]
                lbm%vel(lbm%lmin+6)%coo = [0,0,-1]
                lbm%vel(lbm%lmin+7)%coo = [1,1,0]
                lbm%vel(lbm%lmin+8)%coo = [-1,1,0]
                lbm%vel(lbm%lmin+9)%coo = [1,-1,0]
                lbm%vel(lbm%lmin+10)%coo = [-1,-1,0]
                lbm%vel(lbm%lmin+11)%coo = [1,0,1]
                lbm%vel(lbm%lmin+12)%coo = [-1,0,1]
                lbm%vel(lbm%lmin+13)%coo = [1,0,-1]
                lbm%vel(lbm%lmin+14)%coo = [-1,0,-1]
                lbm%vel(lbm%lmin+15)%coo = [0,1,1]
                lbm%vel(lbm%lmin+16)%coo = [0,-1,1]
                lbm%vel(lbm%lmin+17)%coo = [0,1,-1]
                lbm%vel(lbm%lmin+18)%coo = [0,-1,-1]
                lbm%vel(lbm%lmin+19)%coo = [1,1,1]
                lbm%vel(lbm%lmin+20)%coo = [-1,1,1]
                lbm%vel(lbm%lmin+21)%coo = [1,-1,1]
                lbm%vel(lbm%lmin+22)%coo = [1,1,-1]
                lbm%vel(lbm%lmin+23)%coo = [-1,-1,1]
                lbm%vel(lbm%lmin+24)%coo = [-1,1,-1]
                lbm%vel(lbm%lmin+25)%coo = [1,-1,-1]
                lbm%vel(lbm%lmin+26)%coo = [-1,-1,-1]
            case default
                stop "You ask for a DnQm lattice that is not implemented"
            end select
        end subroutine

!===================================================================================================================================

        subroutine init_weight_factors
            implicit none
            real(dp), parameter :: a_00 = 1.0_dp/3.0_dp ! Weighting factors wi for two different LB models for the discrete velocity vectors e_i
            real(dp), parameter :: a_01 = 1.0_dp/18.0_dp
            real(dp), parameter :: a_02 = 1.0_dp/36.0_dp
            real(dp), parameter :: a_10 = 0.0_dp
            real(dp), parameter :: a_11 = 1.0_dp/6.0_dp
            real(dp), parameter :: a_12 = 1.0_dp/12.0_dp
            real(dp), parameter :: itself=0.0_dp, nn=1.0_dp, nnn=sqrt(2.0_dp) ! nearest and next nearest neighbour distance
            if (lbm%nvel==15) then
                stop "lacks weight factors for D3Q15"
            else if (lbm%nvel==19) then
                lbm%vel%a0 =&
  [ a_00, a_01, a_01, a_01, a_01, a_01, a_01, a_02, a_02, a_02, a_02, a_02, a_02, a_02, a_02, a_02, a_02, a_02, a_02 ]
                lbm%vel%a1 =&
  [ a_10, a_11, a_11, a_11, a_11, a_11, a_11, a_12, a_12, a_12, a_12, a_12, a_12, a_12, a_12, a_12, a_12, a_12, a_12 ]
                lbm%vel%delta = [ itself, nn, nn, nn, nn, nn, nn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn ]
            else if (lbm%nvel==27) then
                stop "lacks weight factors for D3Q27"
            end if
            if ( abs(sum(lbm%vel%a0))-1._dp > epsilon(1.0)) stop "The sum of the weights of the velocities (a0) must be 1."
        end subroutine

!===================================================================================================================================

        subroutine determine_velocity_inverse
            implicit none
            integer(i2b) :: l, li
            do concurrent (l=lbm%lmin:lbm%lmax, li=lbm%lmin:lbm%lmax)
                if (all (lbm%vel(l)%coo == -lbm%vel(li)%coo) ) lbm%vel(l)%inv = li
            end do
        end subroutine

!===================================================================================================================================

end module
