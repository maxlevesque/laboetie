!... Where we define lattices using DnQm classification

module mod_lbmodel

    use precision_kinds, only: dp
    implicit none

    private
    type type_velocity
        integer, dimension(3) :: coo
        real(dp) :: a0, a1, a2, delta
        integer :: inv
    end type
    type type_lbmodel
        character(len=50) :: name="D3Q19"
        integer :: dim=3 ! e.g. 3 in D3Q19
        integer :: nvel=19 ! number of velocities, e.g. 19 in D3Q19
        integer :: lmin=1, lmax=19 ! velmin can be 0 or 1, so is velmax 18 or 19 for instance for D3Q19
        integer, allocatable :: cx(:), cy(:), cz(:), linv(:)
        real(dp) :: csq ! velocity of sound, squared, c²
        real(dp), allocatable :: a0(:), a1(:), a2(:), delta(:)
        type (type_velocity), dimension(:), allocatable :: vel
    end type
    type(type_lbmodel), public, protected :: lbm
    public :: init

contains

!===================================================================================================================================

    SUBROUTINE init
        use module_input, only: getinput
        lbm%name = getinput%char("lbmodel", defaultvalue="D3Q19")
        IF( lbm%name(1:3) /= "D3Q" ) ERROR STOP "I dont understand your lbmodel"
        READ( lbm%name(2:2),'(I4)')  lbm%dim
        READ( lbm%name(4:10),'(I4)') lbm%nvel
        lbm%lmin = 1
        lbm%lmax = lbm%lmin + lbm%nvel -1
        call init_velocities
        call init_weight_factors
        call determine_velocity_inverse
    END SUBROUTINE init

!===================================================================================================================================

    SUBROUTINE init_velocities
        implicit none
        allocate (lbm%vel(lbm%lmin:lbm%lmax))
        if(lbm%lmin/=1)error stop "lbm%lmin/=1 in init_velocities in module_lbm"

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
            lbm%csq = 1.0_dp/3.0_dp ! velocity of sound
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
        allocate(  lbm%cx(lbm%lmin:lbm%lmax), source=lbm%vel(lbm%lmin:lbm%lmax)%coo(1)   )
        allocate(  lbm%cy(lbm%lmin:lbm%lmax), source=lbm%vel(lbm%lmin:lbm%lmax)%coo(2)   )
        allocate(  lbm%cz(lbm%lmin:lbm%lmax), source=lbm%vel(lbm%lmin:lbm%lmax)%coo(3)   )
    end subroutine

!============================================================================================================================

        subroutine init_weight_factors
            implicit none
            ! a0_l = w_l, where w_i are the weights associated to each velocity
            ! To simplify calculations in lattice units, we define a0_l, a1_l and a2_l
            ! a1_l = w_l/cs², where cs is the speed of sound and cs²=1/3
            real(dp), parameter :: csq = 1.0_dp/3.0_dp
            real(dp), parameter :: a_00 = 1.0_dp/3.0_dp ! Weighting factors w_i for two different LB models for the discrete velocity vectors e_i
            real(dp), parameter :: a_01 = 1.0_dp/18.0_dp
            real(dp), parameter :: a_02 = 1.0_dp/36.0_dp
            real(dp), parameter :: a_10 = a_00/csq
            real(dp), parameter :: a_11 = a_01/csq
            real(dp), parameter :: a_12 = a_02/csq
            real(dp), parameter :: a_20 = a_00/(2*csq**2)
            real(dp), parameter :: a_21 = a_01/(2*csq**2)
            real(dp), parameter :: a_22 = a_02/(2*csq**2)
            real(dp), parameter :: itself=0.0_dp, nn=1.0_dp, nnn=sqrt(2.0_dp) ! nearest and next nearest neighbour distance
            integer :: l

            select case (lbm%nvel)
            case (15)
                stop "lacks weight factors for D3Q15"

            case (19)

lbm%vel%a0=[a_00,a_01,a_01,a_01,a_01,a_01,a_01,a_02,a_02,a_02,a_02,a_02,a_02,a_02,a_02,a_02,a_02,a_02,a_02]
lbm%vel%a1=[a_10,a_11,a_11,a_11,a_11,a_11,a_11,a_12,a_12,a_12,a_12,a_12,a_12,a_12,a_12,a_12,a_12,a_12,a_12]
lbm%vel%a2=[a_20,a_21,a_21,a_21,a_21,a_21,a_21,a_22,a_22,a_22,a_22,a_22,a_22,a_22,a_22,a_22,a_22,a_22,a_22]

lbm%vel%delta = [ itself, nn, nn, nn, nn, nn, nn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn, nnn ]

allocate(lbm%a0(lbm%lmin:lbm%lmax))
allocate(lbm%a1(lbm%lmin:lbm%lmax))
allocate(lbm%a2(lbm%lmin:lbm%lmax))
allocate(lbm%delta(lbm%lmin:lbm%lmax))
do l=1,19
    lbm%a0(l)=lbm%vel(l)%a0
    lbm%a1(l)=lbm%vel(l)%a1
    lbm%a2(l)=lbm%vel(l)%a2
    lbm%delta(l)=lbm%vel(l)%delta
end do

            case (27)
                stop "lacks weight factors for D3Q27"
            end select

            if ( abs(sum(lbm%vel%a0))-1._dp > epsilon(1.0)) stop "The sum of the weights of the velocities (a0) must be 1."
        end subroutine

!=============================================================================================================================

        subroutine determine_velocity_inverse
            implicit none
            integer :: l, li
            do concurrent (l=lbm%lmin:lbm%lmax, li=lbm%lmin:lbm%lmax)
                if (all (lbm%vel(l)%coo == -lbm%vel(li)%coo) ) lbm%vel(l)%inv = li
            end do
            allocate(lbm%linv(lbm%lmin:lbm%lmax))
            do l=lbm%lmin,lbm%lmax
                lbm%linv(l)=lbm%vel(l)%inv
            end do
        end subroutine

!=============================================================================================================================

end module
