MODULE POPULATIONS
  use precision_kinds
  use constants, only: x, y, z
  use system, only: supercell
  use mod_lbmodel, only: lbm
  implicit none
  private
  public :: calc_n, calc_n_momprop

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_N
  use system, only: jx, jy, jz, fluid, rho, f_ext, n, solid
  integer(kind=i2b) :: l
  ! apply force on all fluid nodes and update populations
  do l= lbm%lmin, lbm%lmax
    where (supercell%node%nature == fluid)
      n(:,:,:,l) = lbm%vel(l)%a0*rho + lbm%vel(l)%a1*( lbm%vel(l)%coo(x)*(jx + f_ext(x)) + &
          lbm%vel(l)%coo(y)*(jy + f_ext(y)) + lbm%vel(l)%coo(z)*(jz + f_ext(z)) )
    elsewhere
      n(:,:,:,l) = lbm%vel(l)%a0*rho + lbm%vel(l)%a1*( lbm%vel(l)%coo(x)*jx + lbm%vel(l)%coo(y)*jy + lbm%vel(l)%coo(z)*jz )
    end where
  end do

  ! check that no population n < 0
  call check_population(n)

END SUBROUTINE CALC_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_N_MOMPROP
  use system, only: n, f_ext, rho, fluid, jx, jy, jz, elec_slope
  use input, only: input_dp
  integer(kind=i2b) :: l
  real(dp) :: D_tracer, tracer_z
  D_tracer = input_dp('tracer_Db')
  if (D_tracer<=0.0_dp) stop 'D_tracer, ie tracer_Db in input is invalid'
  tracer_z = input_dp('tracer_z')

  ! apply force on all fluid nodes and update populations
  do concurrent( l= lbm%lmin: lbm%lmax )
    where( supercell%node%nature ==fluid )
      n(:,:,:,l) = lbm%vel(l)%a0*rho &
                 + lbm%vel(l)%a1*( lbm%vel(l)%coo(x)*(jx + f_ext(x) - rho*tracer_z *D_tracer *elec_slope(x)) &
                         + lbm%vel(l)%coo(y)*(jy + f_ext(y) - rho*tracer_z *D_tracer *elec_slope(y)) &
                         + lbm%vel(l)%coo(z)*(jz + f_ext(z) - rho*tracer_z *D_tracer *elec_slope(z)) )
    elsewhere
      n(:,:,:,l) = lbm%vel(l)%a0*rho + lbm%vel(l)%a1*( lbm%vel(l)%coo(x)*jx + lbm%vel(l)%coo(y)*jy + lbm%vel(l)%coo(z)*jz )
    end where
  end do

  ! check that no population n < 0
  call check_population(n)
END SUBROUTINE CALC_N_MOMPROP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE CHECK_POPULATION (arrayin)
        real(dp), dimension(:,:,:,:), intent(in) :: arrayin
        real(dp), parameter :: permitted_minimum_value_numerical_noise = epsilon(1.0_dp)
        if (any(arrayin > 1.0_dp + permitted_minimum_value_numerical_noise)) stop 'population n > 1 somewhere'
        if (any(arrayin < permitted_minimum_value_numerical_noise)) stop 'some population n < 0 somewhere.'
    END SUBROUTINE CHECK_POPULATION

END MODULE POPULATIONS
