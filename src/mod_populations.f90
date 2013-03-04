MODULE POPULATIONS
  use precision_kinds
  use constants, only: x, y, z
  implicit none
  private
  public :: calc_n, calc_n_momprop

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_N
  use system, only: jx, jy, jz, a0, a1, c, NbVel, fluid, rho, f_ext, inside, n, solid
  implicit none
  integer(kind=i2b) :: l

  ! apply force on all fluid nodes and update populations
  do concurrent( l= lbound(n,4):ubound(n,4) )
    where( inside==fluid)
      n(:,:,:,l) = a0(l)*rho + a1(l)*( c(x,l)*(jx + f_ext(x)) + c(y,l)*(jy + f_ext(y)) + c(z,l)*(jz + f_ext(z)) )
    elsewhere
      n(:,:,:,l) = a0(l)*rho + a1(l)*( c(x,l)*jx + c(y,l)*jy + c(z,l)*jz )
    end where
  end do

  ! check that no population n < 0
  call check_population(n)

END SUBROUTINE CALC_N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CALC_N_MOMPROP
  use system, only: n, f_ext, rho, inside, fluid,&
                     jx, jy, jz, a0, a1, elec_slope_x, elec_slope_y, elec_slope_z,&
                     c, NbVel
  use input, only: input_dp
  implicit none
  integer(kind=i2b) :: l
  real(dp) :: D_tracer, tracer_z
  D_tracer = input_dp('tracer_Db'); if(D_tracer<=0.0_dp) stop 'D_tracer, ie tracer_Db in input is invalid'
  tracer_z = input_dp('tracer_z')

  ! apply force on all fluid nodes and update populations
  do concurrent( l=1:NbVel )
    where( inside==fluid )
      n(:,:,:,l) = a0(l)*rho(:,:,:) &
                 + a1(l)*( c(x,l)*(jx + f_ext(x) - rho*tracer_z *D_tracer *elec_slope_x) &
                         + c(y,l)*(jy + f_ext(y) - rho*tracer_z *D_tracer *elec_slope_y) &
                         + c(z,l)*(jz + f_ext(z) - rho*tracer_z *D_tracer *elec_slope_z) )
    elsewhere
      n(:,:,:,l) = a0(l)*rho + a1(l)*( c(x,l)*jx + c(y,l)*jy + c(z,l)*jz )
    end where
  end do

  ! check that no population n < 0
  call check_population(n)
END SUBROUTINE CALC_N_MOMPROP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHECK_POPULATION (arrayin)
  implicit none
  real(dp), dimension(:,:,:,:), intent(in) :: arrayin
  real(dp), parameter :: permitted_minimum_value_numerical_noise = epsilon(1.0_dp)
  if (any(arrayin > 1.0_dp + permitted_minimum_value_numerical_noise)) stop 'population n > 1 somewhere'
  if (any(arrayin < permitted_minimum_value_numerical_noise)) stop 'some population n < 0 somewhere.'
END SUBROUTINE CHECK_POPULATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE POPULATIONS
