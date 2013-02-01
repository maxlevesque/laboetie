subroutine print_everything_related_to_charge_equil

  use system, only: lx, ly, lz, phi

  implicit none

  ! print internal potential
  call print4darray(lx,ly,lz,phi,'output/phi_of_x_y_z.dat')

end subroutine print_everything_related_to_charge_equil
