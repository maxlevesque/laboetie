subroutine backup_phi_c_plus_c_minus

  use system, only: phi, phi_old, c_plus, c_plus_old, c_minus, c_minus_old, lx, ly, lz
  implicit none

  ! the original distribution of concentrations should already be initiated before backing it up
  if( .not. allocated( c_plus )) then
    stop 'Solute + concentrations should already been initiated before backing it up. stop in backup_phi_c_plus_c_minus.f90'
  end if
  if( .not. allocated( c_minus )) then
    stop 'Solute - concentrations should already been initiated before backing it up. stop in backup_phi_c_plus_c_minus.f90'
  end if

  ! the electrostatic potential should already be initiated
  if( .not. allocated( phi )) then
    stop 'The electrostatic potential should already be initiated in backup_phi_c_plus_c_minus.f90. stop.'
  end if

  ! if it is the first iteration in find_equilibrium_charge_distribution then allocate backup arrays
  if( .not. allocated( c_plus_old) ) allocate( c_plus_old (lx,ly,lz) )
  if( .not. allocated( c_minus_old) ) allocate( c_minus_old (lx,ly,lz) )
  if( .not. allocated( phi_old) ) allocate ( phi_old (lx,ly,lz) )

  ! backup
  c_plus_old = c_plus
  c_minus_old = c_minus
  phi_old = phi

end subroutine backup_phi_c_plus_c_minus
