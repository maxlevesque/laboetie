SUBROUTINE backup_phi_c_plus_c_minus

    USE system, ONLY: phi, phi_old, c_plus, c_plus_old, c_minus, c_minus_old
    USE myallocations, ONLY: allocateReal3D

    IMPLICIT NONE

    ! if it is the first iteration in find_equilibrium_charge_distribution then allocate backup arrays
    IF( .NOT. ALLOCATED( c_plus_old ) ) CALL allocateReal3D( c_plus_old ) !allocate( c_plus_old (lx,ly,lz) )
    IF( .NOT. ALLOCATED( c_minus_old) ) CALL allocateReal3D( c_minus_old ) !allocate( c_minus_old (lx,ly,lz) )
    IF( .NOT. ALLOCATED( phi_old    ) ) CALL allocateReal3D( phi_old) !allocate ( phi_old (lx,ly,lz) )

    c_plus_old = c_plus
    c_minus_old = c_minus
    phi_old = phi

END SUBROUTINE backup_phi_c_plus_c_minus
