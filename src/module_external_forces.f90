module module_external_forces
    implicit none
    private
    public apply_external_forces
contains
    
    subroutine apply_external_forces( fextx, fexty, fextz, nature)
        use precision_kinds, only: dp
        use system, only: fluid
        use module_input, only: getinput
        implicit none
        real(dp), intent(out), dimension(:,:,:) :: fextx, fexty, fextz
        integer(kind(fluid)), intent(in) :: nature(:,:,:)
        real(dp) :: fext_tmp(3) ! external forces as requested in input file
        integer :: particleDiameter, particleRadius, lx, ly, lz, countFluidNodes, particleCoordinates(3), px, py, pz, countNodesInParticle
        integer :: i, j, k, geometryLabel
        logical :: compensate_f_ext
        real(dp), parameter :: zerodp = 0._dp
        
        lx = ubound(nature,1)
        ly = ubound(nature,2)
        lz = ubound(nature,3)
        fext_tmp = getinput%dp3("f_ext", defaultvalue= [0._dp,0._dp,0._dp] )
        compensate_f_ext = getinput%log( "compensate_f_ext", defaultvalue = .false.)

        ! IN THE CASE I DONT WANT A NEUTRALIZING / COMPENSATING BACKGROUND
        if(.not. compensate_f_ext) then ! the force is applied everywhere with same intensity, i.e., homogeneously
            ! We also apply the extrnal forces on the solid nodes. Nobody cares about constraints on solid nodes since density is zero there.
            where( nature == fluid )
                fextx = fext_tmp(1)
                fexty = fext_tmp(2)
                fextz = fext_tmp(3)
            end where
            ! It is overkill to have a whole array of the same value.
        ! IF I WANT A COMPENSATING BACKGROUND. We call "particle" the ensemble of nodes on which the singular force is applied. All other nodes see a compensating background.
        else if(compensate_f_ext) then ! force applied to a central particle only
            ! start by initiating external force on all nodes to zero
            fextx = zerodp
            fexty = zerodp
            fextz = zerodp
            ! then get informations about the particle
            particleDiameter = getinput%int("dominika_particle_diameter", defaultvalue=1) ! the particle diameter
            if( modulo(particleDiameter, 2) == 0 ) error stop "Dominika's particle have an even diameter. It must be odd."
            if( modulo(lx,2)==0 .or. modulo(ly,2)==0 .or. modulo(lz,2)==0) then
                error stop "With compensate_f_ext, there should be odd number of nodes in all directions"
            end if
            particleRadius = particleDiameter/2 ! nodes of the particle on the right (or left) of the particle center. If particle is diameter 3, we have 1 node on the left and 1 on the right, so pd=3, pdr=3/2=1
            particleCoordinates = getinput%int3("particle_coordinates", defaultvalue=[lx/2+1,ly/2+1,lz/2+1] )
            px = particleCoordinates(1)
            py = particleCoordinates(2)
            pz = particleCoordinates(3)
            open(47, file = "./output/dominika_particle_shape.xyz")
            countNodesInParticle=0 ! ADE: l counts the number of node within the particle
            do k = pz-particleRadius, pz+particleRadius
                do j = py-particleRadius, py+particleRadius
                    do i = px-particleRadius, px+particleRadius
                        if( (i-px)**2+(j-py)**2+(k-pz)**2 > particleRadius**2 ) then
                            cycle ! node (i,j,k) is outside the particle 
                        else
                            if ( nature(i,j,k) /= fluid ) then
                                error stop "The so-called Dominika's particle contains solid nodes"
                            end if
                            countNodesInParticle = countNodesInParticle + 1 ! One more node within the particle
                            fextx(i,j,k) = fext_tmp(1) ! apply the user requestd force to this node that is inside the particle
                            fexty(i,j,k) = fext_tmp(2)
                            fextz(i,j,k) = fext_tmp(3)
                            write(47,*) i, j, k ! use ListPointPlot3D[data,BoxRatios->{1,1,1}] in Mathematica to read this file
                        end if
                    end do
                end do
            end do
            close(47)
            ! ADE: We distribute the total force upon the particle evenly
            ! throughout the various nodes
            geometryLabel = getinput%int("geometryLabel", defaultvalue=0) ! if geometryLabel=-1 =>bulk case
            ! ADE : I modified the following lines
            ! the idea is that whenever we have a slit kind of geometry,
            ! we do not want to add a compensation force in the rest of
            ! the nodes, as we believe that the force will dissipate
            ! within the walls
            if (geometryLabel==-1) then
                countFluidNodes = count( nature == fluid )
                where(fextx==fext_tmp(1) .and. fexty==fext_tmp(2).and.fextz==fext_tmp(3) )
                    fextx = -fext_tmp(1)/(countFluidNodes) + fext_tmp(1) / countNodesInParticle
                    fexty = -fext_tmp(2)/(countFluidNodes) + fext_tmp(2) / countNodesInParticle
                    fextz = -fext_tmp(3)/(countFluidNodes) + fext_tmp(3) / countNodesInParticle
                elsewhere
                    fextx = -fext_tmp(1)/(countFluidNodes)
                    fexty = -fext_tmp(2)/(countFluidNodes)
                    fextz = -fext_tmp(3)/(countFluidNodes)
                end where
                if( any(abs([sum(fextx)/countFluidNodes,sum(fexty)/countFluidNodes,sum(fextz)/countFluidNodes])> epsilon(1._dp) ) ) then
                    print*,"=====  The compensation is not well-implemented."
                    print*,"       sum(fextx)=",sum(fextx)
                    print*,"       sum(fexty)=",sum(fexty)
                    print*,"       sum(fextz)=",sum(fextz)
                    stop
                end if
            else
                where(fextx==fext_tmp(1) .and. fexty==fext_tmp(2) .and.fextz==fext_tmp(3) )
                    fextx = fextx / countNodesInParticle
                    fexty = fexty / countNodesInParticle
                    fextz = fextz / countNodesInParticle
                elsewhere
                    fextx = zerodp
                    fexty = zerodp
                    fextz = zerodp
                end where
            endif
            where(nature/=fluid) ! Here also PAY ATTENTION TO WHAT YOU MEAN. NATURE==FLUID, INTERFACIAL ETC
                fextx = zerodp
                fexty = zerodp
                fextz = zerodp
            end where
        end if ! compensate
    end subroutine apply_external_forces
end module module_external_forces