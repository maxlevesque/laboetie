subroutine close_simu
    use precision_kinds
    use system, only: tic, tic_eq_flux, tac_eq_flux, tic_mp, tac_mp, tic_eq_ads, tac_eq_ads
    use OMP_lib
    implicit none
    ! total execution time
    call timeExecution
    contains
        subroutine timeExecution

            real(sp) :: t_total, t_eq_flux, t_eq_ads, t_mp
            integer(i4b) :: tac, count_rate
            character(len=5) :: unite

            call SYSTEM_CLOCK(tac, count_rate)
            t_total = real(tac - tic)/count_rate
            t_eq_flux = real(tac_eq_flux - tic_eq_flux)/count_rate
            t_eq_ads = real(tac_eq_ads - tic_eq_ads)/count_rate
            t_mp = real(tac_mp - tic_mp)/count_rate

            !if (t_total > 3600) then
                !t_total = t_total/3600
               ! t_eq_flux = t_eq_flux/3600
              !  t_eq_ads = t_eq_ads/3600
             !   t_mp = t_mp/3600
            !    unite = 'h'
            !else if (3600 > t_total .and. t_total > 60) then
                t_total = t_total/60
                t_eq_flux = t_eq_flux/60
                t_eq_ads = t_eq_ads/60
                t_mp = t_mp/60
                unite = 'min'
            !else if (60 > t_total) then
            !    unite = 's'
            !end if

            print*, 'Flux Equilibration time (',trim(unite),') = ',t_eq_flux
            print*, 'Adsorption Equilibration time (',trim(unite),') = ',t_eq_ads
            print*, 'Moment propagation time (',trim(unite),') = ',t_mp
            print*,'Total simulation time (',trim(unite),') = ',t_total  
            open(unit=1, file='output/time_report.dat')
            write(1,*) 'Flux Equilibration time (',trim(unite),') = ',t_eq_flux
            write(1,*) 'Adsorption Equilibration time (',trim(unite),') = ',t_eq_ads
            write(1,*) 'Moment propagation time (',trim(unite),') = ',t_mp
            write(1,*) 'Total simulation time (',trim(unite),') = ',t_total    
            close(1)

        end subroutine timeExecution
end subroutine close_simu
