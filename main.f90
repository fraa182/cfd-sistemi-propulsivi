program CFD
    use variabili
    implicit none
    integer::k


    ! ***********************************************************
    ! INPUT SECTION
    kend=1000   ! Total number of time steps
    kout=100    ! Save solution every kout time steps
    kinf=10     ! Write on screen every kinf time steps



    ! ***********************************************************


    write(*,*)'nc,mc = ',nc,mc

    ! Generate grid
    call grid

    ! Save grid on file
    !call WDKs_tk(0,-3, "grid", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    call init
    ! Save initial condition on file
    !call WDKs_tk(0,1, "sol", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    time=0.

    do k=1,kend

        call CFL
        call flux_right_tilde
        Call flux_up_tilde
        !call bordi
        !call update_cons
        !call decod

        !if(mod(k,kinf).eq.0) write(*,*)'k = ',k
        !if(mod(k,kout).eq.0) call WDKs_tk(k,1, "sol", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    end do

end program

