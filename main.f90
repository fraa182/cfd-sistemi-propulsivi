program CFD
    use variabili
    implicit none
    integer::k
    real::pexit

    ! ***********************************************************
    ! INPUT SECTION
    kend=10   ! Total number of time steps
    kout=100    ! Save solution every kout time steps
    kinf=1    ! Write on screen every kinf time steps

    ! ***********************************************************
    
	call system("del *.exe")
    call system("del *.plt")
    
    write(*,*)'nc,mc = ',nc,mc

    ! Generate grid
    call grid

    ! Save grid on file
    call WDKs_tk(0,-3, "grid", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    call init
    
    ! Save initial condition on file
    call WDKs_tk(0,1, "sol", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    time=0.
    pexit=1

    do k=1,kend

        call CFL
		call flux_right_tilde
		call flux_up_tilde
		call bordi_s1(x,y,nx_right,ny_right,nx_up,ny_up,a,u,v,s,F1right,F2right,F3right,F4right,F1up,F2up,&
	F3up,F4up,nc,mc,nc,mc,ami,pexit)
        call update_cons
		
		
	!write(*,*) 'u1 = ',u1
		
       !call decod

        if(mod(k,kinf).eq.0) write(*,*)'k = ',k
       !if(mod(k,kout).eq.0) call WDKs_tk(k,1, "sol", x,y,nc,mc,nc,mc,xg,yg,a,u,v,s)

    end do

end program