subroutine CFL
    use variabili
    implicit none
    integer::n,m
    real::lambda_max,dt_loc,stab,dx
    
    stab=0.8

    ! Perform a double loop on all the cells ( 1<=n<=nc, 1<=m<=mc)

    dt=100000. !Initial tentative value for the global time step: if you find something more restrictive, correct it

    do m=1,mc
    
        do n=1,nc
            ! In each cell compute the max propagation speed
			
			lambda_max=max(sqrt(u(n,m)**2+v(n,m)**2)+a(n,m),sqrt(u(n,m)**2+v(n,m)**2)+a(n,m))
			
            ! In each cell compute the max allowable time step
            
            dx=x(n,m)-x(n-1,m)
            dt_loc=dx/lambda_max

            ! If local time step is less restrictive then global time step then update the global time step
            
            if (dt_loc<dt) then
            
				dt=dt_loc
				
			end if

        end do
        
    end do
    
    dt=dt*stab

end subroutine
