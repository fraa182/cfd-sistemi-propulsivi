
subroutine update_cons
    use variabili
    implicit none
    integer::n,m


    !Perform a double loop on all cells (1<=n<=nc, 1<=m<=mc) and update the conservative variables using the computed fluxes

	do n=1,nc
	
		do m=1,mc
			
			u1(n,m)=u1(n,m)-dt/area(n,m)*((F1right(n,m)+F1up(n,m))-(F1right(n-1,m)+F1up(n,m-1)))
			u2(n,m)=u2(n,m)-dt/area(n,m)*((F2right(n,m)+F2up(n,m))-(F2right(n-1,m)+F2up(n,m-1)))
			u3(n,m)=0!u3(n,m)-dt/area(n,m)*((F3right(n,m)+F3up(n,m))-(F3right(n-1,m)+F3up(n,m-1)))
			u4(n,m)=u4(n,m)-dt/area(n,m)*((F4right(n,m)+F4up(n,m))-(F4right(n-1,m)+F4up(n,m-1)))
			
		end do
		
	end do




end subroutine
