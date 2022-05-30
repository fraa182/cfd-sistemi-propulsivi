
subroutine update_cons
    use variabili
    implicit none
    integer::n,m


    !Perform a double loop on all cells (1<=n<=nc, 1<=m<=mc) and update the conservative variables using the computed fluxes

	do n=1,nc
	
		do m=1,mc
			
			!write(*,*)'Flussi right'
			!write(*,*) n, m, F1right(n-1,m), F2right(n-1,m), F3right(n-1,m), F4right(n-1,m)
			!write(*,*)'Flussi up'
			!write(*,*) n, m, F1up(n,m-1), F2up(n,m-1), F3up(n,m-1), F4up(n,m-1)
			!write(*,*)'Conservative'
			!write(*,*) n, m, u1(n,m), u2(n,m), u3(n,m), u4(n,m)
			!read(*,*)
			!write(*,*)'var termodinamiche'
			!write(*,*) n, m, a(n,m), u(n,m), v(n,m), S(n,m)
			
			u1(n,m)=u1(n,m)-dt/area(n,m)*((F1right(n,m)+F1up(n,m))-(F1right(n-1,m)+F1up(n,m-1)))
			u2(n,m)=u2(n,m)-dt/area(n,m)*((F2right(n,m)+F2up(n,m))-(F2right(n-1,m)+F2up(n,m-1)))
			u3(n,m)=0!u3(n,m)-dt/area(n,m)*((F3right(n,m)+F3up(n,m))-(F3right(n-1,m)+F3up(n,m-1)))
			u4(n,m)=u4(n,m)-dt/area(n,m)*((F4right(n,m)+F4up(n,m))-(F4right(n-1,m)+F4up(n,m-1)))
			
			!write(*,*)'Update conservative'
			!write(*,*) n, m, u1(n,m), u2(n,m), u3(n,m), u4(n,m)
			!write(*,*) x(n,m), y(n,m), area(n,m)
			!read(*,*)
			
		end do
		
	!write(*,*)n
	!read(*,*)
		
	end do




end subroutine
