subroutine decod
    use variabili
    implicit none
    integer::n,m

    ! Perform the double loop on all the cells and update the primitive variables (p,T,S,a,u,v) starting from the conservative variables (U1,U2,U3,U4)

	do n=1,nc
		
		do m=1,mc
		
		    u(n,m)=u2(n,m)/u1(n,m)
			v(n,m)=u3(n,m)/u1(n,m)
			T(n,m)=(u4(n,m)/u1(n,m)-0.5*(u(n,m)**2+v(n,m)**2))*(gamma-1)
			p(n,m)=u1(n,m)*T(n,m)
			S(n,m)=log((T(n,m))/(p(n,m)**((gamma-1)/gamma)))
			a(n,m)=sqrt(gamma*T(n,m))
			
		!write(*,*) 'n  m   T	p	S	a 	u	v'
		!write(*,*) n, m, T(n,m), p(n,m), S(n,m), a(n,m), u(n,m), v(n,m)
		!read(*,*)
			
		end do
		
	end do

end subroutine
