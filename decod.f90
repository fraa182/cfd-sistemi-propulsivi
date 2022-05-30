subroutine decod
    use variabili
    implicit none
    integer::n,m

    ! Perform the double loop on all the cells and update the primitive variables (p,T,S,a,u,v) starting from the conservative variables (U1,U2,U3,U4)

	do n=0,nc
		
		do m=0,mc
		
			T(n,m)=(u4(n,m)/u1(n,m))*(gamma-1)/R
			p(n,m)=u1(n,m)*R*T(n,m)
			S(n,m)=log((T(n,m)**gamma)/(p(n,m)**(gamma-1)))
			a(n,m)=sqrt(gamma*R*T(n,m))
			u(n,m)=u2(n,m)/u1(n,m)
			v(n,m)=u3(n,m)/u1(n,m)
			
		!write(*,*) 'T	p	S	a 	u	v'
		!write(*,*) T(n,m), p(n,m), S(n,m), a(n,m), u(n,m), v(n,m)
		!read(*,*)
			
		end do
		
	end do

end subroutine
