subroutine init
    use variabili
    implicit none
    integer::n,m
    real::Mach,ttotin,ptotin,q


    ! **********************************************************************
    ! Define initial Mach, inlet total temp, inlet total press

    ptotin=1
    ttotin=1
    Mach=0
	ami=Mach

    ! **********************************************************************
    ! Compute temperature, pressure, density in each cell: T(n,m),p(n,m),u1(n,m) with 1<=n<=nc, 1<=m<=mc

    do n=1,nc
    
        do m=1,mc
        
            T(n,m)=ttotin
            p(n,m)=ptotin
            u1(n,m)=p(n,m)/T(n,m)
            a(n,m)=sqrt(gamma*T(n,m))
            S(n,m)=log((T(n,m)**gamma)/(p(n,m)**(gamma-1)))
            
        end do
        
    end do

    ! **********************************************************************
    ! Compute velocity magnitude (use the scalar temporary variable q, the same in all cells)

    q=Mach*sqrt(gamma) !q=M*a=M*sqrt(gamma*T) -> T=1 perché T=T/Tref, ma T=ttotin (assegno valore iniziale) e Tref=ttotin, quindi T=1

    ! **********************************************************************
    ! Compute cartesian velocity components u(n,m) and v(n,m) in each cell starting from q and normal components nx_up,ny_up with 1<=n<=nc, 1<=m<=mc

    do n=1,nc
    
        do m=1,mc
        
            u(n,m)=q/sqrt(1+(nx_up(n,m)/ny_up(n,m))**2)
            v(n,m)=-u(n,m)*nx_up(n,m)/ny_up(n,m)
            
        end do
        
    end do

    ! **********************************************************************
    ! Compute conservative variables u1(n,m), u2(n,m), u3 (n,m), u4(n,m) with 1<=n<=nc, 1<=m<=mc

    do n=1,nc
    
        do m=1,mc
        
            u2(n,m)=u1(n,m)*u(n,m)
            u3(n,m)=u1(n,m)*v(n,m)
            u4(n,m)=u1(n,m)*(1/(gamma-1)*T(n,m)+0.5*q**2)
            
        end do
        
    end do

end subroutine
