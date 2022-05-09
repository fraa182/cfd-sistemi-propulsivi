subroutine init
    use variabili
    implicit none
    integer::n,m
    real::Mach,ttotin,ptotin,q


    ! **********************************************************************
    ! Define initial Mach, inlet total temp, inlet total press

    ptotin=1
    ttotin=1
    Mach=0.6

    ! **********************************************************************
    ! Compute temperature, pressure, density in each cell: T(n,m),p(n,m),u1(n,m) with 1<=n<=nc, 1<=m<=mc

    do n=1,nc
    
        do m=1,mc
        
            T(n,m)=ttotin/(1+0.5*(gamma-1)*Mach**2)
            p(n,m)=T(n,m)**(gamma/(gamma-1))
            u1(n,m)=p(n,m)/T(n,m)
            
        end do
        
    end do

    ! **********************************************************************
    ! Compute velocity magnitude (use the scalar temporary variable q, the same in all cells)

    q=sqrt(R*ttotin/(1+0.5*(gamma-1)*Mach**2))

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
