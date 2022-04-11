subroutine grid
    use variabili
    implicit none
    integer::n,m
    real::b1, dx, dy

!     Element numbering (1<=n<=nc, 1<=m<=mc)
!            * ----- *
!            |       |
!            |(n,m+1)|
!            |       |
!    * ----- * ----- * ----- *
!    |       |       |       |
!    |(n-1,m)| (n,m) |(n+1,m)|
!    |       |       |       |
!    *-------*-------*-------*
!            |       |
!            |(n,m-1)|
!            |       |
!            * ----- *

    ! *******************************************************************

!     Node numbering (0<=n<=nc, 0<=m<=mc) -> partono da 0 perchè sempre un nodo in più 


!        *Node(n-1,m)-------*Node(n,m)
!        |                  |
!        |                  |
!        |                  |
!        |                  |
!        |     ele(n,m)     |
!        |                  |
!        |                  |
!        |                  |
!        |                  |
!        *Node(n-1,m-1)-----*Node(n,m-1)





    !********************************************************************
    ! Define geometry (x1,x2,h,b(:), c(:)
    ! Compute node coordinates x(n,m), y(n,m)  with 0<= n <= nc, 0<= m <= mc -> matrici

    x1=-1;
    x2=2;
    h=1;
    b1=0.5;
    
    ytop(:)=h;
    
    do n=0,nc
    
        do m=0,mc
        
            x(n,m)=x1+(x2-x1)*(1.*n/nc)
            if (x(n,m)<0 .or. x(n,m)>1) then
            
                ybottom(n)=0;
                
            else
            
                ybottom(n)=b1*sin(pi*x(n,m));
                
            end if
            y(n,m)=ybottom(n)+(ytop(n)-ybottom(n))*(1.*m/mc)
            
        end do
    end do

    !********************************************************************
    ! Compute edge length and normal on right edge nx_right(n,m), ny_right(n,m), length_right(n,m)  with 0<= n <= nc, 1<= m <= mc

    do n=0,nc
    
        do m=1,mc
        
            dx=x(n,m)-x(n,m-1)
            dy=y(n,m)-y(n,m-1)
            
            length_right(n,m)=sqrt(dx**2+dy**2)
            
            nx_right(n,m)=dy/length_right(n,m)
            ny_right(n,m)=-dx/length_right(n,m)
        
        end do
        
    end do

    !********************************************************************
    ! Compute edge length and normal on upper edge nx_up(n,m), ny_up(n,m), length_up(n,m)  with 1<= n <= nc, 0<= m <= mc

    do n=1,nc
    
        do m=0,mc
        
            dx=x(n,m)-x(n-1,m)
            dy=y(n,m)-y(n-1,m)
            
            length_up(n,m)=sqrt(dx**2+dy**2)
            
			nx_up(n,m)=-dy/length_up(n,m)
            ny_up(n,m)=dx/length_up(n,m)
        
        end do
        
    end do

    !********************************************************************
    ! Compute element area(n,m) and center coordinates xg(n,m), yg(n,m) with 1<= n <= nc, 1<= m <= mc

    do n=1,nc
    
        do m=1,mc
        
            area(n,m)=0.5*abs(x(n-1,m-1)*y(n,m-1)-x(n,m-1)*y(n-1,m-1)+x(n,m-1)*y(n,m)-x(n,m)*y(n,m-1)+x(n,m)*y(n-1,m)&
            -x(n-1,m)*y(n,m)+x(n-1,m)*y(n-1,m-1)-x(n-1,m-1)*y(n-1,m))
        
            xg(n,m)=(x(n-1,m-1)+x(n,m-1)+x(n,m)+x(n-1,m))/4
            yg(n,m)=(y(n-1,m-1)+y(n,m-1)+y(n,m)+y(n-1,m))/4
        
        end do
        
    end do

    write(*,*)'nx_up:'
    write(*,*)' '
    do m=1,mc
    
        write(*,*)nx_up(:,m)
        
    end do
    write(*,*)'ny_up:'
    write(*,*)' '
    do m=0,mc
    
        write(*,*)ny_up(:,m)
        
    end do
    write(*,*)' '
    write(*,*)'nx_right:'
    write(*,*)' '
    do m=1,mc
    
        write(*,*)nx_right(:,m)
        
    end do
    write(*,*)'ny_right:'
    write(*,*)' '
    do m=0,mc
    
        write(*,*)ny_right(:,m)
        
    end do
    write(*,*)' '
    write(*,*)'xg:'
	write(*,*)' '
    do m=1,mc
    
		write(*,*)xg(:,m)
		
	end do
	write(*,*)'yg:'
	write(*,*)' '
	do m=1,mc
    
		write(*,*)yg(:,m)
		
	end do


end subroutine
