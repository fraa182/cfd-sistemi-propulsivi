module variabili
    save
    integer,parameter::nc=30,mc=10
    real,dimension(0:nc,0:mc)::x,y,nx_up,ny_up,nx_right,ny_right,length_right,length_up
    real,dimension(0:nc,0:mc)::F1right,F2right,F3right,F4right,F1up,F2up,F3up,F4up
    real,dimension(nc,mc)::area,xg,yg
    real,dimension(nc,mc)::u1,u2,u3,u4,p,T,u,v,a,S

    integer::kend,kinf,kout

    real::time,dt
    real,parameter::pi=DACOS(-1.D0)

    ! Geometry
    real::x1,x2,x3,x4,h
    real,dimension(0:nc)::ybottom,ytop

end module