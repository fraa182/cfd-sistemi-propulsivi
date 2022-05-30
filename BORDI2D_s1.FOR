	subroutine bordi_s1 (xw,yw,anxx,anxy,anyx,anyy
     x                   ,a,u,v,s
     x                   ,f1,f2,f3,f4
     x                   ,g1,g2,g3,g4
     x                   ,kn,km,nc,mc
     x                   ,ami,pexit)
!	include 'fl2d.inc'

        real, dimension (0:kn,0:km)    :: xw,yw,anxx,anxy,anyx,anyy
        real, dimension (0:kn,0:km)    :: f1,f2,f3,f4,g1,g2,g3,g4
        real, dimension (  kn,  km)    :: a,u,v,s
        real, dimension (  *)          :: ami

!	common ami(km),bet0,betc
      dimension Uaf(6),Ubf(6),Ucf(6),Ff(6),Fg(6)
      dimension T0_i(km),s_i(km),sig(km)
       T0_i=1.
	 s_i=0.
	 sig=0.
	 nLE=0
	 nTE=nc
c======== parete m=0===================
	m=0
       call b_wall("b",CMPLX(xw(0:,m) ,yw(0:,m) )
     x                ,anyx(0:,m),anyy(0:,m),nLE+1,nTE         
     x                ,a(:,m+1),u(:,m+1),v(:,m+1),s(:,m+1)            
     x                ,g1(0:,m),g2(0:,m),g3(0:,m),g4(0:,m))

      if(nLE.gt.0) then 
c          ============= periodico ===============

	      call b_peri("b",CMPLX(xw(0:,m) ,yw(0:,m) )
     x                     ,anyx(0:,m),anyy(0:,m),1,nLE        
     x                     ,a(:,m+1),u(:,m+1),v(:,m+1),s(:,m+1)            
     x                     ,a(:,mc),u(:,mc),v(:,mc),s(:,mc)            
     x                     ,g1(0:,m),g2(0:,m),g3(0:,m),g4(0:,m))

	       g1(1:nLE,mc)=g1(1:nLE,m)
	       g2(1:nLE,mc)=g2(1:nLE,m)
	       g3(1:nLE,mc)=g3(1:nLE,m)
	       g4(1:nLE,mc)=g4(1:nLE,m)
	      
		  call b_peri("b",CMPLX(xw(0:,m) ,yw(0:,m) )
     x                     ,anyx(0:,m),anyy(0:,m),nTE+1,nc       
     x                     ,a(:,m+1),u(:,m+1),v(:,m+1),s(:,m+1)            
     x                     ,a(:,mc),u(:,mc),v(:,mc),s(:,mc)            
     x                     ,g1(0:,m),g2(0:,m),g3(0:,m),g4(0:,m))

	       g1(nTE+1:nc,mc)=g1(nTE+1:nc,m)
	       g2(nTE+1:nc,mc)=g2(nTE+1:nc,m)
	       g3(nTE+1:nc,mc)=g3(nTE+1:nc,m)
	       g4(nTE+1:nc,mc)=g4(nTE+1:nc,m)
	end if


c ============ parete m=mc  ============
	m=mc
       call b_wall("a",CMPLX(xw(0:,m) ,yw(0:,m) )
     x                ,anyx(0:,m),anyy(0:,m),nLE+1,nTE        
     x                ,a(:,m),u(:,m),v(:,m),s(:,m)            
     x                ,g1(0:,m),g2(0:,m),g3(0:,m),g4(0:,m))

c ============= uscita =========
	n=nc

      call b_exit("a",CMPLX(xw(nc,0:) ,yw(nc,0:) )
     x               ,anxx(nc,0:),anxy(nc,0:),1,mc,pexit         
     x               ,a(nc,:),u(nc,:),v(nc,:),s(nc,:)            
     x               ,f1(nc,0:),f2(nc,0:),f3(nc,0:),f4(nc,0:))

c ============= ingresso ===============
	n=0
      call b_inlt("b",CMPLX(xw(0,0:) ,yw(nc,0:) )
     x               ,anxx(0,0:),anxy(0,0:),1,mc
     x               ,T0_i,s_i,ami,sig        
     x               ,a(1,:),u(1,:),v(1,:),s(1,:)            
     x               ,f1(0,0:),f2(0,0:),f3(0,0:),f4(0,0:))

	return
	end
