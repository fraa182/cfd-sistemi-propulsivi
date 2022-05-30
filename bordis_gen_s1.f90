     subroutine b_wall(Given,zw,anx,any,nf,nl            &
	                        ,a_lr,u_lr,v_lr,s_lr         &
					        ,f1,f2,f3,f4)

	    !implicit double precision (a,b,d-h,o-y)
        INTEGER, PARAMETER           :: rdp=4
	    
        
		CHARACTER*1 Given
		COMPLEX(rdp) zw(0:*)

		REAL(KIND=rdp), DIMENSION (6)          :: Uag, Ubg, Ucg, Fg
		REAL(KIND=rdp), DIMENSION (*)          :: a_lr, u_lr, v_lr, s_lr !a, u, v, s,
		REAL(KIND=rdp), DIMENSION (0:*)        :: anx,any, f1,f2,f3,f4   !, xw, yw
        
		COMPLEX(rdp) znor
        
		gamma=1.4d0
		gd=(gamma-1.d0)/2.d0

        SELECT CASE (given)
		  CASE ("a")
!       -------------------------
!         PARETE  SUPERIORE
!       -------------------------
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sa=s_lr(n)   !s(n)  +s_lr(n)*(islop)
              aa=a_lr(n)   !a(n)  +a_lr(n)*(islop)
              udum=u_lr(n) !u(n)+u_lr(n)*(islop)
              vdum=v_lr(n) !v(n)+v_lr(n)*(islop)
              ua=udum*anx_+vdum*any_
              va=udum*any_-vdum*anx_

              R3a=aa/gd+ua
	        
	          uc=0.d0
	          ac=(R3a-uc)*gd
              aw=ac
              uw=uc
              vw=va
              sw=sa

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_)
	          dzz=abs(zw(n-1)-zw(n))
	          
			  call flusso (znor,Ucg,Fg)
              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
	         END DO
		  CASE ("b")
!       -------------------------
!          PARETE  INFERIORE
!       -------------------------
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sb=s_lr(n)   !s(n)   +s_lr(n)*(islop)
			  ab=a_lr(n)   !a(n)   +a_lr(n)*(islop)
              udum=u_lr(n) !u(n) +u_lr(n)*(islop)
              vdum=v_lr(n) !v(n) +v_lr(n)*(islop)
              ub=udum*anx_+vdum*any_
              vb=udum*any_-vdum*anx_

              R1b=ab/gd-ub 

              uc=0.d0
	          aw=(R1b+uc)*gd
              uw=uc
              vw=vb
              sw=sb

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_,rdp)
	          dzz=abs(zw(n-1)-zw(n))
 	          call flusso (znor,Ucg,Fg)

              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
             END DO
     END SELECT     

   END SUBROUTINE b_wall 
!
!
!
     subroutine b_exit(Given,zw,anx,any,nf,nl,pexit  &
	                        ,a_lr,u_lr,v_lr,s_lr         &
					        ,f1,f2,f3,f4)

	    !implicit double precision (a,b,d-h,o-y)
        INTEGER, PARAMETER           :: rdp=4
	    
        
		CHARACTER*1 Given
		COMPLEX(rdp) zw(0:*)
		
		REAL(KIND=rdp), DIMENSION (6)          :: Uag, Ubg, Ucg, Fg
		REAL(KIND=rdp), DIMENSION (*)          :: a_lr, u_lr, v_lr, s_lr !a, u, v, s,
		REAL(KIND=rdp), DIMENSION (0:*)        :: anx,any, f1,f2,f3,f4   !, xw, yw
        
		COMPLEX(rdp) znor
        
		gamma=1.4d0
		gd=(gamma-1.d0)/2.d0
        ga=gamma/(gamma-1.d0)


        SELECT CASE (given)
		  CASE ("a")
!       -------------------------
!         EXIT piu
!       -------------------------
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sa=s_lr(n)   !s(n)  +s_lr(n)*(islop)
              aa=a_lr(n)   !a(n)  +a_lr(n)*(islop)
              udum=u_lr(n) !u(n)+u_lr(n)*(islop)
              vdum=v_lr(n) !v(n)+v_lr(n)*(islop)
              ua=udum*anx_+vdum*any_
              va=udum*any_-vdum*anx_

               if (ua.gt.aa)then
                !        .... uscita supersonica
                ac=aa
                uc=ua
                vc=va
                sc=sa
               else
              !         uscita subsonica (pressione imposta)
                R3a=aa/gd+ua
                sc=sa
                vc=va
                pc=pexit
                ac=sqrt(gamma*pc**(1./ga)*exp(sc/gamma))
                uc=R3a-ac/gd
              endif
              aw=ac
              uw=uc
              vw=va
              sw=sa

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_,rdp)
	          dzz=abs(zw(n-1)-zw(n))
	          
			  call flusso (znor,Ucg,Fg)
              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
	         END DO
		  CASE ("b")
!       -------------------------
!          EXIT meno
!       -------------------------
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sb=s_lr(n)   !s(n)   +s_lr(n)*(islop)
			  ab=a_lr(n)   !a(n)   +a_lr(n)*(islop)
              udum=u_lr(n) !u(n) +u_lr(n)*(islop)
              vdum=v_lr(n) !v(n) +v_lr(n)*(islop)
              ub=udum*anx_+vdum*any_
              vb=udum*any_-vdum*anx_

               if (ub.gt.ab)then
                !        .... uscita supersonica
                ac=ab
                uc=ub
                vc=vb
                sc=sb
               else
              !         uscita subsonica (pressione imposta)
                R1b=ab/gd-ub
                sc=sb
                vc=vb
                pc=pexit
                ac=sqrt(gamma*pc**(1./ga)*exp(sc/gamma))
                uc=R1b+ac/gd
              endif
	          aw=ac
              uw=uc
              vw=vb
              sw=sb

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_,rdp)
	          dzz=abs(zw(n-1)-zw(n))
 	          call flusso (znor,Ucg,Fg)

              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
             END DO
     END SELECT     

   END SUBROUTINE b_exit 
!
!
!
   SUBROUTINE b_inlt(Given,zw,anx,any,nf,nl       &
	                      ,T0_i,s_i,am_i,sig      &
	                      ,a_lr,u_lr,v_lr,s_lr    &
					      ,f1,f2,f3,f4)

	    !implicit double precision (a,b,d-h,o-y)
        INTEGER, PARAMETER           :: rdp=4
	    
        
		CHARACTER*1 Given
		COMPLEX(rdp) zw(0:*)

		REAL(KIND=rdp), DIMENSION (*)          :: T0_i, s_i, sig, am_i 
		!REAL(KIND=rdp), DIMENSION (*)          :: a_w, u_w, v_w, s_w 
		REAL(KIND=rdp), DIMENSION (6)          :: Uag, Ubg, Ucg, Fg
		REAL(KIND=rdp), DIMENSION (*)          :: a_lr, u_lr, v_lr, s_lr !a, u, v, s,
		REAL(KIND=rdp), DIMENSION (0:*)        :: anx,any, f1,f2,f3,f4   !, xw, yw
        
		COMPLEX(rdp) znor
        
		gamma=1.4d0
		gd=(gamma-1.d0)/2.d0
        ga=gamma/(gamma-1.d0)


        SELECT CASE (given)
		  CASE ("a")
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sa=s_lr(n)   !s(n)  +s_lr(n)*(islop)
              aa=a_lr(n)   !a(n)  +a_lr(n)*(islop)
              udum=u_lr(n) !u(n)+u_lr(n)*(islop)
              vdum=v_lr(n) !v(n)+v_lr(n)*(islop)
              ua=udum*anx_+vdum*any_
              va=udum*any_-vdum*anx_
	          dum1=any_-sig(n)*anx_
	          dum2=anx_+sig(n)*any_
	          fsig=dum1/dum2
              sc=s_i(n)
			  Tc=T0_i(n)
               if (am_i(n).gt.1.)then
                !        .... ingresso supersonico
                ac=sqrt(gamma*Tc/(1.d0+gd*am_i(n)**2))
                qc=ac*am_i(n)
                uc=qc/sqrt(fsig)
                vc=uc*fsig
               else
              !         ingresso subsonico
                R3a=aa/gd+ua
	            ddddddddddd=exp((sb-sc)/(2.d0*gamma))/gd
	            e1=1.+d**2*gd*(1.d0+fsig)
	            e2=(gamma-1.d0)*r1b*d*(1.d0+fsig)
	            e3=gd*r1b**2*(1.d0+fsig)-gamma*Tc
	             if((e2**2-4.*e1*e3) .le. 0.) then
	              istatus=-10
	              e1=0. !return
	             end if
   	            ac=(e2+sqrt(e2**2-4.*e1*e3))/(2*e1)
	            uccccccccc=d*ac-R3a
               endif
              aw=ac
              uw=uc
              vw=vc
              sw=sc

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_,rdp)
	          dzz=abs(zw(n-1)-zw(n))
	          
			  call flusso (znor,Ucg,Fg)
              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
	         END DO
!
		  CASE ("b")
             DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sb=s_lr(n)   !s(n)   +s_lr(n)*(islop)
			  ab=a_lr(n)   !a(n)   +a_lr(n)*(islop)
              udum=u_lr(n) !u(n) +u_lr(n)*(islop)
              vdum=v_lr(n) !v(n) +v_lr(n)*(islop)
              ub=udum*anx_+vdum*any_
              vb=udum*any_-vdum*anx_
	          dum1=any_-sig(n)*anx_
	          dum2=anx_+sig(n)*any_
	          fsig=dum1/dum2
              sc=s_i(n)
			  Tc=T0_i(n)
               if (am_i(n)*cos(atan(sig(n))).gt.1.)then
                !        .... ingresso supersonico
                ac=sqrt(gamma*Tc/(1.d0+gd*am_i(n)**2))
                qc=ac*am_i(n)
                uc=qc/sqrt(1.d0+fsig**2)
               else
              !         ingresso subsonico
                R1b=ab/gd-ub
	            d=exp((sb-sc)/(2.d0*gamma))/gd
	            e1=1.+d**2*gd*(1.d0+fsig)
	            e2=(gamma-1.d0)*r1b*d*(1.d0+fsig**2)
	            e3=gd*r1b**2*(1.d0+fsig**2)-gamma*Tc
	             if((e2**2-4.d0*e1*e3) .le. 0.) then
	              istatus=-10
	              e1=0. !return
	             end if
   	            ac=(e2+sqrt(e2**2-4.d0*e1*e3))/(2.d0*e1)
	            uc=d*ac-R1b
               endif
	          
			  vc=uc*fsig
	          aw=ac
              uw=uc
              vw=vc
              sw=sc

	          Ucg(1)=Aw
	          Ucg(2)=Uw
	          Ucg(3)=Vw
	          Ucg(4)=Sw

	          znor=cmplx(anx_,any_,rdp)
	          dzz=abs(zw(n-1)-zw(n))
 	          call flusso (znor,Ucg,Fg)

              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz

              ucc=uw*anx_+vw*any_
              vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
			 END DO
     END SELECT     

   END SUBROUTINE b_inlt 
!
!
!
   SUBROUTINE b_peri(Given,zw,anx,any,nf,nl            &
	                      ,a_lr1,u_lr1,v_lr1,s_lr1     &
	                      ,a_lr2,u_lr2,v_lr2,s_lr2     &
	                      ,f1,f2,f3,f4)

	    !implicit double precision (a,b,d-h,o-y)
	    
        !IMPLICIT NONE
        INTEGER, PARAMETER           :: rdp=4
       
		COMPLEX(rdp) zw(0:*)

        INTEGER                                 :: nf, nl, nGauss, kn
		REAL(KIND=rdp), DIMENSION (*)           :: a_lr1, u_lr1, v_lr1, s_lr1 !a, u, v, s,
		REAL(KIND=rdp), DIMENSION (*)           :: a_lr2, u_lr2, v_lr2, s_lr2 !a, u, v, s,
		REAL(KIND=rdp), DIMENSION (0:*)         :: anx,any, f1,f2,f3,f4   !, xw, yw


		REAL(KIND=rdp), DIMENSION (6)        :: Uag, Ubg, Ucg, Fg
        
		CHARACTER*1                          :: Given
        INTEGER                              :: n,n1
		REAL(KIND=rdp)                       :: gamma, gd


		COMPLEX(rdp) znor
        
		gamma=1.4d0
		gd=(gamma-1.d0)/2.d0

           DO n=nf,nl
              anx_=anx(n)
              any_=any(n)
              sa  =s_lr2(n)
              aa=  a_lr2(n)
              udum=u_lr2(n)
              vdum=v_lr2(n)
              ua=udum*anx_+vdum*any_
              va=udum*any_-vdum*anx_

              sb  =s_lr1(n)
              ab  =a_lr1(n)
              udum=u_lr1(n)
              vdum=v_lr1(n)
              ub=udum*anx_+vdum*any_
              vb=udum*any_-vdum*anx_
	        
              Uag(1)=Aa
              Uag(2)=Ua
              Uag(3)=Va
              Uag(4)=Sa

              Ubg(1)=Ab
              Ubg(2)=Ub
              Ubg(3)=Vb
              Ubg(4)=Sb

	          znor=cmplx(anx_,any_,rdp)
	          
              om=0.d0
              Fg=0.
			  !call riema (znor,om,Uag,Ubg,Ucg,Fg)
	          
			 dzz=abs(zw(n-1)-zw(n))
			  
              f1(n)=Fg(1)*dzz
              f2(n)=Fg(2)*dzz
              f3(n)=Fg(3)*dzz
              f4(n)=Fg(4)*dzz
 
	          !Aw = Ucg(1)
	          !Uw = Ucg(2)
	          !Vw = Ucg(3)
	          !Sw = Ucg(4)

              !ucc=uw*anx_+vw*any_
              !vcc=uw*any_-vw*anx_
              !a_W(n)=aw
              !u_W(n)=ucc
              !v_W(n)=vcc
              !s_W(n)=sw
	         END DO
  
  END SUBROUTINE b_peri
! ****
! ****
        SUBROUTINE flusso(znor,Uwg,Fw) 
!        include 'fl2d.inc'
         ga()=gamma/(gamma-1.)
         gb()=1./(gamma-1.)

        real                           ::gamma=1.4
	  dimension Uwg(*),Fw(*)
	  complex znor

		!gamma=1.4d0

	  anx=real(znor)
	  any=aimag(znor)

	  Aw=Uwg(1)
	  Uw=Uwg(2)
	  Vw=Uwg(3)
	  Sw=Uwg(4)

 	  tw=Aw**2/gamma
	  pw=tw**ga()/exp(sw*gb())
	  rhow=pw/tw
	  ew=rhow*(tw*gb() +.5d0*(Uw**2+Vw**2))
	  Uwc=Uw*anx+Vw*any
	  Vwc=Uw*any-Vw*anx
      
 	  Fw(1)=rhow*Uw
	  Fw(2)=pw*anx+Fw(1)*Uwc
	  Fw(3)=Fw(1)*Vwc+pw*any
	  Fw(4)=Uw*(ew+pw)

        return
        end
