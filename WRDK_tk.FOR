        SUBROUTINE WDKs_tk(it,iflg, TITLE, xw,yw,kn,km,nc,mc,x,y,a,u,
     _                v,s)
c        include 'fl2d.inc'
        character titolo*60
        character*8 strin
        character*5 grid
        character*40 junkfl, gridfl

        character*(*) TITLE
        dimension xw(0:kn,0:km), yw(0:kn,0:km)
        real, dimension (kn,km) :: x,y,a,u,v,s

                gamma=1.4
			  ga=gamma/(gamma-1.)

       write (strin,fmt='(i5)') it

        junkfl= TITLE//strin(1:5)//'.plt'
        gridfl=TITLE//'_grd.plt'


	   !if(it.eq.0) then
!	     open(44,file=junkfl,status="unknown") !,access='append')
	   !else
         !  open(44,file=junkfl,status="old",access='append')
	     !open(44,file=junkfl,status="unknown") !,access='append')
	   !end if

	  titolo='TITLE='//'"'//TITLE//'"'

      if(iflg.eq.1) then
	     open(44,file=junkfl,status="unknown") !,access='append')
		 rewind(44)
	     open(44,file=junkfl,status="unknown") !,access='append')
		 rewind(44)
!        if (kflg.eq.0) rewind (unit=44)
	  write(44,'(a60)') titolo
        write(44,*) 'VARIABLES="X","Y","P","M","U","V","TT","S"'

	  ig=it
        write (grid,fmt='(i5)') ig
        write (44,*) 'ZONE T="' ,
     _                GRID , '",' ,
     _               ',I=' , nc,
     _               ',J=' , mc,
     _               ',F="POINT"'

          do m=1,mc
	      do n= 1,nc
              TT=(a(n,m)**2/gamma
     x          +.5*(u(n,m)**2+v(n,m)**2)/ga)

              write (44,*) x(n,m), y(n,m),
     _              (a(n,m)**2/gamma)**ga/exp(s(n,m)/(gamma-1.)),
     _               sqrt(u(n,m)**2+(v(n,m))**2)/a(n,m),
     _               u(n,m),(v(n,m)),
     _               TT,s(n,m)
    	     end do
	    end do
        close(44)
c............................
c
	else if(iflg.eq.-3) then
	  open(45,file=gridfl,status="unknown")
	  rewind(45)
	  ib=1
        write(45,'(a60)') titolo

        write(45,*) 'VARIABLES="X","Y"'

	  ib=ib+1
        write (grid,fmt='(i5)') ib
           write (45,*) 'ZONE T="' ,
     _                   GRID , '",' ,
     _                  ',I=' , nc+1,
     _                  ',J=' , mc+1,
     _                  ',F="POINT"'

            do m=0,mc
	        do n= 0,nc
	          write (45,*)xw(n,m),yw(n,m)
	         end do
	        end do
       close(45)
      end if
        return
        end
c
