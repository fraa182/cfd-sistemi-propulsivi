 subroutine flux_up_tilde
 ! This subroutine is focused on the computation of the fluxes across the upper edge of each element (fluid interfaces) inside the domain
    ! The two elements on the two sides are called A and B: A is (n,m); B is (n,m+1)
    use variabili
    implicit none
    integer::n,m
    real::utilde_A,utilde_B,vtilde_A,vtilde_B, F1A,F2A,F3A,F4A,F1B,F2B,F3B,F4B

    ! Perform a double loop for all fluid interfaces (upper side of the element) inside the domain:
    ! 1<=n<=nc, 1<=m<=mc-1 (BC will be computed in another subroutine)


    do m=1,mc-1
        do n=1,nc

            ! consider the normal of the upper interface n,m: compute the projection of the velocity vector on the normal and tangent for the element A and B
            ! Use a rotation matrix M to get (utilde,vtilde) from (u,v)

			utilde_A=u(n,m)*nx_right(n,m)+v(n,m)*ny_right(n,m)
			vtilde_A=-u(n,m)*ny_right(n,m)+v(n,m)*nx_right(n,m)
			
			utilde_B=u(n,m+1)*nx_right(n,m+1)+v(n,m+1)*ny_right(n,m+1)
			vtilde_B=-u(n,m+1)*ny_right(n,m+1)+v(n,m+1)*nx_right(n,m+1)

            ! Compute the fluxes of the conservative variables across the interface using data from the element A (F1A,F2A,F3A,F4A) and
            ! element B (F1B,F2B,F3B,F4B)(use here the "tilde" variables in the rotated frame of reference)
            ! Example for mass: F1A=u1(n,m)*utilde_A
            !                   F1B=u1(n,m+1)*utilde_B
            ! Go on with the other fluxes

			F1A=u1(n,m)*utilde_A
			F1B=u1(n,m+1)*utilde_B
			
			F2A=u1(n,m)*utilde_A**2+u1(n,m)*T(n,m)
			F2B=u1(n,m+1)*utilde_B**2+u1(n,m+1)*T(n,m+1)
			
			F3A=u1(n,m)*utilde_A*vtilde_A
			F3B=u1(n,m+1)*utilde_B*vtilde_B
			
			F4A=(u4(n,m)+u1(n,m)*T(n,m))*utilde_A
			F4B=(u4(n,m+1)+u1(n,m+1)*T(n,m+1))*utilde_B

            ! Compute the fluxes as average between the fluxes computed from side A and side B of the interface
            ! Example for the mass: F1up(n,m)= 0.5*(F1A+F1B)
            ! Go on with the other fluxes

			F1right(n,m)=0.5*(F1A+F1B)
			F2right(n,m)=0.5*(F2A+F2B)
			F3right(n,m)=0.5*(F3A+F3B)
			F4right(n,m)=0.5*(F4A+F4B)

            ! Rotate the fluxes back from the local to the global frame of reference (only momentum fluxes need rotation!)
            ! Use the inverse of M (computed at the beginning of the subroutine)

			F2right(n,m)=F2right(n,m)*nx_up(n,m)-F3right(n,m)*ny_up(n,m)
            F3right(n,m)=F2right(n,m)*ny_up(n,m)+F3right(n,m)*nx_up(n,m)

        end do
    end do

end subroutine
