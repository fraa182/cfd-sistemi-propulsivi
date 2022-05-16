subroutine flux_right_tilde
    ! This subroutine is focused on the computation of the fluxes across the "vertical" (the right edge of each element) fluid interfaces inside the domain
    ! The two elements on the two sides are called A and B: A is (n,m); B is (n+1,m)
    use variabili
    implicit none
    integer::n,m
    real::utilde_A,utilde_B,vtilde_A,vtilde_B, F1A,F2A,F3A,F4A,F1B,F2B,F3B,F4B

    ! Perform a double loop for all fluid interfaces (right side of the element) inside the domain:
    ! 1<=n<=nc-1, 1<=m<=mc (BC will be computed in another subroutine)

    do m=1,mc
    
        do n=1,nc-1

            ! consider the normal of the upper interface n,m: compute the projection of the velocity vector on the normal and tangent for
            ! the element A and B
            ! Use a rotation matrix M to get (utilde,vtilde) from (u,v)

			utilde_A=u(n,m)*nx_right(n,m)+v(n,m)*ny_right(n,m)
			vtilde_A=-u(n,m)*ny_right(n,m)+v(n,m)*nx_right(n,m)
			
			utilde_B=u(n+1,m)*nx_right(n+1,m)+v(n+1,m)*ny_right(n+1,m)
			vtilde_B=-u(n+1,m)*ny_right(n+1,m)+v(n+1,m)*nx_right(n+1,m)

            ! Compute the fluxes of the conservative variables across the interface using data from the element A (F1A,F2A,F3A,F4A) and
            ! element B (F1B,F2B,F3B,F4B) (use here the "tilde" variables in the rotated frame of reference)
            ! Example for mass: F1A=u1(n,m)*utilde_A
            !                   F1B=u1(n+1,m)*utilde_B
            ! Go on with the other fluxes

			F1A=u1(n,m)*utilde_A*length_right(n,m)
			F1B=u1(n+1,m)*utilde_B*length_right(n+1,m)
			
			F2A=u1(n,m)*utilde_A**2+u1(n,m)*T(n,m)*length_right(n,m)
			F2B=u1(n+1,m)*utilde_B**2+u1(n+1,m)*T(n+1,m)*length_right(n+1,m)
			
			F3A=u1(n,m)*utilde_A*vtilde_A*length_right(n,m)
			F3B=u1(n+1,m)*utilde_B*vtilde_B*length_right(n+1,m)
			
			F4A=(u4(n,m)+u1(n,m)*T(n,m))*utilde_A*length_right(n,m)
			F4B=(u4(n+1,m)+u1(n+1,m)*T(n+1,m))*utilde_B*length_right(n+1,m)

            ! Compute the fluxes as average between the fluxes computed from side A and side B of the interface
            ! Example for the mass: F1right(n,m)= 0.5*(F1A+F1B)
            ! Go on with the other fluxes

			F1right(n,m)=0.5*(F1A+F1B)
			F2right(n,m)=0.5*(F2A+F2B)
			F3right(n,m)=0.5*(F3A+F3B)
			F4right(n,m)=0.5*(F4A+F4B)

            ! Rotate the fluxes back from the local to the global frame of reference (only momentum fluxes need rotation!)
            ! Use the inverse of M (computed at the beginning of the subroutine)
            
            F2right(n,m)=F2right(n,m)*nx_right(n,m)-F3right(n,m)*ny_right(n,m)
            F3right(n,m)=F2right(n,m)*ny_right(n,m)+F3right(n,m)*nx_right(n,m)

        end do
        
    end do

end subroutine
