 subroutine flux_up_tilde
 ! This subroutine is focused on the computation of the fluxes across the upper edge of each element (fluid interfaces) inside the domain
    ! The two elements on the two sides are called A and B: A is (n,m); B is (n,m+1)
    use variabili
    implicit none
    integer::n,m
    real::utilde_A,utilde_B,vtilde_A,vtilde_B, F1A,F2A,F3A,F4A,F1B,F2B,F3B,F4B,lambda
    real::Utilde1_A,Utilde2_A,Utilde3_A,Utilde4_A,Utilde1_B,Utilde2_B,Utilde3_B,Utilde4_B,F1corr,F2corr,F3corr,F4corr

    ! Perform a double loop for all fluid interfaces (upper side of the element) inside the domain:
    ! 1<=n<=nc, 1<=m<=mc-1 (BC will be computed in another subroutine)


    do m=1,mc-1
    
        do n=1,nc

            ! consider the normal of the upper interface n,m: compute the projection of the velocity vector on the normal and tangent for the element A and B
            ! Use a rotation matrix M to get (utilde,vtilde) from (u,v)

			utilde_A=u(n,m)*nx_up(n,m)+v(n,m)*ny_up(n,m)
			vtilde_A=-u(n,m)*ny_up(n,m)+v(n,m)*nx_up(n,m)
			
			utilde_B=u(n,m+1)*nx_up(n,m+1)+v(n,m+1)*ny_up(n,m+1)
			vtilde_B=-u(n,m+1)*ny_up(n,m+1)+v(n,m+1)*nx_up(n,m+1)

            ! Compute the fluxes of the conservative variables across the interface using data from the element A (F1A,F2A,F3A,F4A) and
            ! element B (F1B,F2B,F3B,F4B)(use here the "tilde" variables in the rotated frame of reference)
            ! Example for mass: F1A=u1(n,m)*utilde_A
            !                   F1B=u1(n,m+1)*utilde_B
            ! Go on with the other fluxes

			F1A=u1(n,m)*utilde_A*length_up(n,m)
			F1B=u1(n,m+1)*utilde_B*length_up(n,m+1)
			
			F2A=(u1(n,m)*utilde_A**2+u1(n,m)*T(n,m))*length_up(n,m)
			F2B=(u1(n,m+1)*utilde_B**2+u1(n,m+1)*T(n,m+1))*length_up(n,m+1)
			
			F3A=u1(n,m)*utilde_A*vtilde_A*length_up(n,m)
			F3B=u1(n,m+1)*utilde_B*vtilde_B*length_up(n,m+1)
			
			F4A=(u4(n,m)+u1(n,m)*T(n,m))*utilde_A*length_up(n,m)
			F4B=(u4(n,m+1)+u1(n,m+1)*T(n,m+1))*utilde_B*length_up(n,m+1)

            ! Compute the fluxes as average between the fluxes computed from side A and side B of the interface
            ! Example for the mass: F1up(n,m)= 0.5*(F1A+F1B)
            ! Go on with the other fluxes

			F1up(n,m)=0.5*(F1A+F1B)
			F2up(n,m)=0.5*(F2A+F2B)
			F3up(n,m)=0.5*(F3A+F3B)
			F4up(n,m)=0.5*(F4A+F4B)
			
			! Correct the fluxes with numerical viscosity term
			! Calculate local lambda value, as the max of lambda in cell n,m  and lambda in cell n+1,m
			! Calculate the conservatives vector in normal-tangential SR
			! Calculate the flux correction
			lambda=max((vtilde_A+a(n,m)),(vtilde_B+a(n,m+1)))
			
			Utilde1_A=u1(n,m)*length_up(n,m)
			Utilde2_A=u1(n,m)*utilde_A*length_up(n,m)
			Utilde3_A=u1(n,m)*vtilde_A*length_up(n,m)
			Utilde4_A=u4(n,m)*length_up(n,m)
			
			Utilde1_B=u1(n,m+1)*length_up(n,m+1)
			Utilde2_B=u1(n,m+1)*utilde_B*length_up(n,m+1)
			Utilde3_B=u1(n,m+1)*vtilde_B*length_up(n,m+1)
			Utilde4_B=u4(n,m+1)*length_up(n,m+1)
			
			F1corr=0.5*lambda*(Utilde1_B-Utilde1_A)
			F2corr=0.5*lambda*(Utilde2_B-Utilde2_A)
			F3corr=0.5*lambda*(Utilde3_B-Utilde3_A)
			F4corr=0.5*lambda*(Utilde4_B-Utilde4_A)
			
			F1up(n,m)=F1up(n,m)-F1corr
			F2up(n,m)=F2up(n,m)-F2corr
			F3up(n,m)=F3up(n,m)-F3corr
			F4up(n,m)=F4up(n,m)-F4corr

            ! Rotate the fluxes back from the local to the global frame of reference (only momentum fluxes need rotation!)
            ! Use the inverse of M (computed at the beginning of the subroutine)

			F2up(n,m)=F2up(n,m)*nx_up(n,m)-F3up(n,m)*ny_up(n,m)
            F3up(n,m)=F2up(n,m)*ny_up(n,m)+F3up(n,m)*nx_up(n,m)

        end do
        
    end do

end subroutine
