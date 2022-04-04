subroutine grid
    use variabili
    implicit none
    integer::n,m


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

	

    !********************************************************************
    ! Compute node coordinates x(n,m), y(n,m)  with 0<= n <= nc, 0<= m <= mc -> matrici







    !********************************************************************
    ! Compute edge length and normal on right edge nx_right(n,m), ny_right(n,m), length_right(n,m)  with 0<= n <= nc, 1<= m <= mc






    !********************************************************************
    ! Compute edge length and normal on upper edge nx_up(n,m), ny_up(n,m), length_up(n,m)  with 1<= n <= nc, 0<= m <= mc






    !********************************************************************
    ! Compute element area(n,m) and center coordinates xg(n,m), yg(n,m) with 1<= n <= nc, 1<= m <= mc









end subroutine
