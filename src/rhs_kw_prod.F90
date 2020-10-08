
subroutine rhs_kw_prod
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coordinates

  ! Calculate source terms of kw (low-Reynolds number version, 
  !                               Wilcox, 1994)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  real (kind = rdf) :: pkw1a, pkw2a, pkw2b

!   do k = 2, km-1
!   do j = 2, jm-1
!   do i = 2, im-1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     pkw1a =  Bs(i,j,k) * q(6,i,j,k) * q(5,i,j,k)
     pkw(1,i,j,k) = (-ptke(i,j,k) + pkw1a) / aj(i,j,k)

     pkw2a = -Ac(i,j,k) * q(6,i,j,k) * ptke(i,j,k) / q(5,i,j,k)
     pkw2b =  cw2 * q(6,i,j,k) * q(6,i,j,k)
     pkw(2,i,j,k) = ( pkw2a + pkw2b ) / aj(i,j,k)

  end do
  end do
  end do


end subroutine rhs_kw_prod


