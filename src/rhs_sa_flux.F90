
subroutine rhs_sa_flux ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, non-orthogonal, curvlinear coordinates

  ! Calculate the convective terms of scalar transport equations
  ! using flux difference splitting first-order upwind differences.

  ! input
  !     ucn_j(3,ijk)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     q(1:5,ijk)
  !     dc, de, dz

  ! output
  !     rh(ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local variable
  real (kind = rdf) :: gcu, gev, gzw
  real (kind = rdf) :: eigh

  !
  ! csi direction
  !

  ! calculate the convective flux at i+1/2
  do k = k_mysta,   k_myend
  do j = j_mysta,   j_myend
  do i = i_mysta-1, i_myend

     eigh = pt5 * abs(ucn_j(1,i+1,j,k) + ucn_j(1,i,j,k))

     fv(i,j,k) = pt5 * ( q(5,i+1,j,k) * ucn_j(1,i+1,j,k) + &
                         q(5,i  ,j,k) * ucn_j(1,i  ,j,k) - &
                         eigh * (q(5,i+1,j,k) - q(5,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at i
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(i,j,k) = fv(i,j,k) - fv(i-1,j,k)
  end do
  end do
  end do

  !
  ! eta direction
  !

  ! calculate the convective flux at j+1/2
  do k = k_mysta,   k_myend
  do j = j_mysta-1, j_myend
  do i = i_mysta,   i_myend

     eigh = pt5 * abs(ucn_j(2,i,j+1,k) + ucn_j(2,i,j,k))

     fv(i,j,k) = pt5 * ( q(5,i,j+1,k) * ucn_j(2,i,j+1,k) + &
                         q(5,i,j  ,k) * ucn_j(2,i,j  ,k) - &
                         eigh * (q(5,i,j+1,k) - q(5,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at j
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(i,j,k) = rh(i,j,k) + fv(i,j,k) - fv(i,j-1,k)
  end do
  end do
  end do

  !
  ! zet direction
  !

  ! calculate the convective flux at k+1/2
  do k = k_mysta-1, k_myend
  do j = j_mysta,   j_myend
  do i = i_mysta,   i_myend

     eigh = pt5 * abs(ucn_j(3,i,j,k+1) + ucn_j(3,i,j,k))

     fv(i,j,k) = pt5 * ( q(5,i,j,k+1) * ucn_j(3,i,j,k+1) + &
                         q(5,i,j,k  ) * ucn_j(3,i,j,k  ) - &
                         eigh * (q(5,i,j,k+1) - q(5,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at k
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(i,j,k) = rh(i,j,k) + fv(i,j,k) - fv(i,j,k-1)
  end do
  end do
  end do


end subroutine rhs_sa_flux




