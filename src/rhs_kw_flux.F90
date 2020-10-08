
subroutine rhs_kw_flux ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, non-orthogonal, curvlinear coordinates

  ! Calculate the convective terms of scalar transport equations
  ! using flux difference splitting first-order upwind differences.

  ! input
  !     ucn_j(3,ijk)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     q(5:6,ijk)
  !     dc, de, dz

  ! output
  !     rh(ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local dummy array
  !real (kind = rdf), dimension(2,im,jm,km) :: ft

  ! local variable
  real (kind = rdf) :: gcu, gev, gzw
  real (kind = rdf) :: eigh

  !
  ! csi direction
  !

  ! calculate the convective flux at i+1/2
  !
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 1, im - 1

  do k = k_mysta,     k_myend
  do j = j_mysta,     j_myend
  do i = i_mysta - 1, i_myend

     gcu = (csi(1,i+1,j,k) + csi(1,i,j,k)) * (q(2,i+1,j,k) + q(2,i,j,k))
     gev = (csi(2,i+1,j,k) + csi(2,i,j,k)) * (q(3,i+1,j,k) + q(3,i,j,k))
     gzw = (csi(3,i+1,j,k) + csi(3,i,j,k)) * (q(4,i+1,j,k) + q(4,i,j,k))

     eigh = pt5 * abs((gcu + gev + gzw) / (aj(i+1,j,k) + aj(i,j,k)))

     fv(1,i,j,k) = pt5 * ( q(5,i+1,j,k) * ucn_j(1,i+1,j,k) + &
                           q(5,i  ,j,k) * ucn_j(1,i  ,j,k) - &
                           eigh * (q(5,i+1,j,k) - q(5,i,j,k)) )

     fv(2,i,j,k) = pt5 * ( q(6,i+1,j,k) * ucn_j(1,i+1,j,k) + &
                           q(6,i  ,j,k) * ucn_j(1,i  ,j,k) - &
                           eigh * (q(6,i+1,j,k) - q(6,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at i
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(1,i,j,k) = fv(1,i,j,k) - fv(1,i-1,j,k)
     rh(2,i,j,k) = fv(2,i,j,k) - fv(2,i-1,j,k)
  end do
  end do
  end do

  !
  ! eta direction
  !

  ! calculate the convective flux at j+1/2
!   do k = 2, km - 1
!   do j = 1, jm - 1
!   do i = 2, im - 1

  do k = k_mysta,     k_myend
  do j = j_mysta - 1, j_myend
  do i = i_mysta,     i_myend

     gcu = (eta(1,i,j+1,k) + eta(1,i,j,k)) * (q(2,i,j+1,k) + q(2,i,j,k))
     gev = (eta(2,i,j+1,k) + eta(2,i,j,k)) * (q(3,i,j+1,k) + q(3,i,j,k))
     gzw = (eta(3,i,j+1,k) + eta(3,i,j,k)) * (q(4,i,j+1,k) + q(4,i,j,k))

     eigh = pt5 * abs((gcu + gev + gzw) / (aj(i,j+1,k) + aj(i,j,k)))

     fv(1,i,j,k) = pt5 * ( q(5,i,j+1,k) * ucn_j(2,i,j+1,k) + &
                           q(5,i,j  ,k) * ucn_j(2,i,j  ,k) - &
                           eigh * (q(5,i,j+1,k) - q(5,i,j,k)) )

     fv(2,i,j,k) = pt5 * ( q(6,i,j+1,k) * ucn_j(2,i,j+1,k) + &
                           q(6,i,j  ,k) * ucn_j(2,i,j  ,k) - &
                           eigh * (q(6,i,j+1,k) - q(6,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at j
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(1,i,j,k) = rh(1,i,j,k) + fv(1,i,j,k) - fv(1,i,j-1,k)
     rh(2,i,j,k) = rh(2,i,j,k) + fv(2,i,j,k) - fv(2,i,j-1,k)
  end do
  end do
  end do

  !
  ! zet direction
  !

  ! calculate the convective flux at k+1/2
!   do k = 1, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta - 1, k_myend
  do j = j_mysta,     j_myend
  do i = i_mysta,     i_myend

     gcu = (zet(1,i,j,k+1) + zet(1,i,j,k)) * (q(2,i,j,k+1) + q(2,i,j,k))
     gev = (zet(2,i,j,k+1) + zet(2,i,j,k)) * (q(3,i,j,k+1) + q(3,i,j,k))
     gzw = (zet(3,i,j,k+1) + zet(3,i,j,k)) * (q(4,i,j,k+1) + q(4,i,j,k))

     eigh = pt5 * abs((gcu + gev + gzw) / (aj(i,j,k+1) + aj(i,j,k)))

     fv(1,i,j,k) = pt5 * ( q(5,i,j,k+1) * ucn_j(3,i,j,k+1) + &
                           q(5,i,j,k  ) * ucn_j(3,i,j,k  ) - &
                           eigh * (q(5,i,j,k+1) - q(5,i,j,k)) )

     fv(2,i,j,k) = pt5 * ( q(6,i,j,k+1) * ucn_j(3,i,j,k+1) + &
                           q(6,i,j,k  ) * ucn_j(3,i,j,k  ) - &
                           eigh * (q(6,i,j,k+1) - q(6,i,j,k)) )

  end do
  end do
  end do

  ! calculate the convective flux at k
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     rh(1,i,j,k) = rh(1,i,j,k) + fv(1,i,j,k) - fv(1,i,j,k-1)
     rh(2,i,j,k) = rh(2,i,j,k) + fv(2,i,j,k) - fv(2,i,j,k-1)
  end do
  end do
  end do

end subroutine rhs_kw_flux



