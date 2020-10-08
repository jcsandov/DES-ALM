
subroutine rhs_ns_rstress
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! routine to compute the Reynolds-stress gradient terms       
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local arrays
  ! 
  real (kind = rdf), dimension(1:3,il:iu,jl:ju,kl:ku) :: fv

  ! local, dummy variables
  !
  real (kind = rdf)  :: au1, au2, au3
  real (kind = rdf)  :: av1, av2, av3
  real (kind = rdf)  :: aw1, aw2, aw3
  real (kind = rdf)  :: gx, gy, gz, ajc
 
  !
  ! csi-derivative terms
  !
!   do k = 2, km-1
!   do j = 2, jm-1
!   do i = 1, im-1

  rstr  = zero
  fv = zero

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta-1, i_myend

     au1 = uij(1,i+1,j,k) + uij(1,i,j,k)
     av1 = uij(4,i+1,j,k) + uij(4,i,j,k)
     aw1 = uij(5,i+1,j,k) + uij(5,i,j,k)

     au2 = uij(4,i+1,j,k) + uij(4,i,j,k)
     av2 = uij(2,i+1,j,k) + uij(2,i,j,k)
     aw2 = uij(6,i+1,j,k) + uij(6,i,j,k)

     au3 = uij(5,i+1,j,k) + uij(5,i,j,k)
     av3 = uij(6,i+1,j,k) + uij(6,i,j,k)
     aw3 = uij(3,i+1,j,k) + uij(3,i,j,k)

     gx  = csi(1,i+1,j,k) + csi(1,i,j,k)
     gy  = csi(2,i+1,j,k) + csi(2,i,j,k)
     gz  = csi(3,i+1,j,k) + csi(3,i,j,k)

     ajc = aj(i+1,j,k) + aj(i,j,k)

     fv(1,i,j,k) = pt5 * (gx*au1 + gy*av1 + gz*aw1) / ajc
     fv(2,i,j,k) = pt5 * (gx*au2 + gy*av2 + gz*aw2) / ajc
     fv(3,i,j,k) = pt5 * (gx*au3 + gy*av3 + gz*aw3) / ajc

  end do
  end do
  end do

  ! derivative
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend


     rstr(1,i,j,k) = dc * (fv(1,i,j,k) - fv(1,i-1,j,k))
     rstr(2,i,j,k) = dc * (fv(2,i,j,k) - fv(2,i-1,j,k))
     rstr(3,i,j,k) = dc * (fv(3,i,j,k) - fv(3,i-1,j,k))

  end do
  end do
  end do

  !
  ! Eta-momentum
  !
!   do k = 2, km-1
!   do j = 1, jm-1
!   do i = 2, im-1

  do k = k_mysta, k_myend
  do j = j_mysta-1, j_myend
  do i = i_mysta, i_myend

     au1 = uij(1,i,j+1,k) + uij(1,i,j,k)
     av1 = uij(4,i,j+1,k) + uij(4,i,j,k)
     aw1 = uij(5,i,j+1,k) + uij(5,i,j,k)

     au2 = uij(4,i,j+1,k) + uij(4,i,j,k)
     av2 = uij(2,i,j+1,k) + uij(2,i,j,k)
     aw2 = uij(6,i,j+1,k) + uij(6,i,j,k)

     au3 = uij(5,i,j+1,k) + uij(5,i,j,k)
     av3 = uij(6,i,j+1,k) + uij(6,i,j,k)
     aw3 = uij(3,i,j+1,k) + uij(3,i,j,k)

     gx  = eta(1,i,j+1,k) + eta(1,i,j,k)
     gy  = eta(2,i,j+1,k) + eta(2,i,j,k)
     gz  = eta(3,i,j+1,k) + eta(3,i,j,k)

     ajc = aj(i,j+1,k)+aj(i,j,k)

     fv(1,i,j,k) = pt5 * (gx*au1 + gy*av1 + gz*aw1)/ajc
     fv(2,i,j,k) = pt5 * (gx*au2 + gy*av2 + gz*aw2)/ajc
     fv(3,i,j,k) = pt5 * (gx*au3 + gy*av3 + gz*aw3)/ajc

  end do
  end do
  end do

  ! derivative

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend


!   do k = 2, km-1
!   do j = 2, jm-1
!   do i = 2, im-1

     rstr(1,i,j,k) = rstr(1,i,j,k) + de * (fv(1,i,j,k) - fv(1,i,j-1,k))
     rstr(2,i,j,k) = rstr(2,i,j,k) + de * (fv(2,i,j,k) - fv(2,i,j-1,k))
     rstr(3,i,j,k) = rstr(3,i,j,k) + de * (fv(3,i,j,k) - fv(3,i,j-1,k))

  end do
  end do
  end do

  !
  ! Zet-momentum
  !
!   do k = 1, km-1
!   do j = 2, jm-1
!   do i = 2, im-1

  do k = k_mysta-1, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend


     au1 = uij(1,i,j,k+1) + uij(1,i,j,k)
     av1 = uij(4,i,j,k+1) + uij(4,i,j,k)
     aw1 = uij(5,i,j,k+1) + uij(5,i,j,k)

     au2 = uij(4,i,j,k+1) + uij(4,i,j,k)
     av2 = uij(2,i,j,k+1) + uij(2,i,j,k)
     aw2 = uij(6,i,j,k+1) + uij(6,i,j,k)

     au3 = uij(5,i,j,k+1) + uij(5,i,j,k)
     av3 = uij(6,i,j,k+1) + uij(6,i,j,k)
     aw3 = uij(3,i,j,k+1) + uij(3,i,j,k)

     gx  = zet(1,i,j,k+1) + zet(1,i,j,k)
     gy  = zet(2,i,j,k+1) + zet(2,i,j,k)
     gz  = zet(3,i,j,k+1) + zet(3,i,j,k)

     ajc = aj(i,j,k+1)+aj(i,j,k)

     fv(1,i,j,k) = pt5 * (gx*au1 + gy*av1 + gz*aw1) / ajc
     fv(2,i,j,k) = pt5 * (gx*au2 + gy*av2 + gz*aw2) / ajc
     fv(3,i,j,k) = pt5 * (gx*au3 + gy*av3 + gz*aw3) / ajc

  end do
  end do
  end do

  ! derivative

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

!   do k = 2, km-1
!   do j = 2, jm-1
!   do i = 2, im-1

     rstr(1,i,j,k) = rstr(1,i,j,k) + dz * (fv(1,i,j,k) - fv(1,i,j,k-1))
     rstr(2,i,j,k) = rstr(2,i,j,k) + dz * (fv(2,i,j,k) - fv(2,i,j,k-1))
     rstr(3,i,j,k) = rstr(3,i,j,k) + dz * (fv(3,i,j,k) - fv(3,i,j,k-1))

  end do
  end do
  end do

end subroutine rhs_ns_rstress
