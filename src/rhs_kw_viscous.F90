
subroutine rhs_kw_viscous
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coordinates

  ! Calculate the viscous terms for the turbulence transport equations
  ! in all three directions

  ! input
  ! ren (from input file)
  ! xnut(ijk) (from des or rans model)
  ! csi(3,ijk), eta(3,ijk), zet(3,iji)
  ! aj(ijk)
  ! q(5,ijk)

  ! output
  ! visc(4,ijk) : 1 = x-mom, 2 = y-mom, 3 = z-mom, 4 = des transport

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local arrays
  real (kind = rdf), dimension(2,il:iu,jl:ju,kl:ku) :: ret !, fv

  ! local variables
  real (kind = rdf) :: rei
  real (kind = rdf) :: dc4, de4, dz4

  real (kind = rdf) :: vtj1, vtj2
  real (kind = rdf) :: vcs1, vet1, vze1, vcs2, vet2, vze2

  real (kind = rdf) :: a1, a2, a3
  real (kind = rdf) :: b1, b2, b3
  real (kind = rdf) :: c1, c2, c3
  real (kind = rdf) :: cc1, cc2, cc3

  ! dum --> dc4 or de4 or dz4 
  dc4 = pt25 * dc
  de4 = pt25 * de
  dz4 = pt25 * dz


  ! put in a unchanging module with ren
  rei= one / ren

!   do k = 1, km
!   do j = 1, jm
!   do i = 1, im

  do k = kl, ku
  do j = jl, ju
  do i = il, iu
     ret(1,i,j,k) = rei + xnut(i,j,k) / sk
     ret(2,i,j,k) = rei + xnut(i,j,k) / sw
  end do
  end do
  end do

  ! csi direction

!   do k = 2, km -1
!   do j = 2, jm -1
!   do i = 1, im -1

  do k = k_mysta,     k_myend
  do j = j_mysta,     j_myend
  do i = i_mysta - 1, i_myend

     vcs1 =  dc * (q(5,i+1,j  ,k  ) - q(5,i  ,j  ,k  ))
     vet1 = de4 * (q(5,i+1,j+1,k  ) - q(5,i+1,j-1,k  ) + &
                   q(5,i  ,j+1,k  ) - q(5,i  ,j-1,k  ))
     vze1 = dz4 * (q(5,i+1,j  ,k+1) - q(5,i+1,j  ,k-1) + &
                   q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k-1))

     vcs2 =  dc * (q(6,i+1,j  ,k  ) - q(6,i  ,j  ,k  ))
     vet2 = de4 * (q(6,i+1,j+1,k  ) - q(6,i+1,j-1,k  ) + &
                   q(6,i  ,j+1,k  ) - q(6,i  ,j-1,k  ))
     vze2 = dz4 * (q(6,i+1,j  ,k+1) - q(6,i+1,j  ,k-1) + &
                   q(6,i  ,j  ,k+1) - q(6,i  ,j  ,k-1))

     a1  = csi(1,i+1,j,k) + csi(1,i,j,k)
     a2  = csi(2,i+1,j,k) + csi(2,i,j,k)
     a3  = csi(3,i+1,j,k) + csi(3,i,j,k)

     b1  = eta(1,i+1,j,k) + eta(1,i,j,k)
     b2  = eta(2,i+1,j,k) + eta(2,i,j,k)
     b3  = eta(3,i+1,j,k) + eta(3,i,j,k)

     c1  = zet(1,i+1,j,k) + zet(1,i,j,k)
     c2  = zet(2,i+1,j,k) + zet(2,i,j,k)
     c3  = zet(3,i+1,j,k) + zet(3,i,j,k)

     cc1 = a1 * a1 + a2 * a2 + a3 * a3
     cc2 = a1 * b1 + a2 * b2 + a3 * b3
     cc3 = a1 * c1 + a2 * c2 + a3 * c3

     vtj1 = pt25 * (ret(1,i+1,j,k) + ret(1,i,j,k)) / (aj(i+1,j,k) + aj(i,j,k))
     vtj2 = pt25 * (ret(2,i+1,j,k) + ret(2,i,j,k)) / (aj(i+1,j,k) + aj(i,j,k))

     fv(1,i,j,k) = vtj1 * (cc1 * vcs1 + cc2 * vet1 + cc3 * vze1)
     fv(2,i,j,k) = vtj2 * (cc1 * vcs2 + cc2 * vet2 + cc3 * vze2)

  end do
  end do
  end do

  ! compute d F_v_t / d (csi)
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(1,i,j,k) = fv(1,i,j,k) - fv(1,i-1,j,k)
     visc(2,i,j,k) = fv(2,i,j,k) - fv(2,i-1,j,k)
  end do
  end do
  end do
 
  ! eta direction

!   do k = 2, km - 1
!   do j = 1, jm - 1
!   do i = 2, im - 1

  do k = k_mysta,     k_myend
  do j = j_mysta - 1, j_myend
  do i = i_mysta,     i_myend

     vcs1 = dc4 * (q(5,i+1,j+1,k  ) - q(5,i-1,j+1,k  ) + &
                   q(5,i+1,j  ,k  ) - q(5,i-1,j  ,k  ))
     vet1 = de  * (q(5,i  ,j+1,k  ) - q(5,i  ,j  ,k  ))
     vze1 = dz4 * (q(5,i  ,j+1,k+1) - q(5,i  ,j+1,k-1) + &
                   q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k-1))

     vcs2 = dc4 * (q(6,i+1,j+1,k  ) - q(6,i-1,j+1,k  ) + &
                   q(6,i+1,j  ,k  ) - q(6,i-1,j  ,k  ))
     vet2 = de  * (q(6,i  ,j+1,k  ) - q(6,i  ,j  ,k  ))
     vze2 = dz4 * (q(6,i  ,j+1,k+1) - q(6,i  ,j+1,k-1) + &
                   q(6,i  ,j  ,k+1) - q(6,i  ,j  ,k-1))

     a1  = csi(1,i,j+1,k) + csi(1,i,j,k)
     a2  = csi(2,i,j+1,k) + csi(2,i,j,k)
     a3  = csi(3,i,j+1,k) + csi(3,i,j,k)

     b1  = eta(1,i,j+1,k) + eta(1,i,j,k)
     b2  = eta(2,i,j+1,k) + eta(2,i,j,k)
     b3  = eta(3,i,j+1,k) + eta(3,i,j,k)

     c1  = zet(1,i,j+1,k) + zet(1,i,j,k)
     c2  = zet(2,i,j+1,k) + zet(2,i,j,k)
     c3  = zet(3,i,j+1,k) + zet(3,i,j,k)

     cc1 = b1 * a1 + b2 * a2 + b3 * a3
     cc2 = b1 * b1 + b2 * b2 + b3 * b3
     cc3 = b1 * c1 + b2 * c2 + b3 * c3

     vtj1 = pt25 * (ret(1,i,j+1,k) + ret(1,i,j,k)) / (aj(i,j+1,k) + aj(i,j,k))
     vtj2 = pt25 * (ret(2,i,j+1,k) + ret(2,i,j,k)) / (aj(i,j+1,k) + aj(i,j,k))

     fv(1,i,j,k) = vtj1 * (cc1 * vcs1 + cc2 * vet1 + cc3 * vze1)
     fv(2,i,j,k) = vtj2 * (cc1 * vcs2 + cc2 * vet2 + cc3 * vze2)

  end do
  end do
  end do

  ! compute d F_v_t / d (eta)
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(1,i,j,k) = visc(1,i,j,k) + fv(1,i,j,k) - fv(1,i,j-1,k)
     visc(2,i,j,k) = visc(2,i,j,k) + fv(2,i,j,k) - fv(2,i,j-1,k)
  end do
  end do
  end do
 
  ! zeta-direction terms
!   do k = 1, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta - 1, k_myend
  do j = j_mysta,     j_myend
  do i = i_mysta,     i_myend

     vcs1 = dc4 * (q(5,i+1,j  ,k+1) - q(5,i-1,j  ,k+1) + &
                   q(5,i+1,j  ,k  ) - q(5,i-1,j  ,k  ))
     vet1 = de4 * (q(5,i  ,j+1,k+1) - q(5,i  ,j-1,k+1) + &
                   q(5,i  ,j+1,k  ) - q(5,i  ,j-1,k  ))
     vze1 = dz  * (q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k  ))

     vcs2 = dc4 * (q(6,i+1,j  ,k+1) - q(6,i-1,j  ,k+1) + &
                   q(6,i+1,j  ,k  ) - q(6,i-1,j  ,k  ))
     vet2 = de4 * (q(6,i  ,j+1,k+1) - q(6,i  ,j-1,k+1) + &
                   q(6,i  ,j+1,k  ) - q(6,i  ,j-1,k  ))
     vze2 = dz  * (q(6,i  ,j  ,k+1) - q(6,i  ,j  ,k  ))

     a1  = csi(1,i,j,k+1) + csi(1,i,j,k)
     a2  = csi(2,i,j,k+1) + csi(2,i,j,k)
     a3  = csi(3,i,j,k+1) + csi(3,i,j,k)

     b1  = eta(1,i,j,k+1) + eta(1,i,j,k)
     b2  = eta(2,i,j,k+1) + eta(2,i,j,k)
     b3  = eta(3,i,j,k+1) + eta(3,i,j,k)

     c1  = zet(1,i,j,k+1) + zet(1,i,j,k)
     c2  = zet(2,i,j,k+1) + zet(2,i,j,k)
     c3  = zet(3,i,j,k+1) + zet(3,i,j,k)

     cc1 = c1 * a1 + c2 * a2 + c3 * a3
     cc2 = c1 * b1 + c2 * b2 + c3 * b3
     cc3 = c1 * c1 + c2 * c2 + c3 * c3

     vtj1 = pt25 * (ret(1,i,j,k+1) + ret(1,i,j,k)) / (aj(i,j,k+1) + aj(i,j,k))
     vtj2 = pt25 * (ret(2,i,j,k+1) + ret(2,i,j,k)) / (aj(i,j,k+1) + aj(i,j,k))

     fv(1,i,j,k) = vtj1 * (cc1 * vcs1 + cc2 * vet1 + cc3 * vze1)
     fv(2,i,j,k) = vtj2 * (cc1 * vcs2 + cc2 * vet2 + cc3 * vze2)

  end do
  end do
  end do

  ! compute d F_v_t / d (zet)
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(1,i,j,k) = visc(1,i,j,k) + fv(1,i,j,k) - fv(1,i,j,k-1)
     visc(2,i,j,k) = visc(2,i,j,k) + fv(2,i,j,k) - fv(2,i,j,k-1)
  end do
  end do
  end do
 

end subroutine rhs_kw_viscous


