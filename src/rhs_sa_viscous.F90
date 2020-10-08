
subroutine rhs_sa_viscous ()
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
  real (kind = rdf), dimension(:,:,:), allocatable :: ret

  ! local variables
  real (kind = rdf) :: rei
  real (kind = rdf) :: dc4, de4, dz4

  real (kind = rdf) :: vtj
  real (kind = rdf) :: vcs, vet, vze

  real (kind = rdf) :: a1, a2, a3
  real (kind = rdf) :: b1, b2, b3
  real (kind = rdf) :: c1, c2, c3
  real (kind = rdf) :: cc1, cc2, cc3

  ! dum --> dc4 or de4 or dz4 
  dc4 = pt25 * dc
  de4 = pt25 * de
  dz4 = pt25 * dz

  !
  allocate (ret(il:iu,jl:ju,kl:ku))

  ! put in a unchanging module with ren
  rei= one / ren

  do k = kl, ku
  do j = jl, ju
  do i = il, iu
     ret(i,j,k) = (rei + q(5,i,j,k)) / sig
  end do
  end do
  end do

  ! csi direction

  do k = k_mysta,   k_myend
  do j = j_mysta,   j_myend
  do i = i_mysta-1, i_myend

     vcs =  dc * (q(5,i+1,j  ,k  ) - q(5,i  ,j  ,k  ))
     vet = de4 * (q(5,i+1,j+1,k  ) - q(5,i+1,j-1,k  ) + &
                  q(5,i  ,j+1,k  ) - q(5,i  ,j-1,k  ))
     vze = dz4 * (q(5,i+1,j  ,k+1) - q(5,i+1,j  ,k-1) + &
                  q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k-1))

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

     vtj = pt25 * (ret(i+1,j,k) + ret(i,j,k)) / (aj(i+1,j,k) + aj(i,j,k))

     fv(i,j,k) = vtj * (cc1 * vcs + cc2 * vet + cc3 * vze) 

  end do
  end do
  end do

  ! compute d F_v_t / d (csi)
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(i,j,k) = fv(i,j,k) - fv(i-1,j,k)
  end do
  end do
  end do
 
  ! eta direction

  do k = k_mysta,   k_myend
  do j = j_mysta-1, j_myend
  do i = i_mysta,  i_myend

     vcs = dc4 * (q(5,i+1,j+1,k  ) - q(5,i-1,j+1,k  ) + &
                  q(5,i+1,j  ,k  ) - q(5,i-1,j  ,k  ))
     vet = de  * (q(5,i  ,j+1,k  ) - q(5,i  ,j  ,k  ))
     vze = dz4 * (q(5,i  ,j+1,k+1) - q(5,i  ,j+1,k-1) + &
                  q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k-1))

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

     vtj = pt25 * (ret(i,j+1,k) + ret(i,j,k)) / (aj(i,j+1,k) + aj(i,j,k))

     fv(i,j,k) = vtj * (cc1 * vcs + cc2 * vet + cc3 * vze) 

  end do
  end do
  end do

  ! compute d F_v_t / d (eta)
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(i,j,k) = visc(i,j,k) + fv(i,j,k) - fv(i,j-1,k)
  end do
  end do
  end do
 
  ! zeta-direction terms
  do k = k_mysta-1, k_myend
  do j = j_mysta,   j_myend
  do i = i_mysta,   i_myend

     vcs = dc4 * (q(5,i+1,j  ,k+1) - q(5,i-1,j  ,k+1) + &
                  q(5,i+1,j  ,k  ) - q(5,i-1,j  ,k  ))
     vet = de4 * (q(5,i  ,j+1,k+1) - q(5,i  ,j-1,k+1) + &
                  q(5,i  ,j+1,k  ) - q(5,i  ,j-1,k  ))
     vze = dz  * (q(5,i  ,j  ,k+1) - q(5,i  ,j  ,k  ))

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

     vtj = pt25 * (ret(i,j,k+1) + ret(i,j,k)) / (aj(i,j,k+1) + aj(i,j,k))

     fv(i,j,k) = vtj * (cc1 * vcs + cc2 * vet + cc3 * vze) 

  end do
  end do
  end do

  ! compute d F_v_t / d (zet)
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(i,j,k) = visc(i,j,k) + fv(i,j,k) - fv(i,j,k-1)
  end do
  end do
  end do
 
  deallocate (ret)

end subroutine rhs_sa_viscous


