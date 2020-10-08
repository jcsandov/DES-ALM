subroutine brhs_viscous
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates
  
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local arrays, how will you handle these ?
  real (kind = rdf), dimension(:,:,:), allocatable :: fv
  real (kind = rdf), dimension(  :,:), allocatable :: ret

  ! this could be calculated once and stored too!
  real (kind = rdf) :: rei

  real (kind = rdf) ::  de4, dz4

  ! local, dummy variables for all directions
  real (kind = rdf) :: vtj
  real (kind = rdf) :: b1, b2, b3, c1, c2, c3

  ! local, dummy variables for csi direction
  real (kind = rdf) :: uet, uze
  real (kind = rdf) :: uu1, uu2, uu3

  ! local, dummy variables for eta direction
  real (kind = rdf) :: vet, vze
  real (kind = rdf) :: vv1, vv2, vv3
  real (kind = rdf) :: c22, c23

  ! local, dummy variables for zet direction
  real (kind = rdf) :: wet, wze
  real (kind = rdf) :: ww1, ww2, ww3
  real (kind = rdf) :: c32, c33

  ! put in a unchanging module with ren 
  rei = one / ren

  !
  allocate (fv(1:3,jl:ju,kl:ku), &
               ret(jl:ju,kl:ku) )

  ! the following would also apply to
  ! RANS and URANS codes
  if (turbulence) then
     do k = kl, ku
     do j = jl, ju
        ret(j,k) = rei + xnut(j,k)
     end do
     end do
  else
     do k = kl, ku
     do j = jl, ju
        ret(j,k) = rei
     end do
     end do
  end if

  ! spacings
  de4 = pt25 * de
  dz4 = pt25 * dz

  ! eta direction
  do k = k_mysta,   k_myend
  do j = j_mysta-1, j_myend
     uet =  de * (q(2,j+1,k) - q(2,j,k))
     vet =  de * (q(3,j+1,k) - q(3,j,k))
     wet =  de * (q(4,j+1,k) - q(4,j,k))
     uze = dz4 * (q(2,j+1,k+1) - q(2,j+1,k-1) + q(2,j,k+1) - q(2,j,k-1))
     vze = dz4 * (q(3,j+1,k+1) - q(3,j+1,k-1) + q(3,j,k+1) - q(3,j,k-1))
     wze = dz4 * (q(4,j+1,k+1) - q(4,j+1,k-1) + q(4,j,k+1) - q(4,j,k-1))

     b1 = eta(1,j+1,k) + eta(1,j,k)
     b2 = eta(2,j+1,k) + eta(2,j,k)
     b3 = eta(3,j+1,k) + eta(3,j,k)

     c1 = zet(1,j+1,k) + zet(1,j,k)
     c2 = zet(2,j+1,k) + zet(2,j,k)
     c3 = zet(3,j+1,k) + zet(3,j,k)

     uu1 = b1 * uet + c1 * uze
     vv1 = b1 * vet + c1 * vze
     ww1 = b1 * wet + c1 * wze

     uu2 = b2 * uet + c2 * uze
     vv2 = b2 * vet + c2 * vze
     ww2 = b2 * wet + c2 * wze

     uu3 = b3 * uet + c3 * uze
     vv3 = b3 * vet + c3 * vze
     ww3 = b3 * wet + c3 * wze

     c22 = b1 * b1 + b2 * b2 + b3 * b3
     c23 = b1 * c1 + b2 * c2 + b3 * c3

     vtj = one_fourth * (ret(j+1,k) + ret(j,k)) / (aj(j+1,k) + aj(j,k))

     fv(1,j,k) = vtj * (b1 * uu1 + b2 * vv1 + b3 * ww1 + c22 * uet + c23 * uze)
     fv(2,j,k) = vtj * (b1 * uu2 + b2 * vv2 + b3 * ww2 + c22 * vet + c23 * vze)
     fv(3,j,k) = vtj * (b1 * uu3 + b2 * vv3 + b3 * ww3 + c22 * wet + c23 * wze)
  end do
  end do
  
  ! compute d F_v / d (eta)
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     visc(:,j,k) = de * (fv(:,j,k) - fv(:,j-1,k))
  end do
  end do

  ! zeta-direction
  do k = k_mysta-1, k_myend
  do j = j_mysta,   j_myend
     uet = de4 * (q(2,j+1,k+1) - q(2,j-1,k+1) + q(2,j+1,k) - q(2,j-1,k))
     vet = de4 * (q(3,j+1,k+1) - q(3,j-1,k+1) + q(3,j+1,k) - q(3,j-1,k))
     wet = de4 * (q(4,j+1,k+1) - q(4,j-1,k+1) + q(4,j+1,k) - q(4,j-1,k))
     uze =  dz * (q(2,j,k+1) - q(2,j,k))
     vze =  dz * (q(3,j,k+1) - q(3,j,k))
     wze =  dz * (q(4,j,k+1) - q(4,j,k))

     b1 = eta(1,j,k+1) + eta(1,j,k)
     b2 = eta(2,j,k+1) + eta(2,j,k)
     b3 = eta(3,j,k+1) + eta(3,j,k)

     c1 = zet(1,j,k+1) + zet(1,j,k)
     c2 = zet(2,j,k+1) + zet(2,j,k)
     c3 = zet(3,j,k+1) + zet(3,j,k)

     uu1 = b1 * uet + c1 * uze
     vv1 = b1 * vet + c1 * vze
     ww1 = b1 * wet + c1 * wze

     uu2 = b2 * uet + c2 * uze
     vv2 = b2 * vet + c2 * vze
     ww2 = b2 * wet + c2 * wze

     uu3 = b3 * uet + c3 * uze
     vv3 = b3 * vet + c3 * vze
     ww3 = b3 * wet + c3 * wze

     c32 = c1 * b1 + c2 * b2 + c3 * b3
     c33 = c1 * c1 + c2 * c2 + c3 * c3

     vtj = one_fourth * (ret(j,k+1) + ret(j,k)) / (aj(j,k+1) + aj(j,k))

     fv(1,j,k) = vtj * (c1 * uu1 + c2 * vv1 + c3 * ww1 + c32 * uet + c33 * uze)
     fv(2,j,k) = vtj * (c1 * uu2 + c2 * vv2 + c3 * ww2 + c32 * vet + c33 * vze)
     fv(3,j,k) = vtj * (c1 * uu3 + c2 * vv3 + c3 * ww3 + c32 * wet + c33 * wze)
  end do
  end do

  ! compute d F_v / d (zet)
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     visc(:,j,k) = visc(:,j,k) + dz * (fv(:,j,k) - fv(:,j,k-1))
  end do
  end do


  deallocate (fv, ret)


end subroutine brhs_viscous
