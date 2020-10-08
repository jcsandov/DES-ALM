
subroutine rhs_viscous
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate the viscous terms for the momentum equations
  ! in all three directions

  ! input
  !     ren (from input file)
  !     xnut(ijk) (from les or rans model)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     q(2-4,ijk)

  ! output
  !     visc(3,ijk)
 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local arrays
  !
  real (kind = rdf), dimension(:,:,:), allocatable :: ret

  ! this could be calculated once and stored too!
  !
  real (kind = rdf) :: rei

  real (kind = rdf) :: dc4, de4, dz4

  ! local, dummy variables for all directions
  !
  real (kind = rdf) :: vtj
  real (kind = rdf) :: a1, a2, a3, b1, b2, b3, c1, c2, c3

  ! local, dummy variables for csi direction
  !
  real (kind = rdf) :: ucs, uet, uze
  real (kind = rdf) :: uu1, uu2, uu3
  real (kind = rdf) :: c11, c12, c13

  ! local, dummy variables for eta direction
  !
  real (kind = rdf) :: vcs, vet, vze
  real (kind = rdf) :: vv1, vv2, vv3
  real (kind = rdf) :: c21, c22, c23
  
  ! local, dummy variables for zet direction
  !
  real (kind = rdf) :: wcs, wet, wze
  real (kind = rdf) :: ww1, ww2, ww3
  real (kind = rdf) :: c31, c32, c33

  ! initialize to zero
  ! 
  allocate (ret(il:iu,jl:ju,kl:ku))

  fv   = zero
  visc = zero
  ret  = zero

!  put in a unchanging module with ren 
  
  rei = 1. / ren

  ! the following would also apply to
  ! RANS and URANS codes
  ! 
  if ( turbulence ) then
     if ( nlinc ) then
        do k = kl, ku
        do j = jl, ju
        do i = il, iu
           ret(i,j,k) = rei
        end do
        end do
        end do
     else
        do k = kl, ku
        do j = jl, ju
        do i = il, iu
           ret(i,j,k) = rei + xnut(i,j,k)
        end do
        end do
        end do
     end if
  else
     do k = kl, ku
     do j = jl, ju
     do i = il, iu
        ret(i,j,k) = rei
     end do
     end do
     end do
  end if

  ! dcsi --> dc
  ! dcs2 --> dc4
  ! dum1 --> de or dz
  ! dum2 --> de4 or dz4

  dc4 = one_fourth * dc
  de4 = one_fourth * de
  dz4 = one_fourth * dz

  ! csi-direction
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 1, im - 1

  do k = k_mysta,     k_myend
  do j = j_mysta,     j_myend
  do i = i_mysta - 1, i_myend

     ucs =  dc * (q(2,i+1,j,k) - q(2,i,j,k))
     vcs =  dc * (q(3,i+1,j,k) - q(3,i,j,k))
     wcs =  dc * (q(4,i+1,j,k) - q(4,i,j,k))

     uet = de4 * (q(2,i+1,j+1,k) - q(2,i+1,j-1,k) + q(2,i,j+1,k) - q(2,i,j-1,k))
     vet = de4 * (q(3,i+1,j+1,k) - q(3,i+1,j-1,k) + q(3,i,j+1,k) - q(3,i,j-1,k))
     wet = de4 * (q(4,i+1,j+1,k) - q(4,i+1,j-1,k) + q(4,i,j+1,k) - q(4,i,j-1,k))

     uze = dz4 * (q(2,i+1,j,k+1) - q(2,i+1,j,k-1) + q(2,i,j,k+1) - q(2,i,j,k-1))
     vze = dz4 * (q(3,i+1,j,k+1) - q(3,i+1,j,k-1) + q(3,i,j,k+1) - q(3,i,j,k-1))
     wze = dz4 * (q(4,i+1,j,k+1) - q(4,i+1,j,k-1) + q(4,i,j,k+1) - q(4,i,j,k-1))

     a1  = csi(1,i+1,j,k) + csi(1,i,j,k)
     a2  = csi(2,i+1,j,k) + csi(2,i,j,k)
     a3  = csi(3,i+1,j,k) + csi(3,i,j,k)

     b1  = eta(1,i+1,j,k) + eta(1,i,j,k)
     b2  = eta(2,i+1,j,k) + eta(2,i,j,k)
     b3  = eta(3,i+1,j,k) + eta(3,i,j,k)

     c1  = zet(1,i+1,j,k) + zet(1,i,j,k)
     c2  = zet(2,i+1,j,k) + zet(2,i,j,k)
     c3  = zet(3,i+1,j,k) + zet(3,i,j,k)

     uu1 = a1 * ucs + b1 * uet + c1 * uze
     vv1 = a1 * vcs + b1 * vet + c1 * vze
     ww1 = a1 * wcs + b1 * wet + c1 * wze

     uu2 = a2 * ucs + b2 * uet + c2 * uze
     vv2 = a2 * vcs + b2 * vet + c2 * vze
     ww2 = a2 * wcs + b2 * wet + c2 * wze

     uu3 = a3 * ucs + b3 * uet + c3 * uze
     vv3 = a3 * vcs + b3 * vet + c3 * vze
     ww3 = a3 * wcs + b3 * wet + c3 * wze

     c11 = a1 * a1 + a2 * a2 + a3 * a3
     c12 = a1 * b1 + a2 * b2 + a3 * b3
     c13 = a1 * c1 + a2 * c2 + a3 * c3

     vtj = one_fourth * (ret(i+1,j,k) + ret(i,j,k)) / (aj(i+1,j,k) + aj(i,j,k))

     fv(1,i,j,k) = vtj * ( a1 * uu1 +  a2 * vv1 +  a3 * ww1 + &
                          c11 * ucs + c12 * uet + c13 * uze)
     fv(2,i,j,k) = vtj * ( a1 * uu2 +  a2 * vv2 +  a3 * ww2 + &
                          c11 * vcs + c12 * vet + c13 * vze)
     fv(3,i,j,k) = vtj * ( a1 * uu3 +  a2 * vv3 +  a3 * ww3 + &
                          c11 * wcs + c12 * wet + c13 * wze)
  end do
  end do
  end do

  ! compute d F_v / d (csi), this may need unrolling
  ! No need for zeroing visc array because it is overwritten here
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1
  
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     visc(:,i,j,k) = dc * (fv(:,i,j,k) - fv(:,i-1,j,k))
  end do
  end do
  end do


  ! eta direction
!   do k = 2, km - 1
!   do j = 1, jm - 1
!   do i = 2, im - 1

  do k = k_mysta,   k_myend
  do j = j_mysta-1, j_myend
  do i = i_mysta,   i_myend

     ucs = dc4 * (q(2,i+1,j+1,k) - q(2,i-1,j+1,k) + q(2,i+1,j,k) - q(2,i-1,j,k))
     vcs = dc4 * (q(3,i+1,j+1,k) - q(3,i-1,j+1,k) + q(3,i+1,j,k) - q(3,i-1,j,k))
     wcs = dc4 * (q(4,i+1,j+1,k) - q(4,i-1,j+1,k) + q(4,i+1,j,k) - q(4,i-1,j,k))

     uet =  de * (q(2,i,j+1,k) - q(2,i,j,k))
     vet =  de * (q(3,i,j+1,k) - q(3,i,j,k))
     wet =  de * (q(4,i,j+1,k) - q(4,i,j,k))

     uze = dz4 * (q(2,i,j+1,k+1) - q(2,i,j+1,k-1) + q(2,i,j,k+1) - q(2,i,j,k-1))
     vze = dz4 * (q(3,i,j+1,k+1) - q(3,i,j+1,k-1) + q(3,i,j,k+1) - q(3,i,j,k-1))
     wze = dz4 * (q(4,i,j+1,k+1) - q(4,i,j+1,k-1) + q(4,i,j,k+1) - q(4,i,j,k-1))

     a1 = csi(1,i,j+1,k) + csi(1,i,j,k)
     a2 = csi(2,i,j+1,k) + csi(2,i,j,k)
     a3 = csi(3,i,j+1,k) + csi(3,i,j,k)

     b1 = eta(1,i,j+1,k) + eta(1,i,j,k)
     b2 = eta(2,i,j+1,k) + eta(2,i,j,k)
     b3 = eta(3,i,j+1,k) + eta(3,i,j,k)

     c1 = zet(1,i,j+1,k) + zet(1,i,j,k)
     c2 = zet(2,i,j+1,k) + zet(2,i,j,k)
     c3 = zet(3,i,j+1,k) + zet(3,i,j,k)

     uu1 = a1 * ucs + b1 * uet + c1 * uze
     vv1 = a1 * vcs + b1 * vet + c1 * vze
     ww1 = a1 * wcs + b1 * wet + c1 * wze

     uu2 = a2 * ucs + b2 * uet + c2 * uze
     vv2 = a2 * vcs + b2 * vet + c2 * vze
     ww2 = a2 * wcs + b2 * wet + c2 * wze

     uu3 = a3 * ucs + b3 * uet + c3 * uze
     vv3 = a3 * vcs + b3 * vet + c3 * vze
     ww3 = a3 * wcs + b3 * wet + c3 * wze

     c21 = b1 * a1 + b2 * a2 + b3 * a3
     c22 = b1 * b1 + b2 * b2 + b3 * b3
     c23 = b1 * c1 + b2 * c2 + b3 * c3

     vtj = one_fourth * (ret(i,j+1,k) + ret(i,j,k)) / (aj(i,j+1,k) + aj(i,j,k))

     fv(1,i,j,k) = vtj * ( b1 * uu1 +  b2 * vv1 +  b3 * ww1 + &
                          c21 * ucs + c22 * uet + c23 * uze)
     fv(2,i,j,k) = vtj * ( b1 * uu2 +  b2 * vv2 +  b3 * ww2 + &
                          c21 * vcs + c22 * vet + c23 * vze)
     fv(3,i,j,k) = vtj * ( b1 * uu3 +  b2 * vv3 +  b3 * ww3 + &
                          c21 * wcs + c22 * wet + c23 * wze)
  end do
  end do
  end do
  
  ! compute d F_v / d (eta)
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(:,i,j,k) = visc(:,i,j,k) + de * (fv(:,i,j,k) - fv(:,i,j-1,k))
  end do
  end do
  end do

  ! zeta-direction
!   do k = 1, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta-1, k_myend
  do j = j_mysta,   j_myend
  do i = i_mysta,   i_myend
     ucs = dc4 * (q(2,i+1,j,k+1) - q(2,i-1,j,k+1) + q(2,i+1,j,k) - q(2,i-1,j,k))
     vcs = dc4 * (q(3,i+1,j,k+1) - q(3,i-1,j,k+1) + q(3,i+1,j,k) - q(3,i-1,j,k))
     wcs = dc4 * (q(4,i+1,j,k+1) - q(4,i-1,j,k+1) + q(4,i+1,j,k) - q(4,i-1,j,k))

     uet = de4 * (q(2,i,j+1,k+1) - q(2,i,j-1,k+1) + q(2,i,j+1,k) - q(2,i,j-1,k))
     vet = de4 * (q(3,i,j+1,k+1) - q(3,i,j-1,k+1) + q(3,i,j+1,k) - q(3,i,j-1,k))
     wet = de4 * (q(4,i,j+1,k+1) - q(4,i,j-1,k+1) + q(4,i,j+1,k) - q(4,i,j-1,k))

     uze =  dz * (q(2,i,j,k+1) - q(2,i,j,k))
     vze =  dz * (q(3,i,j,k+1) - q(3,i,j,k))
     wze =  dz * (q(4,i,j,k+1) - q(4,i,j,k))

     a1 = csi(1,i,j,k+1) + csi(1,i,j,k)
     a2 = csi(2,i,j,k+1) + csi(2,i,j,k)
     a3 = csi(3,i,j,k+1) + csi(3,i,j,k)

     b1 = eta(1,i,j,k+1) + eta(1,i,j,k)
     b2 = eta(2,i,j,k+1) + eta(2,i,j,k)
     b3 = eta(3,i,j,k+1) + eta(3,i,j,k)

     c1 = zet(1,i,j,k+1) + zet(1,i,j,k)
     c2 = zet(2,i,j,k+1) + zet(2,i,j,k)
     c3 = zet(3,i,j,k+1) + zet(3,i,j,k)

     uu1 = a1 * ucs + b1 * uet + c1 * uze
     vv1 = a1 * vcs + b1 * vet + c1 * vze
     ww1 = a1 * wcs + b1 * wet + c1 * wze

     uu2 = a2 * ucs + b2 * uet + c2 * uze
     vv2 = a2 * vcs + b2 * vet + c2 * vze
     ww2 = a2 * wcs + b2 * wet + c2 * wze

     uu3 = a3 * ucs + b3 * uet + c3 * uze
     vv3 = a3 * vcs + b3 * vet + c3 * vze
     ww3 = a3 * wcs + b3 * wet + c3 * wze

     c31 = c1 * a1 + c2 * a2 + c3 * a3
     c32 = c1 * b1 + c2 * b2 + c3 * b3
     c33 = c1 * c1 + c2 * c2 + c3 * c3

     vtj = one_fourth * (ret(i,j,k+1) + ret(i,j,k)) / (aj(i,j,k+1) + aj(i,j,k))

     fv(1,i,j,k) = vtj * ( c1 * uu1 +  c2 * vv1 +  c3 * ww1 + &
                          c31 * ucs + c32 * uet + c33 * uze)
     fv(2,i,j,k) = vtj * ( c1 * uu2 +  c2 * vv2 +  c3 * ww2 + &
                          c31 * vcs + c32 * vet + c33 * vze)
     fv(3,i,j,k) = vtj * ( c1 * uu3 +  c2 * vv3 +  c3 * ww3 + &
                          c31 * wcs + c32 * wet + c33 * wze)
  end do
  end do
  end do

  ! compute d F_v / d (zet)
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     visc(:,i,j,k) = visc(:,i,j,k) + dz * (fv(:,i,j,k) - fv(:,i,j,k-1))
  end do
  end do
  end do

  deallocate (ret)

end subroutine rhs_viscous












