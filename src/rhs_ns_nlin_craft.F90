
subroutine rhs_ns_nlin_craft
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Reynolds stresses based on Craft et al.                       
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local parameters
  real (kind = rdf), parameter :: onept33 = 1.33333_rdf
  real (kind = rdf), parameter ::    pt33 = 0.33333_rdf

  ! coefficient for nlinear Craft model
  real (kind = rdf) :: ck, c1, c2, c3, c4, c5, c6, c7

  ! local variables
  real (kind = rdf) :: u11, u22, u33, u12, u13, u23
  real (kind = rdf) :: uc, ue, uz, vc, ve, vz, wc, we, wz
  real (kind = rdf) :: d11, d12, d13, d21, d22, d23, d31, d32, d33
  real (kind = rdf) :: s11, s12, s13, s21, s22, s23, s31, s32, s33
  real (kind = rdf) :: w11, w12, w13, w21, w22, w23, w31, w32, w33
  real (kind = rdf) :: as11, as12, as13, as21, as22, as23, as31, as32, as33
  real (kind = rdf) :: aw11, aw12, aw13, aw21, aw22, aw23, aw31, aw32, aw33
  real (kind = rdf) :: asw11, asw12, asw13
  real (kind = rdf) :: asw21, asw22, asw23
  real (kind = rdf) :: asw31, asw32, asw33
  real (kind = rdf) :: bssw11, bssw12, bssw13
  real (kind = rdf) :: bssw21, bssw22, bssw23
  real (kind = rdf) :: bssw31, bssw32, bssw33

  real (kind = rdf) :: sa, wa, sswa
  real (kind = rdf) :: dk, dsa, dwa, dsswa, xke, xnke
  real (kind = rdf) :: Sinv, Winv, SWmax
  real (kind = rdf) :: tum1, tum2, tum3
  real (kind = rdf) :: cmu, xnke2

  ! dummy variable
  real (kind = rdf) :: dc2, de2, dz2

  integer    &
       ista, &
       jsta, &
       ksta, &
       iend, &
       jend, &
       kend

  dc2 = pt5 * dc
  de2 = pt5 * de
  dz2 = pt5 * dz

  ! coefficient for nonlinear Craft model
  ck = 0.09_rdf
  c1 =-0.10_rdf
  c2 = 0.10_rdf    
  c3 = 0.26_rdf   
  c4 =-1.00_rdf   
  c5 = 0.00_rdf
  c6 =-0.10_rdf
  c7 = 0.10_rdf

  uij = zero

  !
  !
  
  ! include ghost nodes in calculation
  ! 
  ista = i_mysta
  jsta = j_mysta
  ksta = k_mysta
  
  iend = i_myend
  jend = j_myend
  kend = k_myend

  if ( mydown  /= mpi_proc_null ) ksta = ksta - 1
  if ( myup    /= mpi_proc_null ) kend = kend + 1
  if ( myleft  /= mpi_proc_null ) jsta = jsta - 1
  if ( myright /= mpi_proc_null ) jend = jend + 1
  if ( myback  /= mpi_proc_null ) ista = ista - 1
  if ( myfront /= mpi_proc_null ) iend = iend + 1

  do k = ksta, kend
  do j = jsta, jend
  do i = ista, iend

     ! Velocity derivatives
     uc = dc2 * (q(2,i+1,j,k) - q(2,i-1,j,k))
     vc = dc2 * (q(3,i+1,j,k) - q(3,i-1,j,k))
     wc = dc2 * (q(4,i+1,j,k) - q(4,i-1,j,k))

     ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
     ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
     we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

     uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
     vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
     wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

     d11= uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
     d12= uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
     d13= uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)

     d21= vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
     d22= vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
     d23= vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

     d31= wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
     d32= wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
     d33= wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

     !
     ! Sij, Wij tensors
     s11 = d11 + d11                   
     s12 = d12 + d21                 
     s13 = d13 + d31                 
     s21 = s12
     s22 = d22 + d22                  
     s23 = d23 + d32                
     s31 = s13
     s32 = s23
     s33 = d33 + d33                 

     w11 = zero
     w12 = d12 - d21                 
     w13 = d13 - d31                
     w21 =-w12
     w22 = zero
     w23 = d23 - d32                    
     w31 =-w13
     w32 =-w23
     w33 = zero

     as11 = s11*s11 + s12*s21 + s13*s31
     as12 = s11*s12 + s12*s22 + s13*s32
     as13 = s11*s13 + s12*s23 + s13*s33
     as21 = as12
     as22 = s21*s12 + s22*s22 + s23*s32
     as23 = s21*s13 + s22*s23 + s23*s33
     as31 = as13
     as32 = as23
     as33 = s31*s13 + s32*s23 + s33*s33

     asw11 =           w12*s21 + w13*s31
     asw12 =           w12*s22 + w13*s32
     asw13 =           w12*s23 + w13*s33
     asw21 = w21*s11           + w23*s31
     asw22 = w21*s12           + w23*s32
     asw23 = w21*s13           + w23*s33
     asw31 = w31*s11 + w32*s21
     asw32 = w31*s12 + w32*s22
     asw33 = w31*s13 + w32*s23

     aw11 =           w12*w21 + w13*w31
     aw12 =           w12*w22 + w13*w32
     aw13 =           w12*w23 + w13*w33
     aw21 = w21*w11           + w23*w31
     aw22 = w21*w12           + w23*w32
     aw23 = w21*w13           + w23*w33
     aw31 = w31*w11 + w32*w21
     aw32 = w31*w12 + w32*w22
     aw33 = w31*w13 + w32*w23

     bssw11 =            as12*w21 + as13*w31
     bssw12 = as11*w12            + as13*w32
     bssw13 = as11*w13 + as12*w23
     bssw21 =            as22*w21 + as23*w31
     bssw22 = as21*w12            + as23*w32
     bssw23 = as21*w13 + as22*w23
     bssw31 =            as32*w21 + as33*w31
     bssw32 = as31*w12            + as33*w32
     bssw33 = as31*w13 + as32*w23

     sa = as11 + as22 + as33
     wa = aw11 + aw22 + aw33
     sswa= bssw11 + bssw22 + bssw33

     !
     !
     dk    = two * q(5,i,j,k) / three
     dsa   = sa / three
     dwa   = wa / three
     dsswa = two * sswa / three
     xke   = one / q(6,i,j,k) / ck
     xnke  = xnut(i,j,k) * xke

     !
     Sinv = xke * sqrt( sa / two)
     Winv = xke * sqrt(-wa / two)

     !
     SWmax = max(Sinv,Winv)
     tum1 = one + 0.35_rdf * SWmax**onept5
     tum2 = exp(-0.75_rdf*SWmax)
     tum3 = exp(-0.36_rdf/tum2)
 
     cmu  = 0.3_rdf *(one - tum3)/tum1

     xnke2= cmu*xnke*xke

     !
     u11 = dk-xnut(i,j,k)*s11 + &
           xnke* (c1*(as11-dsa)+c2*(asw11+asw11)-c3*(aw11-dwa) ) + &
           xnke2*(c4*(bssw11+bssw11-dsswa)+(c6*sa-c7*wa)*s11)

     u22 = dk-xnut(i,j,k)*s22 + &
           xnke* (c1*(as22-dsa)+c2*(asw22+asw22)-c3*(aw22-dwa) ) + &
           xnke2*(c4*(bssw22+bssw22-dsswa)+(c6*sa-c7*wa)*s22)

     u33 = dk-xnut(i,j,k)*s33 + &
           xnke* (c1*(as33-dsa)+c2*(asw33+asw33)-c3*(aw33-dwa) ) + &
           xnke2*(c4*(bssw33+bssw33-dsswa)+(c6*sa-c7*wa)*s33)

     u12 =   -xnut(i,j,k)*s12 + &
           xnke* (c1*as12+c2*(asw12+asw21)      -c3*aw12 ) + &
           xnke2*(c4*(bssw12+bssw21)+(c6*sa-c7*wa)*s12)

     u13 =   -xnut(i,j,k)*s13 + &
           xnke* (c1*as13+c2*(asw13+asw31)      -c3*aw13 ) + &
           xnke2*(c4*(bssw13+bssw31)+(c6*sa-c7*wa)*s13)

     u23 =   -xnut(i,j,k)*s23 + &
           xnke* (c1*as23+c2*(asw23+asw32)      -c3*aw23 ) + &
           xnke2*(c4*(bssw23+bssw32)+(c6*sa-c7*wa)*s23)

     !
     uij(1,i,j,k) = u11
     uij(2,i,j,k) = u22
     uij(3,i,j,k) = u33
     uij(4,i,j,k) = u12
     uij(5,i,j,k) = u13
     uij(6,i,j,k) = u23

  end do
  end do
  end do

  ! boundary conditions for i extremes only;
  ! using f90 array notation
  ! 
  i = i_mysta - 1
  if (myback == mpi_proc_null) then
     do k = ksta, kend
     do j = jsta, jend
     uij(1:6,i,j,k) = uij(1:6,i+1,j,k)
     end do
     end do
  end if
  
  i = i_myend + 1
  if (myfront == mpi_proc_null) then
     do k = ksta, kend
     do j = jsta, jend
     uij(1:6,i,j,k) = uij(1:6,i-1,j,k)
     !uij(1:6,i,j,k) = onept33 * uij(1:6,i-1,j,k) - pt33 * uij(1:6,i-2,j,k)
     end do
     end do
  end if

!   !
!   ! boundary condition
!   do k = 1, km
!   do j = 1, jm
!      i = 1
!      uij(:,i,j,k) = uij(:,i+1,j,k)
!      i = im
!      uij(:,i,j,k) = onept33 * uij(:,i-1,j,k) - pt33 * uij(:,i-2,j,k)
!   end do
!   end do
  
  ! symmetric
  !   k = km
  !do j = 1, jm
  !do i = 1, im
  !   uij(1,i,j,k) = uij(1,i,j,k-1)
  !   uij(2,i,j,k) = uij(2,i,j,k-1)
  !   uij(3,i,j,k) = uij(3,i,j,k-1)
  !   uij(4,i,j,k) = uij(4,i,j,k-1)
  !   uij(5,i,j,k) = zero
  !   uij(6,i,j,k) = zero
  !end do
  !end do

end subroutine rhs_ns_nlin_craft

