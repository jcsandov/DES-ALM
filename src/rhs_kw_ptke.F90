
subroutine rhs_kw_ptke ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Generalized, nonorthogonal curvilinear coordinates

  ! Calculate the mean rate of rotation and strain tensors
  ! original des needs only the mean rate of rotation 
  ! some modified version of des also needs the mean strain rate

  ! only works if dc = de = dz = 1.0 ; that is, on the fine grid only

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local variables

  real (kind = rdf) :: s11, s22, s33
  real (kind = rdf) :: s12, s23, s13
  real (kind = rdf) :: w12, w23, w13

  real (kind = rdf) :: uc, vc, wc
  real (kind = rdf) :: ue, ve, we
  real (kind = rdf) :: uz, vz, wz
  real (kind = rdf) :: u1, u2, u3
  real (kind = rdf) :: v1, v2, v3
  real (kind = rdf) :: w1, w2, w3
  real (kind = rdf) :: sc

  real (kind = rdf) :: dc2, de2, dz2

  integer :: ic

  dc2 = pt5 * dc
  de2 = pt5 * de
  dz2 = pt5 * dz

  ! interior nodes only
!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  if ( .not. nlinc ) then

     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
     do i = i_mysta, i_myend

        uc = dc2 * (q(2,i+1,j,k) - q(2,i-1,j,k))
        vc = dc2 * (q(3,i+1,j,k) - q(3,i-1,j,k))
        wc = dc2 * (q(4,i+1,j,k) - q(4,i-1,j,k))

        ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
        ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
        we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

        uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
        vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
        wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

        u1 = uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
        u2 = uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
        u3 = uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)

        v1 = vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
        v2 = vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
        v3 = vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

        w1 = wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
        w2 = wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
        w3 = wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

        ptke(i,j,k) = xnut(i,j,k) * (two * (u1*u1 + v2*v2 + w3*w3) + &
                     (u2+v1)*(u2+v1) + &
                     (u3+w1)*(u3+w1) + &
                     (v3+w2)*(v3+w2) )

     end do
     end do
     end do


  else ! nlinc == .true.

     do k = k_mysta, k_myend
     do j = j_mysta, j_myend
     do i = i_mysta, i_myend

        uc = dc2 * (q(2,i+1,j,k) - q(2,i-1,j,k))
        vc = dc2 * (q(3,i+1,j,k) - q(3,i-1,j,k))
        wc = dc2 * (q(4,i+1,j,k) - q(4,i-1,j,k))

        ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
        ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
        we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

        uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
        vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
        wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

        u1 = uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
        u2 = uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
        u3 = uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)

        v1 = vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
        v2 = vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
        v3 = vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

        w1 = wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
        w2 = wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
        w3 = wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

        ptke(i,j,k) = xnut(i,j,k) * (two * (u1*u1 + v2*v2 + w3*w3) + &
                     (u2+v1)*(u2+v1) + &
                     (u3+w1)*(u3+w1) + &
                     (v3+w2)*(v3+w2) )

        s11 = u1
        s22 = v2
        s33 = w3

        s12 = pt5 * (u2 + v1)
        s13 = pt5 * (u3 + w1)
        s23 = pt5 * (v3 + w2)

        w12 = pt5 * (u2 - v1)
        w13 = pt5 * (u3 - w1)
        w23 = pt5 * (v3 - w2)

        ! mean strain rate
        so(i,j,k) = s11 * s11 + s22 * s22 + s33 * s33 + &
                   (s12 * s12 + s13 * s13 + s23 * s23 ) * two

        ! mean rate of rotation
        wo(i,j,k) = two * (w12 * w12 + w13 * w13 + w23 * w23)

     end do
     end do
     end do

     ! csi boundaries
     ! 
     if ( myback == mpi_proc_null ) then

        i = i_mysta - 1
        ic = -1
        sc = -dc

        do k = k_mysta, k_myend
        do j = j_mysta, j_myend

           uc = sc  * (q(2,i  ,j,k) - q(2,i-ic,j,k))
           vc = sc  * (q(3,i  ,j,k) - q(3,i-ic,j,k))
           wc = sc  * (q(4,i  ,j,k) - q(4,i-ic,j,k))

           ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
           ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
           we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

           uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
           vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
           wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

           u1 = uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
           u2 = uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
           u3 = uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)

           v1 = vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
           v2 = vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
           v3 = vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

           w1 = wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
           w2 = wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
           w3 = wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

           s11 = u1
           s22 = v2
           s33 = w3

           s12 = pt5 * (u2 + v1)
           s13 = pt5 * (u3 + w1)
           s23 = pt5 * (v3 + w2)

           w12 = pt5 * (u2 - v1)
           w13 = pt5 * (u3 - w1)
           w23 = pt5 * (v3 - w2)

           ! mean strain rate
           so(i,j,k) = s11 * s11 + s22 * s22 + s33 * s33 + &
                (s12 * s12 + s13 * s13 + s23 * s23 ) * two

           ! mean rate of rotation
           wo(i,j,k) = two * (w12 * w12 + w13 * w13 + w23 * w23)

        end do
        end do
     end if

     if ( myfront == mpi_proc_null ) then

        i = i_myend + 1
        ic = 1
        sc = dc

        do k = k_mysta, k_myend
        do j = j_mysta, j_myend

           uc = sc  * (q(2,i  ,j,k) - q(2,i-ic,j,k))
           vc = sc  * (q(3,i  ,j,k) - q(3,i-ic,j,k))
           wc = sc  * (q(4,i  ,j,k) - q(4,i-ic,j,k))

           ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
           ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
           we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

           uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
           vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
           wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

           u1 = uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
           u2 = uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
           u3 = uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)

           v1 = vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
           v2 = vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
           v3 = vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

           w1 = wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
           w2 = wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
           w3 = wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

           s11 = u1
           s22 = v2
           s33 = w3

           s12 = pt5 * (u2 + v1)
           s13 = pt5 * (u3 + w1)
           s23 = pt5 * (v3 + w2)

           w12 = pt5 * (u2 - v1)
           w13 = pt5 * (u3 - w1)
           w23 = pt5 * (v3 - w2)

           ! mean strain rate
           so(i,j,k) = s11 * s11 + s22 * s22 + s33 * s33 + &
                (s12 * s12 + s13 * s13 + s23 * s23 ) * two

           ! mean rate of rotation
           wo(i,j,k) = two * (w12 * w12 + w13 * w13 + w23 * w23)

        end do
        end do
     end if
  end if

end subroutine rhs_kw_ptke


!   if (nlinc) then

!      do i = 1, im, im-1
!         ic = 1
!         sc = dc
!         if (i == 1) then
!            ic = -1
!            sc = -dc
!         end if
!      do k = 2, km-1
!      do j = 2, jm-1

!         uc = sc  * (q(2,i  ,j,k) - q(2,i-ic,j,k))
!         vc = sc  * (q(3,i  ,j,k) - q(3,i-ic,j,k))
!         wc = sc  * (q(4,i  ,j,k) - q(4,i-ic,j,k))

!         ue = de2 * (q(2,i,j+1,k) - q(2,i,j-1,k))
!         ve = de2 * (q(3,i,j+1,k) - q(3,i,j-1,k))
!         we = de2 * (q(4,i,j+1,k) - q(4,i,j-1,k))

!         uz = dz2 * (q(2,i,j,k+1) - q(2,i,j,k-1))
!         vz = dz2 * (q(3,i,j,k+1) - q(3,i,j,k-1))
!         wz = dz2 * (q(4,i,j,k+1) - q(4,i,j,k-1))

!         u1 = uc * csi(1,i,j,k) + ue * eta(1,i,j,k) + uz * zet(1,i,j,k)
!         u2 = uc * csi(2,i,j,k) + ue * eta(2,i,j,k) + uz * zet(2,i,j,k)
!         u3 = uc * csi(3,i,j,k) + ue * eta(3,i,j,k) + uz * zet(3,i,j,k)
          
!         v1 = vc * csi(1,i,j,k) + ve * eta(1,i,j,k) + vz * zet(1,i,j,k)
!         v2 = vc * csi(2,i,j,k) + ve * eta(2,i,j,k) + vz * zet(2,i,j,k)
!         v3 = vc * csi(3,i,j,k) + ve * eta(3,i,j,k) + vz * zet(3,i,j,k)

!         w1 = wc * csi(1,i,j,k) + we * eta(1,i,j,k) + wz * zet(1,i,j,k)
!         w2 = wc * csi(2,i,j,k) + we * eta(2,i,j,k) + wz * zet(2,i,j,k)
!         w3 = wc * csi(3,i,j,k) + we * eta(3,i,j,k) + wz * zet(3,i,j,k)

!         s11 = u1
!         s22 = v2
!         s33 = w3

!         s12 = pt5 * (u2 + v1)
!         s13 = pt5 * (u3 + w1)
!         s23 = pt5 * (v3 + w2)

!         w12 = pt5 * (u2 - v1)
!         w13 = pt5 * (u3 - w1)
!         w23 = pt5 * (v3 - w2)

!         ! mean strain rate
!         so(i,j,k) = s11 * s11 + s22 * s22 + s33 * s33 + &
!                    (s12 * s12 + s13 * s13 + s23 * s23 ) * two

!         ! mean rate of rotation
!         wo(i,j,k) = two * (w12 * w12 + w13 * w13 + w23 * w23)

!      end do
!      end do
!      end do
!   end if



