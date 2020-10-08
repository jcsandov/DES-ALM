
subroutine rhs_sa_adi_solver ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  ! implicit, diagonal approximate factorization solver
  !
  ! input
  !     rh(4,ijk)
  !     dtau(ijk)
  !     csi(ijk), eta(ijk), zet(ijk)
  !     aj(ijk)
  !     dc, de, dz
  !     dc2, de2, dz2
  !
  ! output
  !     rh(4,ijk)
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  implicit none

  integer :: info

  ! dimension lhs and temporary rhs
  real (kind = rdf), dimension(:), allocatable :: aw, ap, ae, av
  real (kind = rdf), dimension(:), allocatable :: aq
  integer :: va, vb

  ! local variables
  real (kind = rdf) :: rei
  real (kind = rdf) :: dtj
  real (kind = rdf) :: gx, gy, gz, ev, ucna

  integer :: ilength
  integer :: jlength
  integer :: klength

  ! put in a unchanging module with ren 
  !
  rei = one / ren

  ! csi operator
  !
  va=i_mysta
  vb=i_myend
  ilength=vb-va+1
  allocate (aw(va:vb), ap(va:vb), ae(va:vb), aq(va:vb), av(va-1:vb) )
  aw=zero; ap=one; ae=zero; aq=zero; av=zero

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend

     do i = i_mysta - 1, i_myend
        gx   = csi(1,i+1,j,k) + csi(1,i,j,k)
        gy   = csi(2,i+1,j,k) + csi(2,i,j,k)
        gz   = csi(3,i+1,j,k) + csi(3,i,j,k)

        ev  = (rei + pt5 * (q(5,i+1,j,k) + q(5,i,j,k))) / sig

        ucna = pt5 * abs(ucn_j(1,i+1,j,k) + ucn_j(1,i,j,k))

        av(i)= ucna + ev * (gx*gx + gy*gy + gz*gz) / (aj(i+1,j,k) + aj(i,j,k))

     end do

     do i = i_mysta, i_myend
        dtj  = pt5 * dtev(i,j,k) * aj(i,j,k) / dk(i,j,k)

        aw(i) =     - dtj * (ucn_j(1,i-1,j,k) + av(i-1)        )
        ap(i) = one + dtj * (                 + av(i-1) + av(i))
        ae(i) =       dtj * (ucn_j(1,i+1,j,k)           - av(i))

        ! transpose rhs
        aq(i) = rh(i,j,k) / dk(i,j,k)

     end do

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ka(1,nb) <= k .and. k <= li_blk_kb(1,nb) .and. &
            li_blk_ja(1,nb) <= j .and. j <= li_blk_jb(1,nb) ) then
            aw(li_blk_ia(1,nb):li_blk_ib(1,nb)) = zero
            ap(li_blk_ia(1,nb):li_blk_ib(1,nb)) = zero
            ae(li_blk_ia(1,nb):li_blk_ib(1,nb)) = zero
            if (li_blk_ia(1,nb) > i_mysta) ae(li_blk_ia(1,nb)-1) = zero
            if (li_blk_ib(1,nb) < i_myend) aw(li_blk_ib(1,nb)+1) = zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !

     call sgtsv(ilength, 1, aw(va+1:vb), ap(va:vb), ae(va:vb-1), aq(va:vb), &
                ilength, info)

     ! transpose rhs back
     do i = i_mysta, i_myend
        rh(i,j,k) = aq(i)
     end do

  end do
  end do

  deallocate (aw, ap, ae, aq, av)

  ! eta operator
  !
  va=j_mysta
  vb=j_myend
  jlength=vb-va+1
  allocate (aw(va:vb), ap(va:vb), ae(va:vb), aq(va:vb), av(va-1:vb) )
  aw=zero; ap=one; ae=zero; aq=zero; av=zero

  do k = k_mysta, k_myend
  do i = i_mysta, i_myend

     do j = j_mysta - 1, j_myend

        gx = eta(1,i,j+1,k) + eta(1,i,j,k)
        gy = eta(2,i,j+1,k) + eta(2,i,j,k)
        gz = eta(3,i,j+1,k) + eta(3,i,j,k)

        ev = (rei + pt5 * (q(5,i,j+1,k) + q(5,i,j,k))) / sig

        ucna = pt5 * abs(ucn_j(2,i,j+1,k) + ucn_j(2,i,j,k))

        av(j)= ucna + ev * (gx*gx + gy*gy + gz*gz) / (aj(i,j+1,k) + aj(i,j,k))
 
     end do

     do j = j_mysta, j_myend

        dtj  = pt5 * dtev(i,j,k) * aj(i,j,k) / dk(i,j,k)

        aw(j) =     - dtj * (ucn_j(2,i,j-1,k) + av(j-1)        )
        ap(j) = one + dtj * (                 + av(j-1) + av(j))
        ae(j) =       dtj * (ucn_j(2,i,j+1,k)           - av(j))

        ! transpose rhs
        aq(j) = rh(i,j,k)

     end do

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ka(1,nb) <= k .and. k <= li_blk_kb(1,nb) .and. &
            li_blk_ia(1,nb) <= i .and. i <= li_blk_ib(1,nb) ) then
            aw(li_blk_ja(1,nb):li_blk_jb(1,nb)) = zero
            ap(li_blk_ja(1,nb):li_blk_jb(1,nb)) = zero
            ae(li_blk_ja(1,nb):li_blk_jb(1,nb)) = zero
            if (li_blk_ja(1,nb) > j_mysta) ae(li_blk_ja(1,nb)-1) = zero
            if (li_blk_jb(1,nb) < j_myend) aw(li_blk_jb(1,nb)+1) = zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !
     call sgtsv(jlength, 1, aw(va+1:vb), ap(va:vb), ae(va:vb-1), aq(va:vb), &
                jlength, info)

     ! transpose rhs back
     !
     do j = j_mysta, j_myend
        rh(i,j,k) = aq(j)
     end do

  end do
  end do

  deallocate (aw, ap, ae, aq, av)

  ! zet operator
  !
  va=k_mysta
  vb=k_myend
  klength=vb-va+1
     
  allocate (aw(va:vb), ap(va:vb), ae(va:vb), aq(va:vb), av(va-1:vb) )
  aw=zero; ap=one; ae=zero; aq=zero; av=zero

  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     do k = k_mysta - 1, k_myend

        gx = zet(1,i,j,k+1) + zet(1,i,j,k)
        gy = zet(2,i,j,k+1) + zet(2,i,j,k)
        gz = zet(3,i,j,k+1) + zet(3,i,j,k)

        ev = (rei + pt5 * (q(5,i,j,k+1) + q(5,i,j,k))) / sig

        ucna = pt5 * abs(ucn_j(3,i,j,k+1) + ucn_j(3,i,j,k))

        av(k)= ucna + ev * (gx*gx + gy*gy + gz*gz) / (aj(i,j,k+1) + aj(i,j,k))
 
     end do

     do k = k_mysta, k_myend

        dtj  = pt5 * dtev(i,j,k) * aj(i,j,k) / dk(i,j,k)

        aw(k) =     - dtj * (ucn_j(3,i,j,k-1) + av(k-1)        )
        ap(k) = one + dtj * (                 + av(k-1) + av(k))
        ae(k) =       dtj * (ucn_j(3,i,j,k+1)           - av(k))

        ! transpose rhs
        aq(k) = rh(i,j,k)

     end do

     ! blanking area
     !
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ja(1,nb) <= j .and. j <= li_blk_jb(1,nb) .and. &
            li_blk_ia(1,nb) <= i .and. i <= li_blk_ib(1,nb) ) then
            aw(li_blk_ka(1,nb):li_blk_kb(1,nb)) = zero
            ap(li_blk_ka(1,nb):li_blk_kb(1,nb)) = zero
            ae(li_blk_ka(1,nb):li_blk_kb(1,nb)) = zero
            if (li_blk_ka(1,nb) > k_mysta) ae(li_blk_ka(1,nb)-1) = zero
            if (li_blk_kb(1,nb) < k_myend) aw(li_blk_kb(1,nb)+1) = zero
        end if
     end do
     end if

     ! solve linear tridiagonal eqns
     !
     call sgtsv(klength, 1, aw(va+1:vb), ap(va:vb), ae(va:vb-1), aq(va:vb), &
                klength, info)

     ! transpose rhs back
     !
     do k = k_mysta, k_myend
        rh(i,j,k) = aq(k)
     end do

  end do
  end do

  deallocate (aw, ap, ae, aq, av)

end subroutine rhs_sa_adi_solver

