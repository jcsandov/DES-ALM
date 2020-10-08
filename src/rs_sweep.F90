
subroutine rs_sweep ()

  !============================================================
  ! perform residual smoothing (with constant coefficients)
  ! in all directions
  ! 
  !============================================================


  implicit none

  integer :: info

  ! dimension lhs and temporary rhs
  !
  real (kind = rdf), dimension(:), allocatable   :: ax
  real (kind = rdf), dimension(:), allocatable   :: cx
  real (kind = rdf), dimension(:), allocatable   :: bx
  real (kind = rdf), dimension(:,:), allocatable :: dx

  real (kind = rdf), dimension(:), allocatable   :: ay
  real (kind = rdf), dimension(:), allocatable   :: cy
  real (kind = rdf), dimension(:), allocatable   :: by
  real (kind = rdf), dimension(:,:), allocatable :: dy

  real (kind = rdf), dimension(:), allocatable   :: az
  real (kind = rdf), dimension(:), allocatable   :: cz
  real (kind = rdf), dimension(:), allocatable   :: bz
  real (kind = rdf), dimension(:,:), allocatable :: dz

  integer :: i, j, k

  integer :: va, vb

  !
  allocate (ax(il:iu), & 
            cx(il:iu), & 
            bx(il:iu), & 
            dx(il:iu,1:4))
  allocate (ay(jl:ju), & 
            cy(jl:ju), & 
            by(jl:ju), & 
            dy(jl:ju,1:4))
  allocate (az(kl:ku), & 
            cz(kl:ku), & 
            bz(kl:ku), & 
            dz(kl:ku,1:4))

  ax=zero; bx=one; cx=zero; dx=zero
  ay=zero; by=one; cy=zero; dy=zero
  az=zero; bz=one; cz=zero; dz=zero

  ! csi sweep
  !
  va = i_mysta
  vb = i_myend

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend

     do i = i_mysta, i_myend
        dx(i,1:4) = rh(1:4,i,j,k)
     end do
     
     ax(:) = -epn(1)
     bx(:) = one + two * epn(1)
     cx(:) = -epn(1)

     ! blanking area
     if (nblk /= 0) then
     do nb = 1, nblk
        if (li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
            li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) ) then
            do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
               ax(i) = zero
               bx(i) = one
               cx(i) = zero
            end do
            if (li_blk_ia(n,nb) > i_mysta) cx(li_blk_ia(n,nb)-1) = zero
            if (li_blk_ib(n,nb) < i_myend) ax(li_blk_ib(n,nb)+1) = zero
        end if
     end do
     end if
     
     call sgtsv(ilength,4,ax(va+1:vb),bx(va:vb),cx(va:vb-1), dx(va:vb,1:4), &
                ilength, info)

     do i = i_mysta, i_myend
        rh(1:4,i,j,k) = dx(i,1:4)
     end do
  end do
  end do

  ! eta sweep
  !
  va = j_mysta
  vb = j_myend

  do k = k_mysta, k_myend
  do i = i_mysta, i_myend

     do j = j_mysta, j_myend
        dy(j,1:4) = rh(1:4,i,j,k)
     end do
     
     ay(:) = -epn(2)
     by(:) = one + two * epn(2)
     cy(:) = -epn(2)
     
     ! blanking area
     if (nblk /= 0) then
     do nb = 1, nblk
        if(li_blk_ka(n,nb) <= k .and. k <= li_blk_kb(n,nb) .and. &
           li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
           do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
              ay(j) = zero
              by(j) = one
              cy(j) = zero
           end do
           if (li_blk_ja(n,nb) > j_mysta) cy(li_blk_ja(n,nb)-1) = zero
           if (li_blk_jb(n,nb) > j_myend) ay(li_blk_jb(n,nb)+1) = zero
        end if
     end do
     end if

     call sgtsv(jlength,4,ay(va+1:vb),by(va:vb),cy(va:vb-1),dy(va:vb,1:4), &
                jlength, info)

     do j = j_mysta, j_myend
        rh(1:4,i,j,k) = dy(j,1:4)
     end do
  end do
  end do

  ! zet sweep
  !
  va = k_mysta
  vb = k_myend

  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     do k = k_mysta, k_myend
        dz(k,1:4) = rh(1:4,i,j,k)
     end do
     
     az(:) = -epn(3)
     bz(:) = one + two * epn(3)
     cz(:) = -epn(3)

     ! blanking area
     if (nblk /= 0) then
     do nb = 1, nblk
        if(li_blk_ja(n,nb) <= j .and. j <= li_blk_jb(n,nb) .and. &
           li_blk_ia(n,nb) <= i .and. i <= li_blk_ib(n,nb) ) then
           do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
              az(k) = zero
              bz(k) = one
              cz(k) = zero
           end do
           if (li_blk_ka(n,nb) > k_mysta) cz(li_blk_ka(n,nb)-1) = zero
           if (li_blk_kb(n,nb) < k_myend) az(li_blk_kb(n,nb)+1) = zero
        end if
     end do
     end if

     call sgtsv(klength,4,az(va+1:vb),bz(va:vb),cz(va:vb-1),dz(va:vb,1:4), &
                klength, info)

     do k = k_mysta, k_myend
        rh(1:4,i,j,k) = dz(k,1:4)
     end do
  end do
  end do

  !
  deallocate (ax, bx, cx, dx)
  deallocate (ay, by, cy, dy)
  deallocate (az, bz, cz, dz)

end subroutine rs_sweep

