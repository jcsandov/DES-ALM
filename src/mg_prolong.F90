
subroutine mg_prolong (n)

  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  ! **** MPI PROCEDURE ****
  ! 
  ! Calculate dq include 1 ghost point
  !           (need to have included ghost points on qold)
  ! Interpolate dq to even nodes on current grid level (interior only)
  ! Calculate new q
  ! Exchange 1 ghost layer on this grid
  ! Repeat until finest grid
  ! 
  !    This routine prolongs the error vector from
  !    grid (n+1) --> grid n using linear interpolation
  !    along the grid lines

  ! input
  !     q(4,ls(n+1),le(n+1))
  !     qold(4,ls(n+1),le(n+1)) [used to be q0(4,:,2)]
  !                              still calculated in nlevel
  ! output
  !     q(4,ls(n),le(n))
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !

  implicit none 

  ! for now, just leave all in global
  !
 !real (kind = rdf) , dimension(4,le_idx_a(1):le_idx_b(ng)) :: dq
  real (kind = rdf) , dimension(:,:), allocatable :: dq

  real (kind = rdf) :: cp, cm

  integer :: k_mystart
  integer :: j_mystart
  integer :: i_mystart

  integer :: k_myend
  integer :: j_myend
  integer :: i_myend

  integer :: ii, jj, kk, ll, ic, ie, iz
  integer :: i, j, k, l
  integer :: n

  ! compute residual on the coarse grid
  ! and then inject it into the fine grid
  ! 
  ! initialize dq
  ! 
  allocate (dq(1:4,le_idx_a(1):le_idx_b(ng)))
  dq = zero

  do k = gi_ka(n+1), gi_kb(n+1)
  do j = gi_ja(n+1), gi_jb(n+1)
  do i = gi_ia(n+1), gi_ib(n+1)
     kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
     jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
     ii = 2**(icrs(myzone)) * i - icrs(myzone)

     ll = gi_2_le_idx(ii,jj,kk,n)
     l = gi_2_le_idx(i,j,k,n+1)

     dq(1,ll) = q(1,l) - qold_mg(1,l)
     dq(2,ll) = q(2,l) - qold_mg(2,l)
     dq(3,ll) = q(3,l) - qold_mg(3,l)
     dq(4,ll) = q(4,l) - qold_mg(4,l)
  end do
  end do
  end do

!  call mg_cksum ('pr- q', n+1, q)
!  call mg_cksum ('pr-old', n+1, qold_mg)
!  call mg_cksum ('pr-dq', n, dq)

  call mg_exchng3_2d (n, dq)

  ! fill in the blanks in the fine grid
  ! using linear interpolation
  !
  cp = pt5
  cm = pt5
   
  ic = 1
  ie = le_ix(n)
  iz = le_ix(n) * le_jx(n)

  if (kcrs(myzone) == 1) then

     k_mystart = gi_ka(n)
     j_mystart = gi_ja(n)
     i_mystart = gi_ia(n)

     k_myend = gi_kb(n)
     j_myend = gi_jb(n)
     i_myend = gi_ib(n)

     if (mod(gi_ka(n),2) /= 0) k_mystart = gi_ka(n) + 1
     if (mod(gi_kb(n),2) /= 0) k_myend = gi_kb(n) - 1

     do k = k_mystart, k_myend, 2**kcrs(myzone)
     do j = j_mystart, j_myend, 2**jcrs(myzone)
     do i = i_mystart, i_myend, 2**icrs(myzone)
        l = gi_2_le_idx(i,j,k,n)
        dq(1,l) = cp * dq(1,l+iz) + cm * dq(1,l-iz)
        dq(2,l) = cp * dq(2,l+iz) + cm * dq(2,l-iz)
        dq(3,l) = cp * dq(3,l+iz) + cm * dq(3,l-iz)
        dq(4,l) = cp * dq(4,l+iz) + cm * dq(4,l-iz)
     end do
     end do
     end do
  end if


  if (jcrs(myzone) == 1) then

     k_mystart = gi_ka(n)
     j_mystart = gi_ja(n)
     i_mystart = gi_ia(n)

     k_myend = gi_kb(n)
     j_myend = gi_jb(n)
     i_myend = gi_ib(n)
     
     if (mod(gi_ja(n),2) /= 0) j_mystart = gi_ja(n) + 1
     if (mod(gi_jb(n),2) /= 0) j_myend = gi_jb(n) - 1

     do k = k_mystart, k_myend
     do j = j_mystart, j_myend, 2**jcrs(myzone)
     do i = i_mystart, i_myend, 2**icrs(myzone)
        l = gi_2_le_idx(i,j,k,n)
        dq(1,l) = cp * dq(1,l+ie) + cm * dq(1,l-ie)
        dq(2,l) = cp * dq(2,l+ie) + cm * dq(2,l-ie)
        dq(3,l) = cp * dq(3,l+ie) + cm * dq(3,l-ie)
        dq(4,l) = cp * dq(4,l+ie) + cm * dq(4,l-ie)
     end do
     end do
     end do
  end if

  if (icrs(myzone) == 1) then

     k_mystart = gi_ka(n)
     j_mystart = gi_ja(n)
     i_mystart = gi_ia(n)

     k_myend = gi_kb(n)
     j_myend = gi_jb(n)
     i_myend = gi_ib(n)

     if (mod(gi_ia(n),2) /= 0) i_mystart = gi_ia(n) + 1
     if (mod(gi_ib(n),2) /= 0) i_myend = gi_ib(n) - 1

     do k = k_mystart, k_myend
     do j = j_mystart, j_myend
     do i = i_mystart, i_myend, 2**icrs(myzone)
        l = gi_2_le_idx(i,j,k,n)
        dq(1,l) = cp * dq(1,l+ic) + cm * dq(1,l-ic)
        dq(2,l) = cp * dq(2,l+ic) + cm * dq(2,l-ic)
        dq(3,l) = cp * dq(3,l+ic) + cm * dq(3,l-ic)
        dq(4,l) = cp * dq(4,l+ic) + cm * dq(4,l-ic)
     end do
     end do
     end do
  end if

  ! blanking area
  !
  if (nblk /= 0) then
  do nb = 1, nblk
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        l = li_idx(i,j,k,n)
        dq(1:4,l) = zero
     end do
     end do
     end do
  end do
  end if

  ! ============================================================
  !
  ! MPI
  !
  ! Don't include boundaries at all in smoothing
  ! only interior nodes
  !
  ! ============================================================
  !
  la = le_idx_a(n)
  lb = le_idx_b(n)

  call rsmooth ( n, 2, &                   ! 1
                 dc(n), de(n), dz(n), &    ! 2
                 le_ia(n), le_ib(n),  &    ! 3
                 le_ja(n), le_jb(n),  &    ! 4
                 le_ka(n), le_kb(n),  &    ! 5
                 igp(n), jgp(n), kgp(n), & ! 6
                 dq(1:4,la:lb) )           ! 7


!============================================================
!
!
!  MPI ---> sync bc on dq with paik's code
!
!============================================================
!
! boundary conditions for dq need to
! sync this with Joongcheol's code
! 
  if (myback == mpi_proc_null) then
     do k = gi_ka(n), gi_kb(n)
     do j = gi_ja(n), gi_jb(n)
        l = gi_2_le_idx(gi_ia(n),j,k,n)
        dq(1:4,l) = zero
     end do
     end do
  end if

  if (myfront == mpi_proc_null) then
     do k = gi_ka(n), gi_kb(n)
     do j = gi_ja(n), gi_jb(n)
        l = gi_2_le_idx(gi_ib(n),j,k,n)
        dq(1:4,l) = zero
     end do
     end do
  end if
  
  if (myleft == mpi_proc_null) then
     do k = gi_ka(n), gi_kb(n)
     do i = gi_ia(n), gi_ib(n)
        l = gi_2_le_idx(i,gi_ja(n),k,n)
        dq(1:4,l) = zero
     end do
     end do
  end if

  if (myright == mpi_proc_null) then
     do k = gi_ka(n), gi_kb(n)
     do i = gi_ia(n), gi_ib(n)
        l = gi_2_le_idx(i,gi_jb(n),k,n)
        dq(1:4,l) = zero
     end do
     end do
  end if

  if (mydown == mpi_proc_null) then
     do j = gi_ja(n), gi_jb(n)
     do i = gi_ia(n), gi_ib(n)
        l = gi_2_le_idx(i,j,gi_ka(n),n)
        dq(1:4,l) = zero
     end do
     end do
  end if

  if (myup == mpi_proc_null) then
     do j = gi_ja(n), gi_jb(n)
     do i = gi_ia(n), gi_ib(n)
        l = gi_2_le_idx(i,j,gi_kb(n),n)
        dq(1:4,l) = zero
     end do
     end do
  end if

  do l = le_idx_a(n), le_idx_b(n)
     q(1,l) = q(1,l) + dq(1,l)
     q(2,l) = q(2,l) + dq(2,l)
     q(3,l) = q(3,l) + dq(3,l)
     q(4,l) = q(4,l) + dq(4,l)
  end do

  deallocate(dq)

end subroutine mg_prolong

! boundary conditions on dq that were eliminated

!!$      do j=1,jm
!!$      do k=1,km
!!$      l=ln(1,j,k,n)
!!$      ll=ln(im,j,k,n)
!!$      dq(ll,1)=0.0
!!$      dq(l,2)=0.0
!!$      dq(l,3)=0.0
!!$      dq(l,4)=0.0
!!$      end do
!!$      end do


!calculation of residual and inject are done separately

! old prolong stuff, keep for reference

!!$      do l=ls(n+1),le(n+1)
!!$      dq(l,m) = q(l,m) - q0(l,m,2)
!!$      end do
!!$
!!$
!!$      do k=1,kmg(n),2
!!$           kk=(k+1)/2
!!$      do j=1,jmg(n),2
!!$           jj=(j+1)/2
!!$      do i=1,img(n),2
!!$           ii=(i+1)/2
!!$
!!$      l=ln(i,j,k,n)
!!$      ll=ln(ii,jj,kk,n+1)
!!$
!!$      dq(l,m) = dq(ll,m)
!!$      end do
!!$      end do
!!$      end do
!!$      cp=0.5
!!$      cm=0.5
!!$
!!$   do k=2,kmg(n)-1,2
!!$      do j=1,jmg(n)  ,2
!!$         do i=1,img(n)  ,2
!!$
!!$      l=ln(i,j,k,n)
!!$      dq(l,m)=cp*dq(l+iz,m)+cm*dq(l-iz,m)
!!$      end do
!!$      end do
!!$      end do
!!$
!!$      do k=1,kmg(n)
!!$      do j=1,jmg(n),2
!!$      do i=2,img(n)-1,2
!!$
!!$      l=ln(i,j,k,n)
!!$      dq(l,m)=cp*dq(l+ic,m)+cm*dq(l-ic,m)
!!$      end do
!!$      end do
!!$      end do
!!$
!!$      do k=1,kmg(n)
!!$      do j=2,jmg(n)-1,2
!!$      do i=1,img(n)
!!$
!!$      l=ln(i,j,k,n)
!!$      dq(l,m)=cp*dq(l+ie,m)+cm*dq(l-ie,m)
!!$      end do
!!$      end do
!!$      end do
!!$
!!$      end do

! ... Smooth the residual before correcting the fine grid

!!$      sdm=ep(2,1)+ep(2,2)+ep(2,3)
!!$      if(sdm.ne.0) call rsmo(n,i12,i13,i23,nij,nik,njk,dq)








