
subroutine rhs_diss_p
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate dissipation for continutiy equation from fourth-order
  ! derivative of the pressure in all three directions

  ! input
  !     pdiss_coef (from input file) (global_param)
  !     csi(3,ijk), eta(3,ijk), zet(3,ijk)
  !     aj(ijk)
  !     q(1,ijk)
  !     dtau(ijk)

  ! output
  !     diss(1,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global_app
  implicit none

  ! local variables
  real (kind = rdf) :: tmp1, tmp2
  real (kind = rdf) :: dcsq, desq, dzsq
  real (kind = rdf) :: g1, g2, g3

  integer :: &
       ista, iend, &
       jsta, jend, &
       ksta, kend

  fv = zero
  diss = zero

  dcsq = dc * dc
  desq = de * de
  dzsq = dz * dz

  ! inner second derivatives and scaling

  ! csi-derivatives
  ! this scheme requires two ghost points on each side
  !
!   do k = k_mysta - 1, k_myend + 1
!   do j = j_mysta - 1, j_myend + 1
!   do i = i_mysta - 1, i_myend + 1

  ! use ista, iend to keep from having
  ! a array-bound access error when process
  ! owns the boundary

  ista = i_mysta - 1
  iend = i_myend + 1
  if (myback == mpi_proc_null)  ista = i_mysta
  if (myfront == mpi_proc_null) iend = i_myend

  do k = k_mysta-1, k_myend+1
  do j = j_mysta-1, j_myend+1
  do i = ista, iend

     g1   = csi(1,i,j,k) * csi(1,i,j,k) + &
            csi(2,i,j,k) * csi(2,i,j,k) + &
            csi(3,i,j,k) * csi(3,i,j,k)
     tmp1 = g1 * dtau(i,j,k) / aj(i,j,k)
     tmp2 = dcsq * (q(1,i+1,j,k) - two * q(1,i,j,k) + q(1,i-1,j,k))
     fv(1,i,j,k) = tmp1 * tmp2
  end do
  end do
  end do

  jsta = j_mysta - 1
  jend = j_myend + 1
  if (myleft == mpi_proc_null)  jsta = j_mysta
  if (myright == mpi_proc_null) jend = j_myend


  ! eta-derivatives
!   do k = k_mysta - 1, k_myend + 1
!   do j = j_mysta - 1, j_myend + 1
!   do i = i_mysta - 1, i_myend + 1

  do k = k_mysta-1, k_myend+1
  do j = jsta, jend
  do i = i_mysta-1, i_myend+1
     g2   = eta(1,i,j,k) * eta(1,i,j,k) + &
            eta(2,i,j,k) * eta(2,i,j,k) + &
            eta(3,i,j,k) * eta(3,i,j,k)
     tmp1 = g2 * dtau(i,j,k) / aj(i,j,k)
     tmp2 = desq * (q(1,i,j+1,k) - two * q(1,i,j,k) + q(1,i,j-1,k)) 
     fv(2,i,j,k) = tmp1 * tmp2
  end do
  end do
  end do

  ! zeta-derivatives
!   do k = k_mysta - 1, k_myend + 1
!   do j = j_mysta - 1, j_myend + 1
!   do i = i_mysta - 1, i_myend + 1

  ksta = k_mysta - 1
  kend = k_myend + 1
  if (mydown == mpi_proc_null)  ksta = k_mysta
  if (myup   == mpi_proc_null)  kend = k_myend

  do k = ksta, kend
  do j = j_mysta-1, j_myend+1
  do i = i_mysta-1, i_myend+1
     g3   = zet(1,i,j,k) * zet(1,i,j,k) + &
            zet(2,i,j,k) * zet(2,i,j,k) + &
            zet(3,i,j,k) * zet(3,i,j,k)
     tmp1 = g3 * dtau(i,j,k) / aj(i,j,k)
     tmp2 = dzsq * (q(1,i,j,k+1) - two * q(1,i,j,k) + q(1,i,j,k-1))
     fv(3,i,j,k) = tmp1 * tmp2
  end do
  end do
  end do

  ! boundary values
  ! 

  ! Put valid axu for nodes on domain (not process) boundary
  ! 
  if (myback == mpi_proc_null)    fv(1,i_mysta-1,:,:) = fv(1,i_mysta,:,:)
  if (myfront == mpi_proc_null)   fv(1,i_myend+1,:,:) = fv(1,i_myend,:,:)
  if (myleft == mpi_proc_null)    fv(2,:,j_mysta-1,:) = fv(2,:,j_mysta,:)
  if (myright == mpi_proc_null)   fv(2,:,j_myend+1,:) = fv(2,:,j_myend,:)
  if (mydown   == mpi_proc_null)  fv(3,:,:,k_mysta-1) = fv(3,:,:,k_mysta)
  if (myup == mpi_proc_null)      fv(3,:,:,k_myend+1) = fv(3,:,:,k_myend)

  ! boundary treatment of blanking nodes
  if (nblk /= 0) then
  do nb = 1, nblk

     if (blktype(1,nb,myzone) == 0) then
        i = li_blk_ia(n,nb)
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
        fv(1,i,j,k) = fv(1,i-1,j,k)
     end do
     end do
     end if

     if (blktype(2,nb,myzone) == 0) then
        i = li_blk_ib(n,nb)
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
        fv(1,i,j,k) = fv(1,i+1,j,k)
     end do
     end do
     end if

     if (blktype(3,nb,myzone) == 0) then
        j = li_blk_ja(n,nb)
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        fv(2,i,j,k) = fv(2,i,j-1,k)
     end do
     end do
     end if

     if (blktype(4,nb,myzone) == 0) then
        j = li_blk_jb(n,nb)
     do k = li_blk_ka(n,nb), li_blk_kb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        fv(2,i,j,k) = fv(2,i,j+1,k)
     end do
     end do
     end if

     if (blktype(5,nb,myzone) == 0) then
        k = li_blk_ka(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        fv(3,i,j,k) = fv(3,i,j,k-1)
     end do
     end do
     end if

     if (blktype(6,nb,myzone) == 0) then
        k = li_blk_kb(n,nb)
     do j = li_blk_ja(n,nb), li_blk_jb(n,nb)
     do i = li_blk_ia(n,nb), li_blk_ib(n,nb)
        fv(3,i,j,k) = fv(3,i,j,k+1)
     end do
     end do
     end if
  end do
  end if

  ! outer second derivative, interior nodes only
  !

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend
     diss(1,i,j,k) = pdiss_coef * ( &
          (fv(1,i+1,j,k) - two * fv(1,i,j,k) + fv(1,i-1,j,k)) + &
          (fv(2,i,j+1,k) - two * fv(2,i,j,k) + fv(2,i,j-1,k)) + &
          (fv(3,i,j,k+1) - two * fv(3,i,j,k) + fv(3,i,j,k-1)))
  end do
  end do
  end do

end subroutine rhs_diss_p



