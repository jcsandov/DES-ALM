
subroutine mg_inject ( dum )

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! inject (transfer) a vector with 1 dof per node from the
  ! fine grid to all coarser grids
  !
  ! does not fix ghost layers
  ! 
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use global
  use global_mpi

  implicit none

  integer :: kk, jj, ii, ll
  integer :: i, j, k, l
  integer :: n

  real (kind = rdf), dimension(le_idx_a(1):le_idx_b(ng)), intent(inout) :: dum

  real (kind = rdf), dimension(size(dum)) :: aux

  ! transfer vector to all coarser grids, but only interior nodes
  !
  do n = 2, ng
     do k = gi_ka(n), gi_kb(n)
        do j = gi_ja(n), gi_jb(n)
           do i = gi_ia(n), gi_ib(n)
              kk = 2**(kcrs(myzone)) * k - kcrs(myzone) - gi_ka(n-1) + 1
              jj = 2**(jcrs(myzone)) * j - jcrs(myzone) - gi_ja(n-1) + 1
              ii = 2**(icrs(myzone)) * i - icrs(myzone) - gi_ia(n-1) + 1

              kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
              jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
              ii = 2**(icrs(myzone)) * i - icrs(myzone)

              ll = gi_2_le_idx(ii,jj,kk,n-1)
               l = gi_2_le_idx(i,j,k,n)


!               ll = le_idx(ii,jj,kk,n-1)
!               l = le_idx(i,j,k-gi_ka(n)+1,n)

              aux(l) = dum(ll)


           end do
        end do
     end do

     do k = gi_ka(n), gi_kb(n)
        do j = gi_ja(n), gi_jb(n)
           do i = gi_ia(n), gi_ib(n)
!               l = le_idx(i,j,k-gi_ka(n)+1,n)
              l = gi_2_le_idx(i,j,k,n)
              dum(l) = aux(l)
           end do
        end do
     end do

  end do

end subroutine mg_inject


