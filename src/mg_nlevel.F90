
subroutine mg_nlevel (n)
  use global
  use global_param
  use global_mpi

  implicit none
  
  real (kind = rdf), parameter :: one_sixteenth = one / sixteen

  integer :: ic, ie, iz
  integer :: i, j, k, l
  integer :: ii, jj, kk, ll
  integer :: n, m

  ! inject q to coarser grid
  ! 

  do k = gi_ka(n), gi_kb(n)
  do j = gi_ja(n), gi_jb(n)
  do i = gi_ia(n), gi_ib(n)
     kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
     jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
     ii = 2**(icrs(myzone)) * i - icrs(myzone)

     ll = gi_2_le_idx(ii,jj,kk,n-1)
     l  = gi_2_le_idx(i,j,k,n)
     q(1,l) = q(1,ll)
     q(2,l) = q(2,ll)
     q(3,l) = q(3,ll)
     q(4,l) = q(4,ll)
  end do
  end do
  end do

  ! les or rans, transfer xnut or prod here
  !
  if ( turbulence ) then
     do k = gi_ka(n), gi_kb(n)
     do j = gi_ja(n), gi_jb(n)
     do i = gi_ia(n), gi_ib(n)
        kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
        jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
        ii = 2**(icrs(myzone)) * i - icrs(myzone)
        ll = gi_2_le_idx(ii,jj,kk,n-1)
        l = gi_2_le_idx(i,j,k,n)
        xnut(l) = xnut(ll)
     end do
     end do
     end do
  end if
     
  ! Non-linear terms
  ! 
  if ( nlinc ) then
     do k = gi_ka(n), gi_kb(n)
     do j = gi_ja(n), gi_jb(n)
     do i = gi_ia(n), gi_ib(n)
        kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
        jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
        ii = 2**(icrs(myzone)) * i - icrs(myzone)
        ll = gi_2_le_idx(ii,jj,kk,n-1)
        l = gi_2_le_idx(i,j,k,n)
        uij(1,l) = uij(1,ll)
        uij(2,l) = uij(2,ll)
        uij(3,l) = uij(3,ll)
        uij(4,l) = uij(4,ll)
        uij(5,l) = uij(5,ll)
        uij(6,l) = uij(6,ll)
     end do
     end do
     end do
  end if

  !print *, 'nlevel - 3', n, myid
  ! save this solution for computing the
  ! residual on this grid level during
  ! prolongation
  ! 
  do l = le_idx_a(n), le_idx_b(n)
     qold_mg(1,l) = q(1,l)
     qold_mg(2,l) = q(2,l)
     qold_mg(3,l) = q(3,l)
     qold_mg(4,l) = q(4,l)
  end do

  !print *, 'nlevel - 4', n, myid
  ! Compute transfer operator for
  ! the residual vector (pk) from the
  ! right hand side on the finer grid level
  !
  ic = 1
  ie = le_ix(n-1)
  iz = le_ix(n-1) * le_jx(n-1)

  !print *, 'nlevel - 5', n, myid

  do k = gi_ka(n), gi_kb(n)
  do j = gi_ja(n), gi_jb(n)
  do i = gi_ia(n), gi_ib(n)
     kk = 2**(kcrs(myzone)) * k - kcrs(myzone)
     jj = 2**(jcrs(myzone)) * j - jcrs(myzone)
     ii = 2**(icrs(myzone)) * i - icrs(myzone)

     ll = gi_2_le_idx(ii,jj,kk,n-1)
      l = gi_2_le_idx(i,j,k,n)

      pk(1:4,l) = &
                  (pt25 * rh(1:4,ll+ic+ie+iz) &
                 + pt5  * rh(1:4,ll+ic+ie   ) &
                 + pt25 * rh(1:4,ll+ic+ie-iz) &
                                    
                 + pt5  * rh(1:4,ll   +ie+iz) &
                 +        rh(1:4,ll   +ie   ) &
                 + pt5  * rh(1:4,ll   +ie-iz) &
                                    
                 + pt25 * rh(1:4,ll-ic+ie+iz) &
                 + pt5  * rh(1:4,ll-ic+ie   ) &
                 + pt25 * rh(1:4,ll-ic+ie-iz) &
                                    
                 + pt5  * rh(1:4,ll+ic   +iz) &
                 +        rh(1:4,ll+ic      ) &
                 + pt5  * rh(1:4,ll+ic   -iz) &
                                    
                 +        rh(1:4,ll      +iz) &
                 + two  * rh(1:4,ll         ) &
                 +        rh(1:4,ll      -iz) &
                                    
                 + pt5  * rh(1:4,ll-ic   +iz) &
                 +        rh(1:4,ll-ic      ) &
                 + pt5  * rh(1:4,ll-ic   -iz) &
                                    
                 + pt25 * rh(1:4,ll+ic-ie+iz) &
                 + pt5  * rh(1:4,ll+ic-ie   ) &
                 + pt25 * rh(1:4,ll+ic-ie-iz) &
                                    
                 + pt5  * rh(1:4,ll   -ie+iz) &
                 +        rh(1:4,ll   -ie   ) &
                 + pt5  * rh(1:4,ll   -ie-iz) &
                                    
                 + pt25 * rh(1:4,ll-ic-ie+iz) &
                 + pt5  * rh(1:4,ll-ic-ie   ) &
                 + pt25 * rh(1:4,ll-ic-ie-iz)) * one_sixteenth
  end do
  end do
  end do

  !print *, 'nlevel - 6', n, myid

end subroutine mg_nlevel

! old stuff with prolong at the boundaries

!!$      do k=2,kmg(n)-1
!!$      do j=2,jmg(n)-1
!!$         kk=itwo**(kcrs) * k - kcrs
!!$         jj=itwo**(jcrs) * j - jcrs
!!$         l =ln(1,j,k,n)
!!$         ll=ln(1,jj,kk,n-1)

!!$          pk(1:4,l)=(     rh(1:4,ll-ie-iz) &
!!$                    + 2.0*rh(1:4,ll-iz   ) &
!!$                    +     rh(1:4,ll+ie-iz) &
!!$                    + 2.0*rh(1:4,ll-ie   ) &
!!$                    + 4.0*rh(1:4,ll      ) &
!!$                    + 2.0*rh(1:4,ll+ie   ) &
!!$                    +     rh(1:4,ll-ie+iz) &
!!$                    + 2.0*rh(1:4,ll+iz   ) &
!!$                    +     rh(1:4,ll+ie+iz)) *one_sixteenth

!!$      end do
!!$      end do

!!$      do k=2,kmg(n)-1
!!$      do j=2,jmg(n)-1
!!$         kk=itwo**(kcrs) * k - kcrs
!!$         jj=itwo**(jcrs) * j - jcrs
!!$         i=img(n)
!!$         ii=itwo*i-icrs
!!$         l =ln(i,j,k,n)
!!$         ll=ln(ii,jj,kk,n-1)

!!$          pk(1:4,l)=(     rh(1:4,ll-ie-iz) &
!!$                    + 2.0*rh(1:4,ll-iz   ) &
!!$                    +     rh(1:4,ll+ie-iz) &
!!$                    + 2.0*rh(1:4,ll-ie   ) &
!!$                    + 4.0*rh(1:4,ll      ) &
!!$                    + 2.0*rh(1:4,ll+ie   ) &
!!$                    +     rh(1:4,ll-ie+iz) &
!!$                    + 2.0*rh(1:4,ll+iz   ) &
!!$                    +     rh(1:4,ll+ie+iz)) * one_sixteenth

!!$      end do
!!$      end do

!!$ end subroutine new_nlevel









