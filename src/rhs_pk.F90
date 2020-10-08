
subroutine rhs_pk (decide_calc_pk)
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! For n > 1 this routine computes the forcing function Pk
  ! Pk is computed only for istage=1 and stored for subsequent use

  ! The logic should be transfered to the calling routine

  ! input
  !     decide_calc_pk
  !     pk(4,ijk) computed in nlevel
  !     rh(4,ijk) unmodified rhs for this grid level

  ! output
  !     rh(4,ijk) [modified by forcing function]

  ! Logic of fine grid forcing
  !-----------------------------------------------------------
  ! Grid Level         Istage    Result
  !-----------------------------------------------------------
  !   fine (n = 1)      n/a       do nothing
  !   coarse (n > 1)     1        compute pk & add pk to rh
  !   coarse (n > 1)     >1       add pk to rh
  !-----------------------------------------------------------
  ! Note: pk = pk (from nlevel) - rh, that is,
  !       the difference between average of rh on
  !       finer grid and rh calculated on this grid

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ! local decision variable
  !
  integer :: decide_calc_pk

  if (decide_calc_pk == 1) then
!      do k = 1, km
!      do j = 1, jm
!      do i = 1, im

     do k = k_mysta-1, k_myend+1
     do j = j_mysta-1, j_myend+1
     do i = i_mysta-1, i_myend+1
        pk(1,i,j,k) = pk(1,i,j,k) - rh(1,i,j,k)
        pk(2,i,j,k) = pk(2,i,j,k) - rh(2,i,j,k)
        pk(3,i,j,k) = pk(3,i,j,k) - rh(3,i,j,k)
        pk(4,i,j,k) = pk(4,i,j,k) - rh(4,i,j,k)
     end do
     end do
     end do

  end if

  ! add forcing function to the current residual
!     do k = 1, km
!     do j = 1, jm
!     do i = 1, im
       
     do k = k_mysta-1, k_myend+1
     do j = j_mysta-1, j_myend+1
     do i = i_mysta-1, i_myend+1
       rh(1,i,j,k) = rh(1,i,j,k) + pk(1,i,j,k)
       rh(2,i,j,k) = rh(2,i,j,k) + pk(2,i,j,k)
       rh(3,i,j,k) = rh(3,i,j,k) + pk(3,i,j,k)
       rh(4,i,j,k) = rh(4,i,j,k) + pk(4,i,j,k)
    end do
    end do
    end do

end subroutine rhs_pk



!!$  ! compute forcing terms Pk (=0 for n=1)
!!$  if (itr == 1 .and. istage == 1) then
!!$     if (n == ns) then
!!$        do m = 1, 4
!!$           do l = ls(n), le(n)
!!$              pk(l,m) = 0.0
!!$           end do
!!$        end do
!!$     else
!!$        do m = 1, 4
!!$           do l = ls(n), le(n)
!!$              aux(l,m) = pk(l,m) - rh(l,m)
!!$           end do
!!$        end do
!!$
!!$        do m = 1, 4
!!$           do l = ls(n), le(n)
!!$              pk(l,m) = aux(l,m)
!!$           end do
!!$        end do
!!$     end if
!!$  end if

!!$  ! add forcing term to the current residual
!!$  do m = 1, 4
!!$     do l = ls(n), le(n)
!!$        rh(l,m) = rh(l,m) + pk(l,m)
!!$     end do
!!$  end do

