
subroutine rhs_zero_bound_rh (decide_direction, decide_boundary)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Zero rh on all boundaries

  ! input
  !     decide_direction = 1, 2, 3 --> i, j, k
  !     decide_boundary = 1-->1, 2-->i or jm or km)
  !     rh(4,ijk)
  
  ! output
  !     rh(4,ijk)

  ! May want to unroll loops later, instead of using
  ! f90 array notation. I don't how this works with open_mp or apo

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!$  use precision
!!$  use global
  !
  use global_app
  implicit none

  ! decision variables
  !
  integer :: decide_direction
  integer :: decide_boundary

  if (decide_direction == 1) then
     if (decide_boundary == 1) then
        do k = ka, kb
        do j = ja, jb
           rh(:,ia-1,j,k) = zero
        end do
        end do
     else
        do k = ka, kb
        do j = ja, jb
           rh(:,ib+1,j,k) = zero
        end do
        end do
     end if
  else if (decide_direction == 2) then
     if (decide_boundary == 1) then
        do k = ka, kb
        do i = ia, ib
           rh(:,i,ja-1,k) = zero
        end do
        end do
     else
        do k = ka, kb
        do i = ia, ib
           rh(:,i,jb+1,k) = zero
        end do
        end do
     end if
  else if (decide_direction == 3) then
     if (decide_boundary == 1) then
        do j = ja, jb
        do i = ia, ib
           rh(:,i,j,ka-1) = zero
        end do
        end do
     else
        do j = ja, jb
        do i = ia, ib
           rh(:,i,j,kb+1) = zero
        end do
        end do
     end if
  else if (decide_direction ==4) then
        do k = kab1(decide_grid_level), kab2(decide_grid_level)
        do j = jab1(decide_grid_level), jab2(decide_grid_level)
        do i = iab1(decide_grid_level), iab2(decide_grid_level)
           rh(:,i,j,k) = zero
        end do
        end do
        end do
  end if

end subroutine rhs_zero_bound_rh




