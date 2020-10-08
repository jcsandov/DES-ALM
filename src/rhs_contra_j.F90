
subroutine rhs_contra_j
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate the contravariant velocity (U) divided by
  ! the Jacobian (J) at all nodes. U/J appears (slightly) more
  ! often than U so it makes sense to use U/J instead of U

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !
  real (kind = rdf) :: tmp

  do k = kl, ku
  do j = jl, ju
  do i = il, iu
     tmp = one / aj(i,j,k)
     ucn_j(1,i,j,k) = tmp * (csi(1,i,j,k) * q(2,i,j,k) + &
                             csi(2,i,j,k) * q(3,i,j,k) + &
                             csi(3,i,j,k) * q(4,i,j,k))
     ucn_j(2,i,j,k) = tmp * (eta(1,i,j,k) * q(2,i,j,k) + &
                             eta(2,i,j,k) * q(3,i,j,k) + &
                             eta(3,i,j,k) * q(4,i,j,k))
     ucn_j(3,i,j,k) = tmp * (zet(1,i,j,k) * q(2,i,j,k) + &
                             zet(2,i,j,k) * q(3,i,j,k) + &
                             zet(3,i,j,k) * q(4,i,j,k))
  end do
  end do
  end do
  
end subroutine rhs_contra_j



