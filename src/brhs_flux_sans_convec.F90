
subroutine brhs_flux_sans_convec
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Orthogonal, curvlinear, cartesian coordinates

  ! Calculate pressure gradient and convective terms
  ! using second-order central differences.

  ! Calculate divergence of the velocity for the continuity
  ! equation also using second-order central differences.

  ! Unlike previous codes we will add viscous, dissipation, and
  ! unsteady terms in the calling routine not in this one.

  ! input
  !     ucn_j(3,ijk)
  !     csi(3,ijk)
  !     aj(ijk)
  !     q(4,ijk)
  !     dc, de, dz

  ! output
  !     rh(4,ijk)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  implicit none

  ! local arrays
  real (kind = rdf), dimension(:,:), allocatable :: p

  real (kind = rdf) :: de2, dz2
  
  ! local dummy reals
  real (kind = rdf) :: dudx, dvdy, dwdz

  ! dum --> dc2 or de2 or dz2
  de2 = pt5 * de
  dz2 = pt5 * dz
  
  !
  allocate (p(jl:ju,kl:ku))

  ! pressure gradient first
  do k = k_mysta - 1, k_myend + 1
  do j = j_mysta - 1, j_myend + 1
     p(j,k) = q(1,j,k) / aj(j,k)
  end do
  end do

  ! pressure gradient

  ! eta direction
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     rh(2,j,k) = de2 * (p(j+1,k) * eta(1,j+1,k) - p(j-1,k) * eta(1,j-1,k))
     rh(3,j,k) = de2 * (p(j+1,k) * eta(2,j+1,k) - p(j-1,k) * eta(2,j-1,k))
     rh(4,j,k) = de2 * (p(j+1,k) * eta(3,j+1,k) - p(j-1,k) * eta(3,j-1,k))
  end do
  end do
  
  ! zet direction
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     rh(2,j,k) = rh(2,j,k) + &
                 dz2 * (p(j,k+1) * zet(1,j,k+1) - p(j,k-1) * zet(1,j,k-1))
     rh(3,j,k) = rh(3,j,k) + &
                 dz2 * (p(j,k+1) * zet(2,j,k+1) - p(j,k-1) * zet(2,j,k-1))
     rh(4,j,k) = rh(4,j,k) + &
                 dz2 * (p(j,k+1) * zet(3,j,k+1) - p(j,k-1) * zet(3,j,k-1))
  end do
  end do
  
  ! velocity divergence
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     dudx = zero
     dvdy = de2 * (ucn_j(2,j+1,k) - ucn_j(2,j-1,k))
     dwdz = dz2 * (ucn_j(3,j,k+1) - ucn_j(3,j,k-1))

     rh(1,j,k) = dudx + dvdy + dwdz
  end do
  end do

  deallocate (p)

end subroutine brhs_flux_sans_convec


