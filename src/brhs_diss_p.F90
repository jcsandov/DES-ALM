
subroutine brhs_diss_p
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
  implicit none

  ! local arrays
  real (kind = rdf), dimension (:,:,:), allocatable :: aux

  ! local variables
  real (kind = rdf) :: tmp1, tmp2
  real (kind = rdf) :: desq, dzsq
  real (kind = rdf) :: g2, g3

  integer :: &
       jsta, jend, &
       ksta, kend

  allocate (aux(1:2,jl:ju,kl:ku))

  desq = de * de
  dzsq = dz * dz

  jsta = j_mysta - 1
  ksta = k_mysta - 1

  jend = j_myend + 1
  kend = k_myend + 1

  if (myleft == mpi_proc_null)  jsta = j_mysta
  if (mydown == mpi_proc_null)  ksta = k_mysta

  if (myright == mpi_proc_null) jend = j_myend
  if (myup    == mpi_proc_null) kend = k_myend

  ! inner second derivatives and scaling

  ! eta-derivatives
  do k = ksta, kend
  do j = jsta, jend
     g2 = eta(1,j,k)*eta(1,j,k) + eta(2,j,k)*eta(2,j,k) + eta(3,j,k)*eta(3,j,k)
     tmp1 = g2 * dtau(j,k) / aj(j,k)
     tmp2 = desq * (q(1,j+1,k) - two * q(1,j,k) + q(1,j-1,k)) 
     aux(1,j,k) = tmp1 * tmp2
  end do
  end do

  ! zeta-derivatives
  do k = ksta, kend
  do j = jsta, jend
     g3 = zet(1,j,k)*zet(1,j,k) + zet(2,j,k)*zet(2,j,k) + zet(3,j,k)*zet(3,j,k)
     tmp1 = g3 * dtau(j,k) / aj(j,k)
     tmp2 = dzsq * (q(1,j,k+1) - two * q(1,j,k) + q(1,j,k-1))
     aux(2,j,k) = tmp1 * tmp2
  end do
  end do

  ! boundary nodes for inner second derivatives
  if (myleft == mpi_proc_null)    aux(1,j_mysta-1,:) = aux(1,j_mysta,:)
  if (myright == mpi_proc_null)   aux(1,j_myend+1,:) = aux(1,j_myend,:)
  if (mydown   == mpi_proc_null)  aux(2,:,k_mysta-1) = aux(2,:,k_mysta)
  if (myup == mpi_proc_null)      aux(2,:,k_myend+1) = aux(2,:,k_myend)

  ! outer second derivative, interior nodes only
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
     diss(1,j,k) = pdiss_coef * ( &
          (aux(1,j+1,k) - two * aux(1,j,k) + aux(1,j-1,k)) + &
          (aux(2,j,k+1) - two * aux(2,j,k) + aux(2,j,k-1)))
  end do
  end do

  deallocate (aux)

end subroutine brhs_diss_p



