subroutine rhs_sa_unst_visc ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ! Calculate second-order accurate unsteady terms and add
  ! currunt viscous and source terms to rh

  !

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  implicit none

  ! local dummy variables
  real (kind = rdf) :: dtj
  real (kind = rdf) :: dq5dt

  ! Add time derivative, viscouse and source terms to rh
  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     dtj = two * delti * aj(i,j,k)

     dq5dt = e_source*(three*q(5,i,j,k) - four*qn(i,j,k) + qnm1(i,j,k)) / dtj

     rh(i,j,k) = rh(i,j,k) + dq5dt - visc(i,j,k) - psa(i,j,k) / aj(i,j,k)

  end do
  end do
  end do

end subroutine rhs_sa_unst_visc

