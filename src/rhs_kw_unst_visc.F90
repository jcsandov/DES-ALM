subroutine rhs_kw_unst_visc ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  ! Calculate second-order accurate unsteady terms and add
  ! currunt viscous and source terms to rh

  !

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  implicit none

  ! local dummy variables
  real (kind = rdf) :: dtj
  real (kind = rdf) :: dq5dt, dq6dt

  ! Add time derivative, viscouse and source terms to rh

!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 2, im - 1

  do k = k_mysta, k_myend
  do j = j_mysta, j_myend
  do i = i_mysta, i_myend

     dtj = two * delti * aj(i,j,k)

     dq5dt = e_source*(three*q(5,i,j,k)-four*qn(1,i,j,k)+qnm1(1,i,j,k)) /dtj
     dq6dt = e_source*(three*q(6,i,j,k)-four*qn(2,i,j,k)+qnm1(2,i,j,k)) /dtj

     rh(1,i,j,k) = rh(1,i,j,k) + dq5dt - visc(1,i,j,k) + pkw(1,i,j,k)
     rh(2,i,j,k) = rh(2,i,j,k) + dq6dt - visc(2,i,j,k) + pkw(2,i,j,k)

  end do
  end do
  end do

!   ! zero rh on blank nodes
!   if (abut) then
!   do k = 2, km - 1
!   do j = jab1(1), jab2(1)
!   do i = iab1(1), iab2(1)
!      rh(1,i,j,k) = zero
!      rh(2,i,j,k) = zero
!   end do
!   end do
!   end do
!   end if

end subroutine rhs_kw_unst_visc


