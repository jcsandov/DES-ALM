
subroutine rhs_sa_implicit_source ()

  implicit none 

  integer :: &
       ista, iend, &
       jsta, jend, &
       ksta, kend
  

  ! implicit treatment of positive contribution of source term
  !
  !      do k = 1, km
  !      do j = 1, jm
  !      do i = 1, im

  !============================================================
  !
  ! Boundaries need to be included to provide dk & kw on
  ! in the kw_solver_adi routine
  !
  !============================================================
  !
  ista = i_mysta
  jsta = j_mysta
  ksta = k_mysta

  iend = i_myend
  jend = j_myend
  kend = k_myend

  if (myback == mpi_proc_null)  ista = i_mysta - 1
  if (myleft == mpi_proc_null)  jsta = j_mysta - 1
  if (mydown == mpi_proc_null)  ksta = k_mysta - 1

  if (myfront == mpi_proc_null) iend = i_myend + 1
  if (myright == mpi_proc_null) jend = j_myend + 1
  if (myup    == mpi_proc_null) kend = k_myend + 1

  do k = ksta, kend
  do j = jsta, jend
  do i = ista, iend
     
!      dk(i,j,k) = one + dtev(i,j,k) * cw1 * fw(i,j,k) * q(5,i,j,k) &
!                      / wd(i,j,k) / wd(i,j,k)
     dk(i,j,k) = one
     
     rh(i,j,k) = - dtev(i,j,k) * aj(i,j,k) * rh(i,j,k)

  end do
  end do
  end do

end subroutine rhs_sa_implicit_source
   
