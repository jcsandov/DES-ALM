
subroutine mg_residual (n, erp, eru, erv, erw)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! This routine computes the L2-norm of the residual
  ! on the finest grid at the end of each pseudo time step

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  real (kind = rdf) :: erp, eru, erv, erw, erx
  real (kind = rdf) :: tnu

  real (kind = rdf), dimension (5) :: resd
  real (kind = rdf), dimension (5) :: resd_tmp

  integer :: i, j, k, l
  integer :: m
  integer :: n

  integer :: k_mysta
  integer :: j_mysta
  integer :: i_mysta

  integer :: k_myend
  integer :: j_myend
  integer :: i_myend

  integer, dimension (0:nproc-1) :: recvcounts

  tnu = real((img(n) - 2) * (jmg(n) - 2) * (kmg(n) - 2), kind = rdf)

  k_mysta = li_ka(n)
  j_mysta = li_ja(n)
  i_mysta = li_ia(n)

  k_myend = li_kb(n)
  j_myend = li_jb(n)
  I_myend = li_ib(n)

  ! use interior nodes only
  ! 
  if (myback  == mpi_proc_null) i_mysta = i_mysta + 1
  if (myleft  == mpi_proc_null) j_mysta = j_mysta + 1
  if (mydown  == mpi_proc_null) k_mysta = k_mysta + 1

  if (myfront == mpi_proc_null) i_myend = i_myend - 1
  if (myright == mpi_proc_null) j_myend = j_myend - 1
  if (myup    == mpi_proc_null) k_myend = k_myend - 1

  ! calculate L2 norm
  ! 
  resd_tmp = zero
  resd     = zero

  do k = k_mysta, k_myend
     do j = j_mysta, j_myend
        do i = i_mysta, i_myend
           l = le_idx(i,j,k,n)
           do m = 1, me
              resd_tmp(m) = resd_tmp(m) + (q(m,l) - qold_mg(m,l))**2
           end do
        end do
     end do
  end do

  resd_tmp(5) = pt5             ! kludge

  recvcounts(:) = 5

  call mpi_reduce (resd_tmp, resd, 5, mpi_real, mpi_sum,&
       & root, mpi_comm_world, ierr)

  !   call mpi_reduce_scatter (resd_tmp, resd, recvcounts, mpi_real, mpi_sum,&
  !        & mpi_comm_world, ierr)

  ! calculate average residual
  !
  if ( myid == root ) then

     resd(:) = sqrt(resd(:) / tnu)

     erp = log10(resd(1))
     eru = log10(resd(2))
     erv = log10(resd(3))
     erw = log10(resd(4))

     if (me == 5) erx = log10(resd(5))

  else
     erp = zero
     eru = zero
     erv = zero
     erw = zero

  end if

end subroutine mg_residual

