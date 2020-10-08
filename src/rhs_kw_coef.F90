
  subroutine rhs_kw_coef ()

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  implicit none

  ! local parameter
  !
  !real (kind = rdf), parameter :: six  = 6.00_rdf

  real (kind = rdf), parameter :: pt09 = 0.09_rdf
  real (kind = rdf), parameter :: pt15 = 0.15_rdf
  real (kind = rdf), parameter :: pt27 = 0.27_rdf
  real (kind = rdf), parameter :: twopt7 = 2.70_rdf
  real (kind = rdf), parameter :: pt5p9 = 0.5555_rdf
  real (kind = rdf), parameter :: bc1 = 1137.77_rdf
  real (kind = rdf), parameter :: bc2 = 4096.00_rdf

  real (kind = rdf) :: Rt, Rt4

!   do k = 2, km - 1
!   do j = 2, jm - 1
!   do i = 1, im

!   integer :: ista, &
!              iend

!   ista = i_mysta
!   iend = i_myend

!   if ( myback  == mpi_proc_null ) ista = i_mysta - 1
!   if ( myfront == mpi_proc_null ) iend = i_myend + 1 

!   do k = k_mysta, k_myend
!   do j = j_mysta, j_myend
!   do i = ista, iend

  ! include ghost points
  ! 
  do k = kl, ku
  do j = jl, ju
  do i = il, iu

     Rt  = ren * q(5,i,j,k) / q(6,i,j,k)
     Rt4 = Rt * Rt * Rt * Rt

     As(i,j,k) = ( pt15 + Rt ) / ( six + Rt )
     Ac(i,j,k) =(( pt27 + Rt ) / ( twopt7 + Rt )) * pt5p9 / As(i,j,k)
     Bs(i,j,k) =(( bc1 + Rt4 ) / (bc2 + Rt4)) * pt09

  end do
  end do
  end do

  end subroutine rhs_kw_coef


