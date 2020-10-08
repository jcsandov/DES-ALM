

! subroutine new_rsmooth(decide_grid_level, decide_smo_coef, &
!                        dc, de, dz, im, jm, km, dum_rh)


subroutine rsmooth ( n, decide_smo_coef, & ! 1
                     dc, de, dz, &      ! 2
                     il, iu, &          ! 3
                     jl, ju, &          ! 4
                     kl, ku, &          ! 5
                     igp, jgp, kgp, &   ! 6
                     rh )               ! 7
                     
  use global_param  
  use global_mpi
  use global_app
  use checksum
  
  implicit none

  integer :: &
       il, iu, &
       jl, ju, &
       kl, ku

  integer :: n

  integer :: igp, jgp, kgp
  
  integer :: i_mysta, &
             j_mysta, &
             k_mysta, &
             i_myend, &
             j_myend, &
             k_myend

  real (kind = rdf), dimension(1:4,il:iu,jl:ju,kl:ku) :: rh

  real (kind = rdf) :: dc, de, dz

  real (kind = rdf), dimension(3) :: epn

  integer :: i, j, k

  integer :: decide_smo_coef
  integer :: vector_len

  integer :: ilength
  integer :: jlength
  integer :: klength
  integer :: ks_klength

  !============================================================
  !
  !
  ! setup ij-direction sweeps
  !
  !
  !============================================================
  !
  ! Definitions below repeated in serveral places, e.g.
  ! solver_daf, kw_eddy, bcond_fm
  !
  i_mysta = il + igp
  j_mysta = jl + jgp
  k_mysta = kl + kgp

  i_myend = iu - igp
  j_myend = ju - jgp
  k_myend = ku - kgp

  ! processes on the domain boundaries
  ! 
  if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

  epn(:) = ep(decide_smo_coef,:,myzone)

  ilength = i_myend - i_mysta + 1
  jlength = j_myend - j_mysta + 1
  klength = k_myend - k_mysta + 1
  
  !============================================================
  !
  ! sweeps
  !
  !============================================================

  call rs_sweep ()

#ifdef DEBUG
  call init_cksum (igp, jgp, kgp)
  call cksum_4d_par ('rs.rh. ijk', rh)
#endif


contains
  
  include 'rs_sweep.F90'

end subroutine rsmooth


