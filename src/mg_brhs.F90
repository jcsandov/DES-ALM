
!============================================================
  subroutine mg_brhs (jl, ju, kl, ku, &
                    dc, de, dz, igp, jgp, kgp, &
                    eta, &
                    zet, &
                    ucn_j, &
                    aj, &
                    dtau, &
                    xnut, &
                    q, &
                    qn, &
                    qnm1, &
                    rh)

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! General curvilinear coordinates

  ! Calculate the rhs of the equations for use in bcond_out

  ! This routine is restricted to the finest grid only.
  !    i.e., no pk forcing (coarser grids only)
  !    i.e., always quick (2nd order dissipation)

  ! This routine is restricted to the i=1 or i=im boundary only.

  ! output

  ! rh(1:4) on (i=1 or i=im) boundary only
!============================================================

  use global_param
  use global_mpi

  implicit none

  ! extents
  ! -------
  integer :: il, jl, kl
  integer :: iu, ju, ku

  ! ghost points
  ! ------------
  integer :: igp, jgp, kgp

  ! main variables
  ! --------------
  real (kind = rdf), dimension(1:4,jl:ju,kl:ku), intent(in) :: q, qn, qnm1
  real (kind = rdf), dimension(1:4,jl:ju,kl:ku), intent(inout) :: rh
  real (kind = rdf), dimension(1:3,jl:ju,kl:ku), intent(in) :: eta, zet, ucn_j
  real (kind = rdf), dimension(jl:ju,kl:ku), intent(in) :: aj, dtau, xnut

  real (kind = rdf), dimension(:,:,:), allocatable :: visc
  real (kind = rdf), dimension(:,:,:), allocatable :: diss

  integer :: i_mysta, j_mysta, k_mysta
  integer :: i_myend, j_myend, k_myend

  real (kind = rdf) :: dc, de, dz

  integer :: i, j, k

  allocate (visc(1:3,jl:ju,kl:ku) , &
            diss(1:4,jl:ju,kl:ku) )

  ! process boundaries
  ! ------------------
  i_mysta=il+igp
  j_mysta=jl+jgp
  k_mysta=kl+kgp

  i_myend=iu-igp
  j_myend=ju-jgp
  k_myend=ku-kgp

  ! processes on the domain boundaries
  ! ----------------------------------
  if (myback == mpi_proc_null)  i_mysta = il + igp + 1
  if (myleft == mpi_proc_null)  j_mysta = jl + jgp + 1
  if (mydown == mpi_proc_null)  k_mysta = kl + kgp + 1

  if (myfront == mpi_proc_null) i_myend = iu - igp - 1
  if (myright == mpi_proc_null) j_myend = ju - jgp - 1
  if (myup    == mpi_proc_null) k_myend = ku - kgp - 1

  call brhs_viscous

  call brhs_diss_p

  call brhs_flux_sans_convec

  call brhs_convec_quick(2)

  call brhs_unst_visc_diss

  deallocate (visc, diss)

contains
  
  include 'brhs_viscous.F90'
  include 'brhs_diss_p.F90'
  include 'brhs_flux_sans_convec.F90'
  include 'brhs_convec_quick.F90'
  include 'brhs_unst_visc_diss.F90'

end subroutine mg_brhs


