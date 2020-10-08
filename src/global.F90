
module global
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! global variables for mainly for multigrid

  ! includes flow solution and grid transformation for all
  ! grids in the multigrid method

  ! Later, we can change this to allocate everything
  ! using dynamic memory allocation

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use precision
  implicit none

  ! grid parameters and variables
  integer, parameter :: ng    =  3

  ! variabes related to overset grid
  integer, parameter :: nzone   = 1 ! No. of zone
  integer, parameter :: nba_max = 1 ! Max no. of blanking area

  ! number of equations
  integer, parameter :: me = 5 ! s-a turbulence equations
 !integer, parameter :: me = 6 ! k-* turbulence equations

  integer :: nmxg

  ! both of these should be in input files

  integer, dimension(ng,nzone) :: img, jmg, kmg

  ! multigrid quantitites
  integer :: ns, nstop
  integer, dimension (ng) :: iter, itm

  integer, dimension (:,:,:,:,:), allocatable :: ln
  integer, dimension (:,:),       allocatable :: ls, le, icg, ieg, izg

  ! grid quantities
  !
  integer, dimension (nzone) :: icrs, jcrs, kcrs

  ! grid
  !
  real (kind = rdf), dimension(:), allocatable :: x, y, z

  ! metrics
  !
  real (kind = rdf), dimension(:,:), allocatable :: csi, eta, zet
  real (kind = rdf), dimension(:), allocatable :: aj

  ! now global variables:
  real (kind = rdf), dimension(:), allocatable :: xc, xe, xz, yc, ye, yz, zc, ze ,zz

  ! local time step
  !
  ! real (kind = rdf), dimension(:), allocatable :: dtau

  ! solution vectors
  !
  real (kind = rdf), dimension(:,:), allocatable :: q, qn, qnm1
  real (kind = rdf), dimension(:,:), allocatable :: qold_mg

  ! right-hand side quantities
  !
  real (kind = rdf), dimension(:,:), allocatable :: rh, pk
  
  ! vector for reynolds number
  ! 
  ! real (kind = rdf), dimension(:), allocatable :: ret

  ! daf, model matrices of Jacobian matrices A, B, C
  !
  ! real (kind = rdf), dimension(:,:,:), allocatable :: mai, n1i, n2i, mc

  ! daf, spectral radius
  !
  ! real (kind = rdf), dimension(:,:), allocatable :: spr

  ! les & des, right-hand-side information
  !
  real (kind = rdf), dimension(:), allocatable :: xnut
  real (kind = rdf), dimension(:), allocatable :: wd, dtev

  ! nonlinear turbulence closure
  ! 
  real (kind = rdf), dimension(:,:), allocatable :: uij

  ! mpi data structures
  ! 

  ! l = number according to local processor
  ! g = number according to global domain
  ! i = interior i.e., incl. ghost points
  ! e = exterior i.e., incl. ghost points
  ! a = start
  ! b = end
  ! x = number of points ---> b - a + 1
  !
  integer, dimension (:), allocatable, save :: li_ia, li_ja, li_ka
  integer, dimension (:), allocatable, save :: le_ia, le_ja, le_ka  
  integer, dimension (:), allocatable, save :: gi_ia, gi_ja, gi_ka
  integer, dimension (:), allocatable, save :: ge_ia, ge_ja, ge_ka
  
  integer, dimension (:), allocatable, save :: li_ib, li_jb, li_kb
  integer, dimension (:), allocatable, save :: le_ib, le_jb, le_kb  
  integer, dimension (:), allocatable, save :: gi_ib, gi_jb, gi_kb
  integer, dimension (:), allocatable, save :: ge_ib, ge_jb, ge_kb

  integer, dimension (:), allocatable, save :: li_ix, li_jx, li_kx
  integer, dimension (:), allocatable, save :: le_ix, le_jx, le_kx
  integer, dimension (:), allocatable, save :: gi_ix, gi_jx, gi_kx
  integer, dimension (:), allocatable, save :: ge_ix, ge_jx, ge_kx

  ! number of ghost points on each grid level
  ! 
  integer, dimension(:), allocatable, save :: igp
  integer, dimension(:), allocatable, save :: jgp
  integer, dimension(:), allocatable, save :: kgp

  ! idx = index for data structure
  ! 
  integer, dimension (:,:,:,:), allocatable, save :: li_idx
  integer, dimension (:,:,:,:), allocatable, save :: le_idx
  integer, dimension (:,:,:,:), allocatable, save :: gi_idx
  integer, dimension (:,:,:,:), allocatable, save :: gi_2_le_idx
 !integer, dimension (:,:,:,:), allocatable, save :: ge_2_le_idx

  ! a = beginning
  ! b = end
  ! 
  integer, dimension (:), allocatable, save :: li_idx_a
  integer, dimension (:), allocatable, save :: li_idx_b
  integer, dimension (:), allocatable, save :: le_idx_a
  integer, dimension (:), allocatable, save :: le_idx_b
  integer, dimension (:), allocatable, save :: gi_idx_a
  integer, dimension (:), allocatable, save :: gi_idx_b
  
  ! mx = max
  ! 
  integer, dimension (:), allocatable, save :: li_idx_mx
  integer, dimension (:), allocatable, save :: le_idx_mx
  integer, dimension (:), allocatable, save :: gi_idx_mx
  integer, dimension (:), allocatable, save :: ge_idx_mx

  ! convective scheme
  !
  character (len = 25), dimension (:), allocatable :: convec

  ! multigrid spacings
  ! 
  real (kind = rdf), dimension(:), allocatable :: dc
  real (kind = rdf), dimension(:), allocatable :: de
  real (kind = rdf), dimension(:), allocatable :: dz

  
end module global



