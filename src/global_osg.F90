
module global_osg
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! global variables for overset grid
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  use precision
  implicit none

  ! multiblock interpolation method
  ! =0, simple interpolation; =1, mass flux based interpolation method
  integer :: mbim
  real (kind=rdf) :: coefv

  real (kind=rdf), dimension(:,:), allocatable :: qtemp

  integer, dimension(:), allocatable :: n_int

  ! nodes(a,b): a=i,j,k dir, b=b'th interface grid node, c=nzone
  ! hosts(a,b): a=i,j,k,nz,  b=
  !  coef(a,d): a=i,j,k, b=b'th grid node
  !
  integer, dimension(:,:,:), allocatable :: nodes
  integer, dimension(:,:,:), allocatable :: hosts
  real (kind=rdf), dimension(:,:,:), allocatable :: coef

  ! mpi version
  integer, dimension(:), allocatable :: donor
  integer, dimension(:), allocatable :: donee


end module global_osg


